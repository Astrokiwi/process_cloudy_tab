module constants
    real(kind=8), parameter :: mean_molecular_mass = 1.22
    real(kind=8), parameter :: amu = 1.66054d-24 ! grams
    real(kind=8), parameter :: solar_mass = 1.9891e33 ! grams
    real(kind=8), parameter :: pc = 3.08567758e18 ! cm
    
    real(kind=8), parameter :: plancks_constant = 6.62607d-27 ! cm**2 g/s
    real(kind=8), parameter :: c_speed = 29979245800. ! cm/s
    real(kind=8), parameter :: boltzmann_constant = 1.38d-16 ! erg/K

    real(kind=8), parameter :: blackbody_constant_cgs = plancks_constant*c_speed/boltzmann_constant ! cm K
    real(kind=8), parameter :: blackbody_constant_micronK = blackbody_constant_cgs*1.d4 ! um K
    
    real(kind=8), parameter :: astro_density_to_cgs = solar_mass/pc**3 ! solar masses/pc**3 to g/cm**3
    real(kind=8), parameter :: astro_surf_density_to_cgs = solar_mass/pc**2 ! solar masses/pc**2 to g/cm**2

    real(kind=8),  parameter :: PI_8  = 4 * atan (1.0_8)


end module constants

module cross_section_calc
    integer :: n_wavelengths
    real(kind=8), allocatable, dimension(:) :: wavelength,opacity,blackbody
    
    real(kind=8) :: bb_norm
    
    contains
        subroutine load(opacity_file,input_density)
            use constants
            implicit none
            real(kind=8),intent(in) :: input_density
            character(len=*),intent(in) :: opacity_file
            
            integer :: ierr,iw
            real(kind=8) :: dummy
            real(kind=8) :: volume_opacity_to_mass_opacity
        
            n_wavelengths = 0
            ierr=0
            open(file=opacity_file,unit=15)
            read(15,*) ! title line
            do while (ierr==0)
                read(15,*,iostat=ierr)
                if ( ierr==0 ) n_wavelengths=n_wavelengths+1
            end do
            close(15)
        
            if ( allocated(wavelength) ) deallocate(wavelength)
            if ( allocated(opacity) ) deallocate(opacity)
            if ( allocated(blackbody) ) deallocate(blackbody)
        
            allocate(wavelength(n_wavelengths))
            allocate(opacity(n_wavelengths))
            allocate(blackbody(n_wavelengths))
        
            ! reverse as we read in to get monotonically increasing values
            open(file=opacity_file,unit=15)
            read(15,*) ! title line
            do iw=1,n_wavelengths
                read(15,*) dummy,wavelength(n_wavelengths-iw+1),opacity(n_wavelengths-iw+1)
            end do
            close(15)
        
            volume_opacity_to_mass_opacity = 1./(mean_molecular_mass*amu*input_density) ! cm**3/g
            opacity=opacity*volume_opacity_to_mass_opacity
        end subroutine load
        
        subroutine set_temp(T_K)
            use constants
            implicit none
            real(kind=8), intent(in) :: T_K
            integer :: iw
            blackbody = (wavelength**(-5))/(dexp(blackbody_constant_micronK/(T_K*wavelength))-1.d0)
            bb_norm = 0.
            do iw=2,n_wavelengths
                bb_norm = bb_norm + (blackbody(iw-1)+blackbody(iw))*(wavelength(iw)-wavelength(iw-1))/2.d0
            end do
!             print *,"bb_norm=",bb_norm
            
!             open(file="data/dump",unit=15)
!             do iw=1,n_wavelengths
!                 write(15,*) wavelength(iw),opacity(iw),blackbody(iw)
!             end do
!             close(15)
        end subroutine set_temp
        
        real(kind=16) function cross_section(m_Msun,h_pc,density_func,surf_func,nr,nz) result(cross)
            use constants
            implicit none

            real(kind=8), intent(in) :: m_Msun,h_pc
            integer,intent(in) :: nr,nz
            real(kind=8) :: density_func,surf_func
            
            integer :: ir,iz,iw
            real(kind=8) :: r,z,dv,dr,dz
            
            real(kind=8) :: surf_norm,dens_cgs
            real(kind=16) :: dcross
            
            integer :: threadid
            integer, external :: OMP_GET_THREAD_NUM
            
            surf_norm = 3./4./PI_8 * m_Msun/h_pc**2 * astro_surf_density_to_cgs
            dens_cgs = 3./4./PI_8 * m_Msun/h_pc**3 * astro_density_to_cgs
            dr = 1.d0/nr
            dz = 2.d0/nz
            dv=dr*dz
!             print *,surf_func(0.,0.),density_func(0.,0.)
!             print *,surf_norm,dens_cgs

            ! dunno how to do trapezium in 3D, let's just do rectangles for integral
!$OMP  PARALLEL DO private(ir,iz,iw,r,z,threadid,dcross)&
!$OMP  shared(nr,nz,n_wavelengths,dz,dr,dv,blackbody,opacity,wavelength,surf_norm)&
!$OMP  default(none)&
!$OMP  schedule(dynamic)&
!$OMP  reduction(+:cross)
            do iz=1,nz-1
!                 threadid=OMP_GET_THREAD_NUM()
!                 if ( threadid==0 ) then
!                     print *,iz,"/",nz,cross
!                 endif
                z = 1.-iz*dz
                do ir=1,nr-1
                    r = ir*dr
!                       cross=cross+                                    &
!                                     r*density_func(r,z)

                    do iw=1,n_wavelengths-1
                        dcross=                                    &
                                    blackbody(iw)*opacity(iw)*         &
                                    density_func(r,z)*                  &
                                    r*                                    &
                                    dexp(-opacity(iw)*surf_func(r,z)*surf_norm)*  &
                                    (wavelength(iw+1)-wavelength(iw))
                        if ( isNaN(dcross) ) then
                            print *,"adding nan!"
                            print *,blackbody(iw),opacity(iw),density_func(r,z),r,&
                            surf_func(r,z),dexp(-opacity(iw)*surf_func(r,z)*surf_norm),&
                            (wavelength(iw+1)-wavelength(iw))
                            stop
                        endif
                        cross = cross + dcross
                    end do
                end do
            end do
!$OMP  END PARALLEL DO

            cross=cross*dens_cgs*dv*2.d0*PI_8*(h_pc)**3*pc/bb_norm ! in pc**2

!             cross=cross*(dens_cgs/astro_density_to_cgs)*dv*2.d0*PI_8*h_pc**3 ! mass in Msun

            
!             cross=cross*dens_cgs/bb_norm*dv*2.d0*PI_8
!             cross=cross*dens_cgs*dv*2.d0*PI_8*(h_pc*pc)**3/bb_norm ! in cm**2
!             print *,cross,"cm**2",log10(cross)
!             print *,cross/pc**2,"pc**2",log10(cross/pc**2)
!             print *,cross/(PI_8*(h_pc*pc)**2)," fraction absorbed"
!             print *,surf_norm,dens_cgs

            return
        end function

end module cross_section_calc

module constant_density_functions
    contains
    
        real(kind=8) function constant_density(r,z)
            implicit none
            real(kind=8), intent(in) :: r,z

            if ( z**2+r**2>=1.d0 ) then
                constant_density = 0.
            else
                constant_density = 1.
            endif
            return

        end function constant_density

        real(kind=8) function constant_density_surf(r,z)
            implicit none
            real(kind=8), intent(in) :: r,z

            if ( z**2+r**2>=1.d0) then
                constant_density_surf = 0.
                return
            endif
            constant_density_surf = sqrt(1.d0-r**2)-z
            return
        end function constant_density_surf

end module constant_density_functions


module sph_density_functions
    real(kind=8), dimension(:,:), allocatable :: surf_table
    integer :: nz_sph,nr_sph
    real(kind=8) :: dr,dz

    contains
        subroutine setup_surf(nr_in,nz_in)
            use constants, only: PI_8
            integer, intent(in) :: nr_in,nz_in
            integer :: ir,iz
            real(kind=8) :: r,z
            
!             real(kind=8) :: test_val
            
            nz_sph=nz_in
            nr_sph=nr_in
            if ( allocated(surf_table) ) then
                deallocate(surf_table)
            endif
            allocate(surf_table(nr_in,nz_in))

            dr = 1.d0/nr_sph
            dz = 2.d0/nz_sph
            
            surf_table(:,1) = 0.d0
            
            do ir=1,nr_sph
                r = (ir-1)*dr
                do iz=2,nz_sph
                    z = 1.-iz*dz
                    surf_table(ir,iz)=surf_table(ir,iz-1)+(sph_density(r,z)+sph_density(r,z+dz))/2.d0*dz
                end do
            end do
            
!             test_val = 0.d0
!             do ir=1,nr_sph
!                 r=(ir-1)*dr
!                 test_val=test_val+sph_surf(r,0.d0)*r*PI_8*dr*2.
!             end do
!             print *,test_val,test_val*2./(4./3.*PI_8)
!             
!             stop
            
        end subroutine setup_surf

        real(kind=8) function sph_density(r,z)
            use constants, only: PI_8
            implicit none
            real(kind=8), intent(in) :: r,z
            real(kind=8) :: x
            
            x=sqrt(r**2+z**2)

            if ( x<0.5 ) then
                sph_density = 1.-6.*x**2+6.*x**3
            else if ( x<=1. ) then
                sph_density = 2.*(1.-x)**3
            else
                sph_density = 0.
            endif
            sph_density=sph_density*(32.d0/3.d0)


            return

        end function sph_density

        real(kind=8) function sph_surf(r,z)
            implicit none
            real(kind=8), intent(in) :: r,z
            integer :: ir,iz
            real(kind=8) :: lowr_lowz_val,lowr_highz_val,highr_lowz_val,highr_highz_val
            real(kind=8) :: r_weight,z_weight
            
            iz = (1.d0-z)*nz_sph/2
            if ( iz<1 ) then
                sph_surf=0.d0
                return
            endif
            if ( iz>nz_sph ) then
                iz=nz_sph
            endif
            ir = r*nr_sph+1
            if ( ir<1 .or. r>nr_sph ) then
                sph_surf=0.d0
                return
            endif
            
            z_weight = (z-(1.d0-iz*dz))/dz
            r_weight = r/dr-ir+1
            if ( z_weight<0. ) z_weight=0.
            if ( z_weight>1. ) z_weight=1.
            if ( r_weight<0. ) r_weight=0.
            if ( r_weight>1. ) r_weight=1.
            if ( ir==nr_sph ) then
                if ( iz==nz_sph ) then
                    lowr_lowz_val=surf_table(ir,iz)
                    lowr_highz_val=surf_table(ir,iz)
                    highr_lowz_val=surf_table(ir,iz)
                    highr_highz_val=surf_table(ir,iz)
                else
                    lowr_lowz_val=surf_table(ir,iz)
                    lowr_highz_val=surf_table(ir,iz+1)
                    highr_lowz_val=surf_table(ir,iz)
                    highr_highz_val=surf_table(ir,iz+1)
                endif
            else
                if ( iz==nz_sph ) then
                    lowr_lowz_val=surf_table(ir,iz)
                    lowr_highz_val=surf_table(ir,iz)
                    highr_lowz_val=surf_table(ir+1,iz)
                    highr_highz_val=surf_table(ir+1,iz)
                else
                    lowr_lowz_val=surf_table(ir,iz)
                    lowr_highz_val=surf_table(ir,iz+1)
                    highr_lowz_val=surf_table(ir+1,iz)
                    highr_highz_val=surf_table(ir+1,iz+1)
                endif
            endif
            
            sph_surf=   lowr_lowz_val*(1.d0-r_weight)*(1.d0-z_weight)+&
                        highr_lowz_val*(r_weight)*(1.d0-z_weight)+&
                        lowr_highz_val*(1.d0-r_weight)*(z_weight)+&
                        lowr_lowz_val*(r_weight)*(z_weight)
            
            
            
            return
        end function sph_surf

end module sph_density_functions


program test
    use cross_section_calc
    use constants
    use constant_density_functions
    use sph_density_functions
    implicit none

!     integer, parameter :: nz=21,nr=11 ! quick test
!     integer, parameter :: nz=201,nr=101 ! medium test
    integer, parameter :: nz=2001,nr=1001 ! minimum for production
    
    real(kind=8), parameter :: logTmin=0.,logTmax=4. ! temperature range for blackbodies
    integer, parameter :: Tsteps = 10

    real(kind=8), parameter :: logRmin=-4.,logRmax=0. ! log radius range for clouds (pc)
    integer, parameter :: Rsteps = 10

    real(kind=8), parameter :: logMmin=-9.,logMmax=5. ! log mass range for clouds (Msun)
    integer, parameter :: Msteps = 30

    logical, parameter :: sph_mode = .true. ! if false, then assume constant density, if .true. then use SPH density


    

!     real(kind=8), external :: constant_density,constant_density_surf ! functions for constant density
!     real(kind=8), external :: sph_density,sph_surf ! functions for sph smoothed density, assuming cutoff at r=h
    real(kind=8) :: val
    
    real(kind=8) :: this_T,this_R,this_M
    integer :: i_t,i_r,i_m
    
    real(kind=8), dimension(:,:), allocatable :: sph_surf_grid
    
    if ( sph_mode ) then
        call setup_surf(nr,nz)
!         call sph_calculate_surf_grid(sph_surf_grid)
    endif
    
    
    call load("data/kappa_sample_n4_In7.7_Te3.75.txt",1.d4)
    print *,"#T(K) R(pc) M(Msun) CrossSection(cm**2)"
    
    call set_temp(1.d1)
!     this_M=1.d-1
!     this_R=1.d-1
!     do i_r=1,Rsteps
!         this_R = 10.d0**((i_r-1)*(logRmax-logRmin)/(Rsteps-1)+logRMin)
! !         this_M = 0.1d0*this_R**2
!         do i_m=1,Msteps
!             this_M = 10.d0**((i_m-1)*(logMmax-logMmin)/(Msteps-1)+logMmin)
!             val=cross_section(this_M,this_R,sph_density,sph_surf,nr,nz)
!             print *,"mass=",val,this_M,val/this_M
!         end do
!     end do

    do i_t=1,Tsteps
        this_T = 10.d0**((i_t-1)*(logTmax-logTmin)/(Tsteps-1)+logTmin)
        call set_temp(this_T)
        do i_r=1,Rsteps
            this_R = 10.d0**((i_r-1)*(logRmax-logRmin)/(Rsteps-1)+logRMin)
    !         this_M = 0.1d0*this_R**2
            do i_m=1,Msteps
                this_M = 10.d0**((i_m-1)*(logMmax-logMmin)/(Msteps-1)+logMmin)
                if ( sph_mode ) then
                    val=cross_section(this_M,this_R,sph_density,sph_surf,nr,nz)
                else
                    val=cross_section(this_M,this_R,constant_density,constant_density_surf,nz,nr)
                endif
                print *,this_T,this_R,this_M,val,val/PI_8/this_R**2
            end do
        end do
    end do
    
end program test

