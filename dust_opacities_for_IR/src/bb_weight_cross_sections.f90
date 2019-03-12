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
        
        real(kind=16) function cross_section(m_Msun,h_pc,density_func,surf_func,nz,nr) result(cross)
            use constants
            implicit none

            real(kind=8), intent(in) :: m_Msun,h_pc
            integer,intent(in) :: nz,nr
            real(kind=8) :: density_func,surf_func
            
            integer :: iz,ir,iw
            real(kind=8) :: r,z,dv,dr,dz
            
            real(kind=8) :: surf_norm,dens_cgs
            
            integer :: threadid
            integer, external :: OMP_GET_THREAD_NUM
            
            surf_norm = 3./4./PI_8 * m_Msun/h_pc**2 * astro_surf_density_to_cgs
            dens_cgs = 3./4./PI_8 * m_Msun/h_pc**3 * astro_density_to_cgs
            dr = 1.d0/nr
            dz = 2.d0/nz
            dv=dr*dz
!             print *,surf_func(0.,0.),density_func(0.,0.)

            ! dunno how to do trapezium in 3D, let's just do rectangles for integral
!$OMP  PARALLEL DO private(iz,ir,iw,r,z,threadid)&
!$OMP  shared(nz,nr,n_wavelengths,dz,dr,dv,blackbody,opacity,wavelength,surf_norm)&
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
!                     cross=cross+density_func(r,z)*r
                    do iw=1,n_wavelengths-1
                        cross=cross+                                    &
                                    blackbody(iw)*opacity(iw)*         &
                                    density_func(r,z)*                  &
                                    r*                                    &
                                    dexp(-opacity(iw)*surf_func(r,z)*surf_norm)*  &
                                    (wavelength(iw+1)-wavelength(iw))
                    end do
                end do
            end do
!$OMP  END PARALLEL DO
            
!             cross=cross*dens_cgs/bb_norm*dv*2.d0*PI_8
!             cross=cross*dens_cgs*dv*2.d0*PI_8*(h_pc*pc)**3/bb_norm ! in cm**2
            cross=cross*dens_cgs*dv*2.d0*PI_8*(h_pc)**3*pc/bb_norm ! in pc**2
!             print *,cross,"cm**2",log10(cross)
!             print *,cross/pc**2,"pc**2",log10(cross/pc**2)
!             print *,cross/(PI_8*(h_pc*pc)**2)," fraction absorbed"
!             print *,surf_norm,dens_cgs

            return
        end function

end module cross_section_calc

program test
    use cross_section_calc
    use constants
    implicit none
    real(kind=8), external :: constant_density,constant_density_surf ! functions
    real(kind=8) :: val
    
    real(kind=8), parameter :: logTmin=0.,logTmax=4. ! temperature range for blackbodies
    integer, parameter :: Tsteps = 5

    real(kind=8), parameter :: logRmin=-4.,logRmax=0. ! log radius range for clouds (pc)
    integer, parameter :: Rsteps = 10

    real(kind=8), parameter :: logMmin=-9.,logMmax=5. ! log mass range for clouds (Msun)
    integer, parameter :: Msteps = 30
    
    real(kind=8) :: this_T,this_R,this_M
    integer :: i_t,i_r,i_m
    
    call load("data/kappa_sample_n4_In7.7_Te3.75.txt",1.d4)
    print *,"#T(K) R(pc) M(Msun) CrossSection(cm**2)"
    do i_t=1,Tsteps
        this_T = 10.d0**((i_t-1)*(logTmax-logTmin)/(Tsteps-1)+logTmin)
        call set_temp(this_T)
        do i_r=1,Rsteps
            this_R = 10.d0**((i_r-1)*(logRmax-logRmin)/(Rsteps-1)+logRMin)
    !         this_M = 0.1d0*this_R**2
            do i_m=1,Msteps
                this_M = 10.d0**((i_m-1)*(logMmax-logMmin)/(Msteps-1)+logMmin)
                val=cross_section(this_M,this_R,constant_density,constant_density_surf,2001,1001)
                print *,this_T,this_R,this_M,val
            end do
        end do
    end do
    
end program test

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
