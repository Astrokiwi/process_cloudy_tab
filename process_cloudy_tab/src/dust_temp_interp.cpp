#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "dust_temp_interp.h"

/*// helper functions
double square(double x) {
    return x*x;
}


int square(int x) {
    return x*x;
}*/



int AGN_heat_table::agn_tab_index(int id, int it, int ii, int is) {
    return  is +
            agn_ncolumn_in*ii +
            agn_ncolumn_in*agn_nintensity*it +
            agn_ncolumn_in*agn_nintensity*agn_ntemp*id;
}

void AGN_heat_table::setupTable(const char* labelFile,const char* tableFile, bool convertLines) {

    int i,ntab,iline;
    
    // read in tables
    std::ifstream f_tablab (labelFile, std::ifstream::in);
    
    if ( !f_tablab ) {
        std::cout << "AGN table labels file not found!" << std::endl;
        std::cout << labelFile << std::endl;
        exit(0);            
    }

    f_tablab >> agn_ndense;
    agn_dense_vals = new double[agn_ndense];
    for ( i=0 ; i<agn_ndense ; i++ ) {
        f_tablab >> agn_dense_vals[i];
    }

    f_tablab >> agn_ntemp;
    agn_temp_vals = new double[agn_ntemp];
    for ( i=0 ; i<agn_ntemp ; i++ ) {
        f_tablab >> agn_temp_vals[i];
    }

    f_tablab >> agn_nintensity;
    agn_intensity_vals = new double[agn_nintensity];
    for ( i=0 ; i<agn_nintensity ; i++ ) {
        f_tablab >> agn_intensity_vals[i];
    }

    f_tablab >> agn_ncolumn_in;
    agn_column_in_vals = new double[agn_ncolumn_in];
    for ( i=0 ; i<agn_ncolumn_in ; i++ ) {
        f_tablab >> agn_column_in_vals[i];
    }
    //All.AGNOpticalDepthCutoff = agn_tau_vals[agn_ntau-1];
    
    f_tablab.close();
    
    ntab=agn_ndense*agn_ntemp*agn_nintensity*agn_ncolumn_in;
    agn_heat_tab = new double[ntab];
    agn_cool_tab = new double[ntab];
    agn_dust_tab = new double[ntab];
    agn_dg_tab = new double[ntab];
    agn_opac_scat_tab = new double[ntab];
    agn_opac_abs_tab = new double[ntab];
    agn_arad_tab = new double[ntab];
    agn_column_out_tab = new double[ntab];
    agn_line_co1 = new double[ntab];
    agn_line_co2 = new double[ntab];
    agn_line_hcn1 = new double[ntab];
    agn_line_hcn2 = new double[ntab];
    agn_line_h2_1 = new double[ntab];
    agn_line_h2_2 = new double[ntab];
    agn_line_h2_3 = new double[ntab];
    
    std::ifstream f_tabvals (tableFile, std::ifstream::in);

    if ( !f_tabvals ) {
        std::cout << "AGN table file not found!" << std::endl;
        std::cout << tableFile << std::endl;
        exit(0);            
    }

    
    for ( i=0 ; i<ntab ; i++ ) {
        f_tabvals >> agn_heat_tab[i] >> agn_cool_tab[i] >> agn_dust_tab[i] >> agn_arad_tab[i] >> agn_dg_tab[i] >>
                     agn_opac_abs_tab[i] >> agn_opac_scat_tab[i] >> agn_column_out_tab[i] >> agn_line_co1[i] >> agn_line_co2[i] >>
                     agn_line_hcn1[i] >> agn_line_hcn2[i] >> agn_line_h2_1[i] >> agn_line_h2_2[i] >> agn_line_h2_3[i];
    }
    f_tabvals.close();
    
    lineArrays[0]=agn_line_co1;
    lineArrays[1]=agn_line_co2;
    lineArrays[2]=agn_line_hcn1;
    lineArrays[3]=agn_line_hcn2;
    lineArrays[4]=agn_line_h2_1;
    lineArrays[5]=agn_line_h2_2;
    lineArrays[6]=agn_line_h2_3;
    
    if ( convertLines ) { // if we are preprocessing the tables, convert lines to erg/s/g. Output tables have this already done
            // divide lines by density to get emission per gram
            //     molecular_mass = 4./(1.+3.*.76)
            //     proton_mass_cgs = 1.6726e-24
            // log10(proton_mass_cgs*molecular_mass)=-23.6904
        int id;
        double min_line = -20;
        for ( id=0 ; id<agn_ndense ; id++ ) {
    //         double physical_density = pow(10.,agn_dense_vals[id]-23.6904);
            double physical_density = agn_dense_vals[id]-23.6904;
            for ( i=agn_tab_index(id,0,0,0) ; i<agn_tab_index(id,agn_ntemp-1,agn_nintensity-1,agn_ncolumn_in-1) ; i++ ) {
                for ( iline=0 ; iline<7 ; iline++ ) {
                    lineArrays[iline][i]=log10(lineArrays[iline][i])-physical_density;
                    if ( !isfinite(lineArrays[iline][i]) ) lineArrays[iline][i] = min_line;
                }
//                 agn_line_co1[i]=log10(agn_line_co1[i])-physical_density;
//                 agn_line_co2[i]=log10(agn_line_co2[i])-physical_density;
//                 agn_line_hcn1[i]=log10(agn_line_hcn1[i])-physical_density;
//                 agn_line_hcn2[i]=log10(agn_line_hcn2[i])-physical_density;
//                 if ( !isfinite(agn_line_co1[i]) ) agn_line_co1[i] = min_line;
//                 if ( !isfinite(agn_line_co2[i]) ) agn_line_co2[i] = min_line;
//                 if ( !isfinite(agn_line_hcn1[i])) agn_line_hcn1[i] = min_line;
//                 if ( !isfinite(agn_line_hcn2[i])) agn_line_hcn2[i] = min_line;
    //             agn_line_co1[i]/=physical_density;
    //             agn_line_co2[i]/=physical_density;
    //             agn_line_hcn1[i]/=physical_density;
    //             agn_line_hcn2[i]/=physical_density;
            }
        }
    } else { // log for better interpolation
        for ( i=0 ; i<ntab ; i++ ) {
            for ( iline=0 ; iline<7 ; iline++ ) {
                lineArrays[iline][i]=log10(lineArrays[iline][i]);
            }
//             agn_line_co1[i]=log10(agn_line_co1[i]);
//             agn_line_co2[i]=log10(agn_line_co2[i]);
//             agn_line_hcn1[i]=log10(agn_line_hcn1[i]);
//             agn_line_hcn2[i]=log10(agn_line_hcn2[i]);
        }
    }
    
    // set up pointers for better looping
    tables[0] = agn_dense_vals;
    ntabs[0] = &agn_ndense;
    tables[1] = agn_temp_vals;
    ntabs[1] = &agn_ntemp;
    tables[2] = agn_intensity_vals;
    ntabs[2] = &agn_nintensity;
    tables[3] = agn_column_in_vals;
    ntabs[3] = &agn_ncolumn_in;

}

/*! find index in (short) table where value fits to the right of
return -1 if it's to the left of everything, and n-1 if it's to the right of everything
*/
int CoolHeatTab::value_to_index(double edges[], int nedges, double value) {
    int i;
    
    if ( value<edges[0] ) return -1;
    for ( i=1 ; i<nedges ; i++ ) {
        if ( value<edges[i] ) {
            return i-1;
        }
    }
    return nedges-1;
}


/*! weighting for interpolation
No error checking!
*/
double CoolHeatTab::right_weight(double edges[], int nedges, int index, double value) {
    if ( index<0 ) return 1.;
    if ( index>=nedges-1) return 0.;
    
    return (value-edges[index])/(edges[index+1]-edges[index]);
}

// int CoolHeatTab::agn_tab_index(int id, int it, int ii, int is) {
//     return  is +
//             agn_ncolumn_in*ii +
//             agn_ncolumn_in*agn_nintensity*it +
//             agn_ncolumn_in*agn_nintensity*agn_ntemp*id;
// };
// 
// struct[] coolHeatDust CoolHeatTab::interpTabVec(double density[], double temperature[], double intensity[], double coldens[], int n) {
//     struct outp[n];
//     
//     for ( int ii=0 ; ii<n ; ii++ ) {
//         outp[ii] = interpTab(density[ii],temperature[ii],intensity[ii],coldens[ii]);
//     }
//     
//     return outp;
// } 

// struct coolHeatDust* CoolHeatTab::interpTabArray(int n0, double* density, int n1, double* temperature, int n2, double* intensity, int n3, double* column_in) {
struct coolHeatDustArray CoolHeatTab::interpTabArray(int n0, double* density, int n1, double* temperature, int n2, double* intensity, int n3, double* column_in) {
//     if ( n0!=n1 || n1!=n2 || n2!=n3 ) {
//         return NULL;
//     }

    int i;
    
    coolHeatDustArray dustStructs;
    
    dustStructs.dCool = new double[n0];
    dustStructs.dHeat = new double[n0];
    dustStructs.dustT = new double[n0];
    dustStructs.opac_abs = new double[n0];
    dustStructs.opac_scat = new double[n0];
    dustStructs.dg = new double[n0];
    dustStructs.column_out = new double[n0];
    dustStructs.arad = new double[n0];
    
    dustStructs.line_co1 = new double[n0];
    dustStructs.line_co2 = new double[n0];
    dustStructs.line_hcn1 = new double[n0];
    dustStructs.line_hcn2 = new double[n0];
    dustStructs.line_h2_1 = new double[n0];
    dustStructs.line_h2_2 = new double[n0];
    dustStructs.line_h2_3 = new double[n0];
    
//     struct coolHeatDust *dustStructs = new coolHeatDust[n0];
//     struct coolHeatDust *dustStructs = (struct coolHeatDust *)malloc(sizeof(struct coolHeatDust)*n0);
    for ( i=0 ; i<n0 ; i++ ) {
        struct coolHeatDust nextEntry = interpTab(density[i],temperature[i],intensity[i],column_in[i]);

        dustStructs.dCool[i] = nextEntry.dCool;
        dustStructs.dHeat[i] = nextEntry.dHeat;
        dustStructs.dustT[i] = nextEntry.dustT;
        dustStructs.opac_abs[i] = nextEntry.opac_abs;
        dustStructs.opac_scat[i] = nextEntry.opac_scat;
        dustStructs.dg[i] = nextEntry.dg;
        dustStructs.column_out[i] = nextEntry.column_out;
        dustStructs.arad[i] = nextEntry.arad;
        dustStructs.line_co1[i] = nextEntry.line_co1;
        dustStructs.line_co2[i] = nextEntry.line_co2;
        dustStructs.line_hcn1[i] = nextEntry.line_hcn1;
        dustStructs.line_hcn2[i] = nextEntry.line_hcn2;
        dustStructs.line_h2_1[i] = nextEntry.line_h2_1;
        dustStructs.line_h2_2[i] = nextEntry.line_h2_2;
        dustStructs.line_h2_3[i] = nextEntry.line_h2_3;
    }
    return dustStructs;
}

// takes input values in *cgs units*
struct coolHeatDust CoolHeatTab::interpTab(double density, double temperature, double intensity, double column_in) {
    if ( !this->data_loaded ) {
        struct coolHeatDust garbage;
        return garbage;
    }

    int jj,kk;
    
    int tab_index[4];
    int interp_indices[4][2];
    double interp_weights[4][2];
    double values[4];

    bool interpBetweenTables = false; // only use one table at a time
    bool dusty;
    
    AGN_heat_table *thisTable; 
    
    if ( temperature<=sputtering_temperature ) {
        thisTable = &mainTable;
        dusty = true;
    } else {
        thisTable = &dustlessTable;
        dusty = false;
    }
    
    values[0] = log10(density);
    values[1] = log10(temperature);
    values[2] = log10(intensity);
    values[3] = column_in;
    
    if ( this->cold_dense_loaded ) {
        if ( values[0]>thisTable->agn_dense_vals[thisTable->agn_ndense-1] && values[1]<densecoldTable.agn_temp_vals[densecoldTable.agn_ndense-1] && dusty ) {
            if (values[0]>=densecoldTable.agn_dense_vals[0] ) {
                // only need to use cold low density table
                thisTable = &densecoldTable;
            } else {
                // need to interpolate between tables
                interpBetweenTables = true;
            }    
        }
    }
        
//     interpBetweenTables = false; // for test
//     thisTable = &mainTable; // also for test
    

    for (jj=interpBetweenTables?1:0 ; jj<4 ; jj++ ) {
        tab_index[jj] = value_to_index(thisTable->tables[jj],*thisTable->ntabs[jj],values[jj]);
        interp_indices[jj][0] = std::max(tab_index[jj],0);
        interp_indices[jj][1] = std::min(tab_index[jj]+1,*thisTable->ntabs[jj]-1);
        interp_weights[jj][1] = right_weight(thisTable->tables[jj],*thisTable->ntabs[jj],tab_index[jj],values[jj]);
        interp_weights[jj][0] = 1.-interp_weights[jj][1];
    }
    double heat_interp=0.,cool_interp=0.,dust_interp=0.,dg_interp=0.,opac_abs_interp=0.,
           opac_scat_interp=0.,column_out_interp=0.,arad_interp=0.;           
//            ,co1_interp=0.,co2_interp=0.,
//            hcn1_interp=0.,hcn2_interp=0.;
    double line_interp[7];
    int iline;
    for ( iline=0 ; iline<7 ; iline++ ) {
        line_interp[iline]=0.;
    }

    
    if ( interpBetweenTables ) {
//         if ( ThisTask==0 ) std::cout << "interpbetweentables" << std::endl;
        interp_indices[0][0] = thisTable->agn_ndense-1;
        interp_indices[0][1] = 0;
        interp_weights[0][1] = (values[0]-thisTable->agn_dense_vals[thisTable->agn_ndense-1])
                                /(densecoldTable.agn_dense_vals[0]-thisTable->agn_dense_vals[thisTable->agn_ndense-1]);
        interp_weights[0][0] = 1.-interp_weights[0][1];
        
        int interp_indices_lowtemp[2];
        if ( interp_indices[1][1]>densecoldTable.agn_ntemp-1 ) {
            for (jj=0 ; jj<2 ; jj++) {
                interp_indices_lowtemp[jj] = densecoldTable.agn_ntemp-1;
            }
        } else {
            for (jj=0 ; jj<2 ; jj++) {
                interp_indices_lowtemp[jj] = interp_indices[1][jj];
            }
        }
        
        AGN_heat_table *curTable;
        for (jj=0 ; jj<16 ; jj++ ) {
            int idex;
            
            if ( jj%2==0 ) {
                curTable = thisTable;
                idex = curTable->agn_tab_index(interp_indices[0][jj%2],interp_indices[1][(jj/2)%2],interp_indices[2][(jj/4)%2],interp_indices[3][(jj/8)%2]);
            } else {
                curTable = &densecoldTable;
                idex = curTable->agn_tab_index(interp_indices[0][jj%2],interp_indices_lowtemp[(jj/2)%2],interp_indices[2][(jj/4)%2],interp_indices[3][(jj/8)%2]);
            }

            double weight = 1.;
            for ( kk=0 ; kk<4 ; kk++ ) {
                weight*=interp_weights[kk][(jj/(1<<kk))%2];
            }

            heat_interp+=weight*curTable->agn_heat_tab[idex];
            cool_interp+=weight*curTable->agn_cool_tab[idex];
            dust_interp+=weight*curTable->agn_dust_tab[idex];
            dg_interp+=weight*curTable->agn_dg_tab[idex];
            opac_scat_interp+=weight*curTable->agn_opac_scat_tab[idex];
            opac_abs_interp+=weight*curTable->agn_opac_abs_tab[idex];
            column_out_interp+=weight*curTable->agn_column_out_tab[idex];
            arad_interp+=weight*curTable->agn_arad_tab[idex];
            for ( iline=0 ; iline<7 ; iline++ ) {
                line_interp[iline]+=weight*curTable->lineArrays[iline][idex];
            }

//             co1_interp+=weight*curTable->agn_line_co1[idex];
//             co2_interp+=weight*curTable->agn_line_co2[idex];
//             hcn1_interp+=weight*curTable->agn_line_hcn1[idex];
//             hcn2_interp+=weight*curTable->agn_line_hcn2[idex];
        }
    } else {
        for (jj=0 ; jj<16 ; jj++ ) {
            int idex = thisTable->agn_tab_index(interp_indices[0][jj%2],interp_indices[1][(jj/2)%2],interp_indices[2][(jj/4)%2],interp_indices[3][(jj/8)%2]);
            double weight = 1.;
            for ( kk=0 ; kk<4 ; kk++ ) {
                weight*=interp_weights[kk][(jj/(1<<kk))%2];
            }
            heat_interp+=weight*thisTable->agn_heat_tab[idex];
            cool_interp+=weight*thisTable->agn_cool_tab[idex];
            dust_interp+=weight*thisTable->agn_dust_tab[idex];
            dg_interp+=weight*thisTable->agn_dg_tab[idex];
            opac_scat_interp+=weight*thisTable->agn_opac_scat_tab[idex];
            opac_abs_interp+=weight*thisTable->agn_opac_abs_tab[idex];
            column_out_interp+=weight*thisTable->agn_column_out_tab[idex];
            arad_interp+=weight*thisTable->agn_arad_tab[idex];
//             co1_interp+=weight*thisTable->agn_line_co1[idex];
//             co2_interp+=weight*thisTable->agn_line_co2[idex];
//             hcn1_interp+=weight*thisTable->agn_line_hcn1[idex];
//             hcn2_interp+=weight*thisTable->agn_line_hcn2[idex];
            for ( iline=0 ; iline<7 ; iline++ ) {
                line_interp[iline]+=weight*thisTable->lineArrays[iline][idex];
            }
        }
    
    }
    
    // do extrapolation for density above max
    if ( values[0]>thisTable->agn_dense_vals[thisTable->agn_ndense-1] && !interpBetweenTables) {
//         if ( ThisTask==0 ) std::cout << "extrapolation" << std::endl;
        heat_interp += values[0]-thisTable->agn_dense_vals[thisTable->agn_ndense-1]; // heating per unit volume is proportional to density - i.e. heating per unit mass is constant
        cool_interp += 2.*(values[0]-thisTable->agn_dense_vals[thisTable->agn_ndense-1]); // cooling per volume is proportional to n^2
    }
    
    // do extrapolation for intensity above max - assume heating is proportional to intensity outside of the table
    // do extrapolation for intensity below minimum
    if ( values[2]<thisTable->agn_intensity_vals[0] ) {
        heat_interp += values[2] - thisTable->agn_intensity_vals[0];
    }
    // do extrapolation for intensity above maximum
    if ( values[2]>thisTable->agn_intensity_vals[thisTable->agn_nintensity-1] ) {
        heat_interp += values[2] - thisTable->agn_intensity_vals[thisTable->agn_nintensity-1];
    }


    // NO EXTRAPOLATION yet
    struct coolHeatDust outp;
    outp.dCool = cool_interp;
    outp.dHeat = heat_interp;
    outp.dustT = dust_interp;
    outp.opac_abs = opac_abs_interp;
    outp.opac_scat = opac_scat_interp;
    outp.dg = dg_interp;
    outp.column_out = column_out_interp;
    outp.arad = arad_interp;
//     outp.line_co1 = pow(10.,co1_interp);
//     outp.line_co2 = pow(10.,co2_interp);
//     outp.line_hcn1 = pow(10.,hcn1_interp);
//     outp.line_hcn2 = pow(10.,hcn2_interp);
    outp.line_co1 = pow(10.,line_interp[0]);
    outp.line_co2 = pow(10.,line_interp[1]);
    outp.line_hcn1 = pow(10.,line_interp[2]);
    outp.line_hcn2 = pow(10.,line_interp[3]);
    outp.line_h2_1 = pow(10.,line_interp[4]);
    outp.line_h2_2 = pow(10.,line_interp[5]);
    outp.line_h2_3 = pow(10.,line_interp[6]);
//     outp.line_co1 = co1_interp;
//     outp.line_co2 = co2_interp;
//     outp.line_hcn1 = hcn1_interp;
//     outp.line_hcn2 = hcn2_interp;
    return outp;
}

//CoolHeatTab::CoolHeatTab(std::string flabels,std::string ftab) {
CoolHeatTab::CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab, bool convertLines) {
    this->data_loaded = false;
    this->cold_dense_loaded = false;

    mainTable.setupTable(flabels,ftab,convertLines);
    dustlessTable.setupTable(dustlessflabels,dustlessftab,convertLines);

    this->data_loaded = true;

}

CoolHeatTab::CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab) {
    this->data_loaded = false;
    this->cold_dense_loaded = false;

    mainTable.setupTable(flabels,ftab,false);
    dustlessTable.setupTable(dustlessflabels,dustlessftab,false);

    this->data_loaded = true;

}

CoolHeatTab::CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab,const char* denseflabels,const char* denseftab, bool convertLines) {
    this->data_loaded = false;
    this->cold_dense_loaded = false;

    mainTable.setupTable(flabels,ftab,convertLines);
    dustlessTable.setupTable(dustlessflabels,dustlessftab,convertLines);
    densecoldTable.setupTable(denseflabels,denseftab,convertLines);

    this->cold_dense_loaded = true;
    this->data_loaded = true;
}


CoolHeatTab::CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab,const char* denseflabels,const char* denseftab) {
    this->data_loaded = false;
    this->cold_dense_loaded = false;

    mainTable.setupTable(flabels,ftab,false);
    dustlessTable.setupTable(dustlessflabels,dustlessftab,false);
    densecoldTable.setupTable(denseflabels,denseftab,false);

    this->cold_dense_loaded = true;
    this->data_loaded = true;
}


// int main(int argc, char **argv) {
//     CoolHeatTab *testTable = new CoolHeatTab("coolheat_tab_marta/shrunk_table_labels_281118tau.dat",
//                                         "coolheat_tab_marta/shrunk_table_281118_m0.0001_hsmooth_tau.dat",
//                                         "coolheat_tab_marta/shrunk_table_labels_281118taunodust.dat",
//                                         "coolheat_tab_marta/shrunk_table_281118_m0.0001_hsmooth_taunodust.dat",
//                                         "coolheat_tab_marta/shrunk_table_labels_281118taudense.dat",
//                                         "coolheat_tab_marta/shrunk_table_281118_m0.0001_hsmooth_taudense.dat");
//     struct coolHeatDust outp;
//     
//     outp = testTable->interpTab(1.e3,3.e1,1.e8,5.);
//     
//     std::cout << "Heat:" << outp.dHeat << " Cool:" << outp.dCool << " Dust:" << outp.dustT << std::endl;
//     std::cout << "co1:" << outp.line_co1 << " hcn1:" << outp.line_hcn1 << " h2_1:" << outp.line_h2_1 << std::endl;
//     std::cout << "co2:" << outp.line_co2 << " hcn2:" << outp.line_hcn2 << " h2_2:" << outp.line_h2_2 << std::endl;
// }


