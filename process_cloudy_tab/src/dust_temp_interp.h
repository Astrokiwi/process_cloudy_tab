struct coolHeatDust {
    double dCool;
    double dHeat;
    double dustT;
    double opac_abs;
    double opac_scat;
    double dg;
    double column_out;
    double arad;
    
    double line_co1;
    double line_co2;
    double line_hcn1;
    double line_hcn2;
    double line_h2_1;
    double line_h2_2;
    double line_h2_3;
    
    /*int id,it,ii,ic;
    double fd,ft,fi,fc;*/
};

struct coolHeatDustArray {
    double *dCool;
    double *dHeat;
    double *dustT;
    double *opac_abs;
    double *opac_scat;
    double *dg;
    double *column_out;
    double *arad;
    
    double *line_co1;
    double *line_co2;
    double *line_hcn1;
    double *line_hcn2;
    double *line_h2_1;
    double *line_h2_2;
    double *line_h2_3;
    
    /*int id,it,ii,ic;
    double fd,ft,fi,fc;*/
};


struct AGN_heat_table {

    // heating/cooling tables under AGN
//     double *agn_heat_tab,*agn_cool_tab,*agn_arad_tab,*agn_opac_tab;
//     int agn_ntemp,agn_ndense,agn_nintensity,agn_ntau;
//     double *agn_temp_vals,*agn_dense_vals,*agn_intensity_vals,*agn_tau_vals;
//     double *tables[4];
//     int *ntabs[4];

    // heating/cooling tables under AGN (table values)
    double  *agn_heat_tab,
            *agn_cool_tab,
            *agn_dust_tab,
            *agn_dg_tab,
            *agn_opac_scat_tab,
            *agn_opac_abs_tab,
            *agn_arad_tab,
            *agn_line_co1,
            *agn_line_co2,
            *agn_line_hcn1,
            *agn_line_hcn2,
            *agn_line_h2_1,
            *agn_line_h2_2,
            *agn_line_h2_3,
            *agn_column_out_tab;
    double *tables[4];
    int *ntabs[4];
    double *lineArrays[7];


    // table axes
    int agn_ntemp,agn_ndense,agn_nintensity,agn_ncolumn_in;
    double *agn_temp_vals,*agn_dense_vals,*agn_intensity_vals,*agn_column_in_vals;


    void setupTable(const char* labelFile,const char* tableFile,bool convertLines);
    
    int agn_tab_index(int id, int it, int ii, int is);
    
};

class CoolHeatTab {
    public:
        struct coolHeatDust interpTab(double density, double temperature, double intensity, double column_in);

//         struct coolHeatDust* interpTabArray(int n0, double* density, int n1, double* temperature, int n2, double* intensity, int n3, double* column_in);
        struct coolHeatDustArray interpTabArray(int n0, double* density, int n1, double* temperature, int n2, double* intensity, int n3, double* column_in);


        //CoolHeatTab(std::string flabels,std::string ftab);
        CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab);
        CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab,const char* denseflabels,const char* denseftab);
        CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab,bool convertLines);
        CoolHeatTab(const char* flabels,const char* ftab,const char* dustlessflabels,const char* dustlessftab,const char* denseflabels,const char* denseftab,bool convertLines);
        
    private:
        static const double sputtering_temperature = 1.e5;

        int agn_tab_index(int id, int it, int ii, int is);
        int value_to_index(double edges[], int nedges, double value);
        double right_weight(double edges[], int nedges, int index, double value);
        
        // heating/cooling tables under AGN (table values)
//         double  *agn_heat_tab,
//                 *agn_cool_tab,
//                 *agn_dust_tab,
//                 *agn_dg_tab,
//                 *agn_opac_scat_tab,
//                 *agn_opac_abs_tab,
//                 *agn_arad_tab,
//                 *agn_column_out_tab;
//         double *tables[4];
//         int *ntabs[4];
// 
//         // table axes
//         int agn_ntemp,agn_ndense,agn_nintensity,agn_ncolumn_in;
//         double *agn_temp_vals,*agn_dense_vals,*agn_intensity_vals,*agn_column_in_vals;
        
        bool data_loaded,cold_dense_loaded;
        
        struct AGN_heat_table mainTable,dustlessTable,densecoldTable;
        
};
