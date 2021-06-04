%module tab_interp

%{
    #define SWIG_FILE_WITH_INIT

    #include <stdlib.h>
    #include <string>
    #include <fstream>
    #include <stdio.h>
    #include <iostream>
    #include <math.h>

    #include "dust_temp_interp.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int n0, double* density), (int n1, double* temperature), (int n2, double* intensity), (int n3, double* column_in)};

# %include <carrays.i>
# %array_class(struct coolHeatDust, coolHeatDustArray);


%include "dust_temp_interp.h"
