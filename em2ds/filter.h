#ifndef __FILTER__
#define __FILTER__

#include "grid2d.h"


void cscalar2d_r2c_filter( t_cscalar_grid2d * grid, const float cutoff[] );
void cvfld2d_r2c_filter( t_cvfld_grid2d * grid, const float cutoff[] );


#endif
