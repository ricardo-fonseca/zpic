/**
 * @file particles.c
 * @author Ricardo Fonseca
 * @brief Particle species
 * @version 0.2
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "grid2d.h"

#include "particles.h"

#include "random.h"
#include "emf.h"
#include "current.h"

#include "zdf.h"
#include "timer.h"

static double _spec_time = 0.0;
static uint64_t _spec_npush = 0;

void spec_sort( t_species *spec );

/**
 * @brief Returns the total time spent pushing particles (includes boundaries and moving window)
 * 
 * @return  Total time in seconds
 */
double spec_time( void )
{
    return _spec_time;
}

/**
 * @brief Returns the total number of particle pushes
 * 
 * @return  Number of particle pushes
 */
uint64_t spec_npush( void )
{
    return _spec_npush;
}

/**
 * @brief Returns the performance achieved by the code (push time)
 * 
 * @return  Performance in seconds per particle
 */
double spec_perf( void )
{
    return ( _spec_npush > 0 )? _spec_time / _spec_npush : -1.0;
}

/*********************************************************************************************
 
 Initialization
 
 *********************************************************************************************/

/**
 * @brief Dummy custom density function that always return 1
 * 
 * This is used when a custom density profile is chosen, but no function is supplied
 * 
 * @param   x       Position
 * @param   data    Pointer to custom data (used when the custom function is defined in Python)
 * @return          Always returns 1.0
 */
float one( float x, void *data ) {
    return 1.0;
}

/**
 * @brief Sets the momentum of the range of particles supplieds using a thermal distribution
 * 
 * @param spec  Particle species
 * @param start Index of the first particle to set the momentum
 * @param end   Index of the last particle to set the momentum
 */
void spec_set_u( t_species* spec, const int start, const int end )
{

#if 0

    for (int i = start; i <= end; i++) {
        spec->part[i].ux = spec -> ufl[0] + spec -> uth[0] * rand_norm(); 
        spec->part[i].uy = spec -> ufl[1] + spec -> uth[1] * rand_norm(); 
        spec->part[i].uz = spec -> ufl[2] + spec -> uth[2] * rand_norm(); 
    }

#else
    // Initialize thermal component
    for (int i = start; i <= end; i++) {
        spec->part[i].ux = spec -> uth[0] * rand_norm(); 
        spec->part[i].uy = spec -> uth[1] * rand_norm(); 
        spec->part[i].uz = spec -> uth[2] * rand_norm(); 
    }

    // Calculate net momentum in each cell
    const int size = spec->nx[0] * spec->nx[1];
    const int stride = spec->nx[1];

    float3 * restrict net_u = (float3 *) malloc( size * sizeof(float3));
    int * restrict    npc   = (int *) malloc( size * sizeof(int));

    // Zero momentum grids
    memset(net_u, 0, size * sizeof(float3) );
    memset(npc, 0, size * sizeof(int) );

    // Accumulate momentum in each cell
    for (int i = start; i <= end; i++) {
        const int idx  = spec -> part[i].ix + stride * spec -> part[i].iy ;

        net_u[ idx ].x += spec->part[i].ux;
        net_u[ idx ].y += spec->part[i].uy;
        net_u[ idx ].z += spec->part[i].uz;

        npc[ idx ] += 1;
    }

    // Normalize to the number of particles in each cell to get the
    // average momentum in each cell
    for(int i =0; i< size; i++ ) {
        const float norm = (npc[ i ] > 0) ? 1.0f/npc[i] : 0;

        net_u[ i ].x *= norm;
        net_u[ i ].y *= norm;
        net_u[ i ].z *= norm;
    }

    // Subtract average momentum and add fluid component
    for (int i = start; i <= end; i++) {
        const int idx  = spec -> part[i].ix + stride * spec -> part[i].iy ;

        spec->part[i].ux += spec -> ufl[0] - net_u[ idx ].x;
        spec->part[i].uy += spec -> ufl[1] - net_u[ idx ].y;
        spec->part[i].uz += spec -> ufl[2] - net_u[ idx ].z;
    }

    // Free temporary memory
    free( npc );
    free( net_u );

#endif

}

/**
 * @brief Sets initial position of particles according to density profile
 * 
 * Note that particle momentum is not set
 * 
 * @param spec      Particle species
 * @param range     Cell range in which to inject the particles
 */
void spec_set_x( t_species* spec, const int range[][2] )
{

    int i, j, k, ip;
    
    float* poscell;
    float start, end;
    
    // Calculate particle positions inside the cell
    const int npc = spec->ppc[0]*spec->ppc[1];
    
    poscell = malloc( 2 * npc * sizeof( float ) );
    ip = 0;
    for (j =0; j<spec->ppc[1]; j++) {
        for (i=0; i<spec->ppc[0]; i++) {
            poscell[ip]   = (1 + 2*i - spec->ppc[0]) / (2.0*spec->ppc[0]);
            poscell[ip+1] = (1 + 2*j - spec->ppc[1]) / (2.0*spec->ppc[1]);
            ip+=2;
        }
    }

    ip = spec -> np;
    
    // Set position of particles in the specified grid range according to the density profile
    switch ( spec -> density.type ) {
    case STEP: // Step like density profile
        
        // Get edge position normalized to cell size;
        start = spec -> density.start / spec -> dx[0];

        for (j = range[1][0]; j <= range[1][1]; j++) {
            for (i = range[0][0]; i <= range[0][1]; i++) {

                for (k=0; k<npc; k++) {
                    if ( i + poscell[2*k] > start ) {
                        spec->part[ip].ix = i;
                        spec->part[ip].iy = j;
                        spec->part[ip].x = poscell[2*k];
                        spec->part[ip].y = poscell[2*k+1];
                        ip++;
                    }
                }
            }
        }
        break;

    case SLAB: // Slab like density profile
        
        // Get edges positions normalized to cell size;
        start = spec -> density.start / spec -> dx[0] - 0.5;
        end   = spec -> density.end / spec -> dx[0] - 0.5;

        for (j = range[1][0]; j <= range[1][1]; j++) {
            for (i = range[0][0]; i <= range[0][1]; i++) {

                for (k=0; k<npc; k++) {
                    if ( i + poscell[2*k] > start &&  i + poscell[2*k] < end ) {
                        spec->part[ip].ix = i;
                        spec->part[ip].iy = j;
                        spec->part[ip].x = poscell[2*k];
                        spec->part[ip].y = poscell[2*k+1];
                        ip++;
                    }
                }
            }
        }
        break;
    case CUSTOM:
        {
            
            const double dx = spec -> dx[0];
            const double dy = spec -> dx[1];

            const double cppx = 1.0 / spec->ppc[0];
            const double cppy = 1.0 / spec->ppc[1];

            // Threshold density for particle injection
            double thresh = 4 * cppx * cppy;
            
            int ix0 = range[0][0];
            int iy0 = range[1][0];

            // Initial injection parameters along x

            // Either start from the bottom of the box or continue
            // from previous injection in a moving window
            unsigned long kx;
            double d0x, d1x;
            kx  = spec -> density.custom_x_total_part;
            d1x = spec -> density.custom_x_total_q;

            // X density at the lower injection point 
            // will be copied onto n0x later
            double n0x, n1x;
            n1x = (*spec -> density.custom_x)
                (ix0 * dx, spec -> density.custom_data_x);

            // Initial injection parameters along y

            // Always starts from the bottom of the box
            unsigned long const ky0 = 0;
            double const d0y0 = 0;
                        
            // Y density at the lower injection point
            double n0y0 = (*spec -> density.custom_y)
                (iy0 * dy, spec -> density.custom_data_y);
            
            // Loop over x cells
            for( int ix = ix0; ix <= range[0][1]; ix++ ) {
                // Get x density on the edges of current cell
                n0x = n1x;
                n1x = (*spec -> density.custom_x)
                    ((ix + 1) * dx, spec -> density.custom_data_x);
                // Get cumulative x density on the edges of current cell
                d0x = d1x;
                d1x += 0.5 * (n0x+n1x);

                double Rsx;
                while( (Rsx = (kx+0.5)*cppx) < d1x) {
                    // x position of particles to inject
                    double x = 2 * (Rsx-d0x) /( sqrt( n0x*n0x + 2 * (n1x-n0x) * (Rsx-d0x) ) + n0x );

                    double nx = (0.5-x)*n0x + (0.5+x)*n1x;

                    // Find y position of particles to inject
                    int ky = ky0;

                    double n0y;
                    double n1y = n0y0;

                    double d0y;
                    double d1y = d0y0;
                    
                    for( int iy = iy0; iy <= range[1][1]; iy++ ) {
                        n0y = n1y;
                        n1y = (*spec -> density.custom_y)
                                ((iy+1) * dy, spec -> density.custom_data_y);
                        d0y = d1y;
                        d1y += 0.5*(n0y+n1y);

                        double Rsy;
                        while( (Rsy = (ky+0.5)*cppy) < d1y) {
                            double y = 2 * (Rsy-d0y) /( sqrt( n0y*n0y + 2 * (n1y-n0y) * (Rsy-d0y) ) + n0y );
                            
                            double ny = (0.5-y)*n0y + (0.5+y)*n1y;
                           
                            if ( nx*ny > thresh ) {
                                spec->part[ip].ix = ix;
                                spec->part[ip].iy = iy;
                                spec->part[ip].x = x;
                                spec->part[ip].y = y;
                                ip++;
                            }
                            // Move to next y position
                            ky++;
                        }
                    }
                    // Move to next x position
                    kx++;
                }
            }

            spec -> density.custom_x_total_q = d1x;
            spec -> density.custom_x_total_part = kx;

        }
        break;

    case EMPTY: // Empty profile
        break;

    default: // Uniform density
        for (j = range[1][0]; j <= range[1][1]; j++) {
            for (i = range[0][0]; i <= range[0][1]; i++) {

                for (k=0; k<npc; k++) {
                    spec->part[ip].ix = i;
                    spec->part[ip].iy = j;
                    spec->part[ip].x = poscell[2*k];
                    spec->part[ip].y = poscell[2*k+1];
                    ip++;
                }
            }
        }
    }
    
    spec -> np = ip;
    
    free(poscell);
    
}

/**
 * @brief Gets number of particles to be injected.
 *
 * Calculates the number of particles to be injected in the specified range according
 * to the specified density profile. The returned value is not exact but it is
 * guaranteed to be larger than the actual number of particles to be injected
 *
 * @param spec          Particle species
 * @param range         Range of cells in which to inject [x,y][lower,upper]
 * @return              Number of particles to be injected
 */
int spec_np_inj( t_species* spec, const int range[][2] )
{
    int np_inj;

    switch ( spec -> density.type ) {
    case STEP: // Step like density profile
        {
            int i0 = spec -> density.start / spec -> dx[0];

            if ( i0 > range[0][1] ) {
                np_inj = 0;
            } else {
                if ( i0 < range[0][0] ) i0 = range[0][0];
                np_inj = ( range[0][1] - i0 + 1 ) * spec -> ppc[0];
            }

            np_inj *= ( range[1][1] - range[1][0] + 1 ) * spec -> ppc[1];
        }
        break;

    case SLAB: // Slab like density profile
        {
            int i0 = spec -> density.start / spec -> dx[0];
            int i1 = spec -> density.end / spec -> dx[1];

            if ( (i0 > range[0][1]) || (i1 < range[0][0]) ) {
                np_inj = 0;
            } else {
                if ( i0 < range[0][0] ) i0 = range[0][0];
                if ( i1 > range[0][1] ) i1 = range[0][1];
                np_inj = ( i1 - i0 + 1 ) * spec -> ppc[0];
            }

            np_inj *= ( range[1][1] - range[1][0] + 1 ) * spec -> ppc[1];
        }
        break;

    case CUSTOM: // custom density profile
        {
            // Integrate total charge along x
            double x = range[0][0] * spec->dx[0];
            double qx = (*spec -> density.custom_x)(x,spec -> density.custom_data_x);
            
            x = (range[0][1] + 1) * spec->dx[0];
            qx += (*spec -> density.custom_x)(x,spec -> density.custom_data_x);
            
            qx *= 0.5;

            for( int i = range[0][0]+1; i <= range[0][1]; i++) {
                x = i * spec->dx[0];
                qx += (*spec -> density.custom_x)(x,spec -> density.custom_data_x);
            }

            // Integrate total charge along y
            double y = range[1][0] * spec->dx[1];
            double qy = (*spec -> density.custom_y)(y,spec -> density.custom_data_y);
            
            y = (range[1][1]+1) * spec->dx[1];
            qy += (*spec -> density.custom_y)(y,spec -> density.custom_data_y);
            
            qy *= 0.5;

            for( int j = range[1][0]+1; j <= range[1][1]; j++) {
                y = j * spec->dx[1];
                qy += (*spec -> density.custom_y)(y,spec -> density.custom_data_y);
            }

            // Get corresponding number of simulation particles, rounding up
            np_inj = ceil(qx * spec -> ppc[0]) * ceil(qy * spec -> ppc[1]);
        }
        break;

    case EMPTY: // Empty profile
        np_inj = 0;
        break;

    default: // Uniform density
        np_inj = ( range[0][1] - range[0][0] + 1 ) * spec -> ppc[0] *
                 ( range[1][1] - range[1][0] + 1 ) * spec -> ppc[1];
    }

    return np_inj;
}

/**
 * @brief Grows particle buffer to specified size.
 * 
 * If the new size is smaller than the previous size the buffer size is not changed
 * and the function returns silently.
 * 
 * @param   spec    Particle species
 * @param   size    New buffer size (will be rounded up to next multiple of 1024)
 */
void spec_grow_buffer( t_species* spec, const int size ) {
    if ( size > spec -> np_max ) {
        // Increase by chunks of 1024 particles
        spec -> np_max = ( size/1024 + 1) * 1024;
        spec -> part = realloc( (void*) spec -> part, spec -> np_max * sizeof(t_part) );
    }
}

/**
 * @brief Inject new particles into the specified grid range.
 * 
 * The particle buffer will be grown if required
 * 
 * @param   spec    Particle species
 * @param   range   Grid range [2][i0, i1] where to inject particles
 **/
void spec_inject_particles( t_species* spec, const int range[][2] )
{
    int start = spec -> np;

    // Get maximum number of particles to inject
    int np_inj = spec_np_inj( spec, range );

    // Check if buffer is large enough and if not reallocate
    spec_grow_buffer( spec, spec -> np + np_inj );

    // Set particle positions
    spec_set_x( spec, range );
    
    // Set momentum of injected particles
    spec_set_u( spec, start, spec -> np - 1 );

}

/**
 * @brief Initialize particle Species object
 * 
 * This routine will also inject the initial particle distribution,
 * setting thermal/fluid velocities.
 *  
 * @param spec      Particle species
 * @param name      Name for the species (used for diagnostic output)
 * @param m_q       Mass over charge ratio for species, in simulation units
 * @param ppc       Reference number of particles per cell [x,y]
 * @param ufl       Initial fluid momentum of particles, may be set to NULL 
 * @param uth       Initial thermal momentum of particles, may be set to NULL
 * @param nx        Number of grid points [x,y]
 * @param box       Simulation box size [x,y] in simulation units
 * @param dt        Simulation time step, in simulation units
 * @param density   Density profile for particle injection, may be set to NULL
 */
void spec_new( t_species* spec, char name[], const float m_q, const int ppc[], 
              const float *ufl, const float * uth,
              const int nx[], float box[], const float dt, t_density* density )
{

    int i, npc;
    
    // Species name
    strncpy( spec -> name, name, MAX_SPNAME_LEN );
    
    npc = 1;
    // Store species data
    for (i=0; i<2; i++) {
        spec->nx[i] = nx[i];
        spec->ppc[i] = ppc[i];
        npc *= ppc[i];
        
        spec->box[i] = box[i];
        spec->dx[i] = box[i] / nx[i];
    }
    
    spec -> m_q = m_q;
    spec -> q = copysign( 1.0f, m_q ) / npc;

    spec -> dt = dt;
    
    // Initialize particle buffer
    spec->np_max = 0;
    spec->part = NULL;
    
    
    // Initialize density profile
    if ( density ) {
        spec -> density = *density;
        if ( spec -> density.n == 0. ) spec -> density.n = 1.0;
    } else {
        // Default values
        spec -> density = (t_density) { .type = UNIFORM, .n = 1.0 };
    }

    // Density multiplier
    spec ->q *= fabsf( spec -> density.n );

    // Initialize temperature profile
    if ( ufl ) {
        for(i=0; i<3; i++) spec -> ufl[i] = ufl[i];
    } else {
        for(i=0; i<3; i++) spec -> ufl[i] = 0;
    }

    if ( uth ) {
        for(i=0; i<3; i++) spec -> uth[i] = uth[i];
    } else {
        for(i=0; i<3; i++) spec -> uth[i] = 0;
    }

    // Reset iteration number
    spec -> iter = 0;

    // Inject initial particle distribution
    spec -> np = 0;
    
    const int range[][2] = {{0, nx[0]-1},
                            {0, nx[1]-1}};

    spec_inject_particles( spec, range );

    // Set default sorting frequency
    spec -> n_sort = 16;

}

/**
 * @brief Frees dynamic memory from particle species
 * 
 * @param spec Particle species
 */
void spec_delete( t_species* spec )
{
    free(spec->part);
    spec->np = -1;
}


/*********************************************************************************************
 
 Charge / Current deposition
 
 *********************************************************************************************/

/**
 * @brief Returns number of cells moved
 * 
 * Note that the particle will move at most 1 cell in either direction
 * 
 * @param x         End particle position, normalized to cell size
 * @return ltrim    Number of cells moved, {-1,0,1}
 */
int ltrim( float x )
{
    return (( x >= 0.5f )?1:0) - (( x < -0.5f )?1:0);
}

/**
 * @brief Deposit single particle charge
 * 
 * @param rho   Charge density grid
 * @param part  Particle data
 * @param q     Species charge per particle
 */
void deposit_charge( t_scalar_grid2d * rho, const t_part* restrict const part, const float q )
{
    float s0x, s1x;
    float s0y, s1y;

    const int nrow = rho -> nrow;
    const int idx  = part->ix + part->iy * nrow;
    
    s0x = 0.5f - part->x;
    s1x = 0.5f + part->x; 

    s0y = 0.5f - part->y;
    s1y = 0.5f + part->y; 

    rho->s[ idx            ] += s0y * s0x * q;
    rho->s[ idx        + 1 ] += s0y * s1x * q;
    rho->s[ idx + nrow     ] += s1y * s0x * q;
    rho->s[ idx + nrow + 1 ] += s1y * s1x * q;

}

/**
 * @brief Deposit single particle current
 * 
 * The particle position is advanced half time step and current is deposited
 * using that position
 * 
 * @param J     Current density grid
 * @param part  Particle data
 * @param q     Species charge per particle
 * @param rg    Particle $1 / \gamma$
 * @param dx    x cell size
 * @param dy    y cell size
 */
void deposit_current( t_float3_grid2d* J, const t_part* restrict const part, const float q, const float rg, 
    const float dx, const float dy )
{
    int i, di;
    int j, dj;
    float x, y;
    float s0x, s1x, s0y, s1y;
    float jx, jy, jz;

    const int nrow = J -> nrow;

    // Find position time centered with velocity
    i = part->ix;
    x = part->x + 0.5f * dx;
    di = ltrim(x);
    i += di;
    x -= di;

    j = part->iy;
    y = part->y + 0.5f * dy;
    dj = ltrim(y);
    j += dj;
    y -= dj;
    
    s0x = 0.5f - x;
    s1x = 0.5f + x; 

    s0y = 0.5f - y;
    s1y = 0.5f + y; 

    jx = q * part -> ux * rg;
    jy = q * part -> uy * rg;
    jz = q * part -> uz * rg;

    int idx = i + j * nrow;

    J->x[idx           ] += s0y * s0x * jx;
    J->y[idx           ] += s0y * s0x * jy;
    J->z[idx           ] += s0y * s0x * jz;

    J->x[idx        + 1] += s0y * s1x * jx;
    J->y[idx        + 1] += s0y * s1x * jy;
    J->z[idx        + 1] += s0y * s1x * jz;

    J->x[idx + nrow    ] += s1y * s0x * jx;
    J->y[idx + nrow    ] += s1y * s0x * jy;
    J->z[idx + nrow    ] += s1y * s0x * jz;

    J->x[idx + nrow + 1] += s1y * s1x * jx;
    J->y[idx + nrow + 1] += s1y * s1x * jy;
    J->z[idx + nrow + 1] += s1y * s1x * jz;

}

/*********************************************************************************************
 
 Sorting
 
 *********************************************************************************************/

/**
 * @brief Sorts particle buffer.
 * 
 * Sorts particles inside the particle buffer according to their cell
 * index to optimize memory cache use. Note: this is a performance
 * optimization only and is not required by the algorithm
 * 
 * @param spec      Particle species
 */
void spec_sort( t_species* spec )
{
    int *idx, *npic;

    int ncell = spec->nx[0]*spec->nx[1];
    
    // Allocate index memory
    idx  = malloc(spec->np*sizeof(int));

    // Allocate temp. array with number of particles in cell
    npic = malloc( ncell * sizeof(int));
    memset( npic, 0, ncell * sizeof(int));

    // Generate sorted index
    int i;
    for (i=0; i<spec->np; i++) {
        idx[i] = spec->part[i].ix + spec->part[i].iy * spec->nx[0];
        npic[idx[i]]++;
    }
    
    int isum = 0, j;
    for (i=0; i<ncell; i++) {
        j = npic[i];
        npic[i] = isum;
        isum += j;
    }
    
    for (i=0; i< spec->np; i++) {
        j = idx[i];
        idx[i] = npic[j]++;
    }
    
    // free temp. array
    free(npic);

    // low mem
    for (i=0; i < spec->np; i++) {
        register t_part tmp;
        register int k;
        
        k = idx[i];
        while ( k > i ) {
            register int t;
            
            tmp = spec->part[k];
            spec->part[k] = spec->part[i];
            spec->part[i] = tmp;
            
            t = idx[k];
            idx[k] = -1;
            k = t;
        }
    }

    free(idx);

}


/*********************************************************************************************
 
 Particle advance
 
 *********************************************************************************************/

/**
 * @brief Interpolates EM fields at particle position
 * 
 * Routine uses linear interpolation and expects quantities to be defined
 * at the lower boundary of the cell
 *
 * @param E     Electric field grid
 * @param B     Magnetic field grid
 * @param part  Particle data
 * @param Ep    E-field interpolated at particle position
 * @param Bp    B-field interpolated at particle position
 */
void interpolate_fld( t_float3_grid2d * E, t_float3_grid2d * B, 
              const t_part* restrict const part, float3* restrict const Ep, float3* restrict const Bp )
{
    float s0x, s0y, s1x, s1y;

    const int nrow = E -> nrow;
    
    const int idx = part->ix + part->iy * nrow;
    
    s0x = 0.5f - part->x;
    s1x = 0.5f + part->x;

    s0y = 0.5f - part->y;
    s1y = 0.5f + part->y;

    Ep -> x = s0y * s0x * E->x[idx           ] + 
              s0y * s1x * E->x[idx        + 1] + 
              s1y * s0x * E->x[idx + nrow    ] + 
              s1y * s1x * E->x[idx + nrow + 1];

    Ep -> y = s0y * s0x * E->y[idx           ] +
              s0y * s1x * E->y[idx        + 1] + 
              s1y * s0x * E->y[idx + nrow    ] + 
              s1y * s1x * E->y[idx + nrow + 1];

    Ep -> z = s0y * s0x * E->z[idx           ] + 
              s0y * s1x * E->z[idx        + 1] + 
              s1y * s0x * E->z[idx + nrow    ] +
              s1y * s1x * E->z[idx + nrow + 1];

    Bp -> x = s0y * s0x * B->x[idx           ] + 
              s0y * s1x * B->x[idx        + 1] + 
              s1y * s0x * B->x[idx + nrow    ] + 
              s1y * s1x * B->x[idx + nrow + 1];

    Bp -> y = s0y * s0x * B->y[idx           ] + 
              s0y * s1x * B->y[idx        + 1] + 
              s1y * s0x * B->y[idx + nrow    ] + 
              s1y * s1x * B->y[idx + nrow + 1];

    Bp -> z = s0y * s0x * B->z[idx           ] + 
              s0y * s1x * B->z[idx        + 1] + 
              s1y * s0x * B->z[idx + nrow    ] +
              s1y * s1x * B->z[idx + nrow + 1];

}

/**
 * @brief Advance Particle species 1 timestep
 * 
 * Particles are advanced in time using a leap-frog method; the velocity
 * advance is done using a relativistic Boris pusher. Current deposition
 * is done at the trajectory mid point; charge deposition is done at the
 * end point.
 * 
 * The routine will also:
 * 1. Calculate total time-centered kinetic energy for the Species
 * 2. Apply boundary conditions
 * 3. Sort particle buffer
 * 
 * @param spec      Particle species
 * @param emf       EM fields
 * @param charge    Charge density
 * @param current   Current density
 */
void spec_advance( t_species* spec, t_emf* emf, t_charge* charge, t_current* current )
{
   
    uint64_t t0;
    t0 = timer_ticks();
    
    const float tem   = 0.5 * spec->dt/spec -> m_q;
    const float dt_dx = spec->dt / spec->dx[0]; 
    const float dt_dy = spec->dt / spec->dx[1]; 

    const int nx0 = spec -> nx[0];
    const int nx1 = spec -> nx[1];

    double energy = 0;

    // Advance particles
    for (int i=0; i<spec->np; i++) {
                
        float3 Ep, Bp;
        float utx, uty, utz;
        float ux, uy, uz, u2;
        float gamma, rg, gtem, otsq;
        
        float x1, y1;
        
        int di, dj;
        float dx, dy; 

        // Load particle momenta
        ux = spec -> part[i].ux;
        uy = spec -> part[i].uy;
        uz = spec -> part[i].uz;

        // interpolate fields
        interpolate_fld( &emf -> E, &emf -> B, &spec -> part[i], &Ep, &Bp );
        
        // advance u using Boris scheme
        Ep.x *= tem;
        Ep.y *= tem;
        Ep.z *= tem;
        
        utx = ux + Ep.x;
        uty = uy + Ep.y;
        utz = uz + Ep.z;

        // Perform first half of the rotation

        // Get time centered gamma
        u2 = utx*utx + uty*uty + utz*utz;
        gamma = sqrtf( 1 + u2 );

        // Accumulate time centered energy
        energy += u2 / ( 1 + gamma );

        gtem = tem / gamma;
        
        Bp.x *= gtem;
        Bp.y *= gtem;
        Bp.z *= gtem;
        otsq = 2.0f / ( 1.0f + Bp.x*Bp.x + Bp.y*Bp.y + Bp.z*Bp.z );

        ux = utx + uty*Bp.z - utz*Bp.y;
        uy = uty + utz*Bp.x - utx*Bp.z;
        uz = utz + utx*Bp.y - uty*Bp.x;
        
        // Perform second half of the rotation
        
        Bp.x *= otsq;
        Bp.y *= otsq;
        Bp.z *= otsq;
        
        utx += uy*Bp.z - uz*Bp.y;
        uty += uz*Bp.x - ux*Bp.z;
        utz += ux*Bp.y - uy*Bp.x;
        
        // Perform second half of electric field acceleration
        ux = utx + Ep.x;
        uy = uty + Ep.y;
        uz = utz + Ep.z;
        
        // Store new momenta
        spec -> part[i].ux = ux;
        spec -> part[i].uy = uy;
        spec -> part[i].uz = uz;
        
        // push particle
        rg = 1.0f / sqrtf(1.0f + ux*ux + uy*uy + uz*uz);
                
        dx = dt_dx * rg * ux;
        dy = dt_dy * rg * uy;

        // Deposit current for the particle (t+1/2)
        deposit_current( &current -> J, &spec -> part[i], spec->q, rg, dx, dy );
        
        x1 = spec -> part[i].x + dx; 
        y1 = spec -> part[i].y + dy;
        
        di = ltrim(x1);
        dj = ltrim(y1);

        x1 -= di;
        y1 -= dj;
        
        
        // Store results
        spec -> part[i].x = x1;
        spec -> part[i].y = y1;
        spec -> part[i].ix += di;
        spec -> part[i].iy += dj;

        // Deposit charge for the particle (t+1)
        deposit_charge( &charge -> rho, &spec -> part[i], spec->q );

        
    }

    // Store energy
    spec -> energy = spec-> q * spec -> m_q * energy * spec -> dx[0] * spec -> dx[1];

    // Advance internal iteration number
    spec -> iter += 1;

    // Use periodic boundaries in both directions
    for (int i=0; i<spec->np; i++) {
        spec -> part[i].ix += (( spec -> part[i].ix < 0 ) ? nx0 : 0 ) - (( spec -> part[i].ix >= nx0 ) ? nx0 : 0);
        spec -> part[i].iy += (( spec -> part[i].iy < 0 ) ? nx1 : 0 ) - (( spec -> part[i].iy >= nx1 ) ? nx1 : 0);
    }
    
    // Sort species at every n_sort time steps
    if ( spec -> n_sort > 0 ) {
        if ( ! (spec -> iter % spec -> n_sort) ) spec_sort( spec );
    }
        
    // Timing info
    _spec_npush += spec -> np;
    _spec_time += timer_interval_seconds( t0, timer_ticks() );
}


/*********************************************************************************************
 
 Charge Deposition
 
 *********************************************************************************************/

/**
 * @brief Deposits particle species charge density
 * 
 * Deposition is done using linear interpolation. Used for diagnostics
 * purpose only.
 * 
 * @param spec      Particle species
 * @param charge    Electric charge density
 */
void spec_deposit_charge( const t_species* spec, float* charge )
{
    int i,j;
    
    // Charge array is expected to have 1 guard cell at the upper boundary
    const int nrow = spec -> nx[0] + 1;
    const float q = spec -> q;
    
    for (i=0; i<spec->np; i++) {
        int idx = spec->part[i].ix + nrow*spec->part[i].iy;

        float s0x = 0.5f - spec->part[i].x;
        float s1x = 0.5f + spec->part[i].x; 

        float s0y = 0.5f - spec->part[i].y;
        float s1y = 0.5f + spec->part[i].y; 
        
        charge[ idx            ] += s0y * s0x * q;
        charge[ idx + 1        ] += s0y * s1x * q;
        charge[ idx     + nrow ] += s1y * s0x * q;
        charge[ idx + 1 + nrow ] += s1y * s1x * q;
    }

    // Correct boundary values

    // x - Periodic boundaries
    for (j = 0; j < spec -> nx[1] + 1; j++) {
        charge[ 0 + j*nrow ] += charge[ spec -> nx[0] + j*nrow ];
    }
    
    // y - Periodic boundaries
    for (i = 0; i < spec->nx[0]+1; i++) {
        charge[ i + 0 ] += charge[ i + spec -> nx[1] * nrow ];
    }

}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/

/**
 * @brief Saves raw particle data information to disk
 * 
 * Saves all particle positions and momenta. Positions are converted to
 * distance from simulation box corner before saving. Data is saved in the
 * "PARTICLES" directory.
 *
 * @param spec 		Particle species
 */
void spec_rep_particles( const t_species *spec )
{
    
    t_zdf_file part_file;

    int i;
    
    const char * quants[] = {
        "x","y",
        "ux","uy","uz"
    };

    const char * qlabels[] = {
        "x", "y",
        "u_x","u_y","u_z"
    };


    const char * qunits[] = {
        "c/\\omega_p", "c/\\omega_p",
        "c","c","c"
    };

    t_zdf_iteration iter = {
        .name = "ITERATION",
        .n = spec->iter,
        .t = spec -> iter * spec -> dt,
        .time_units = "1/\\omega_p"
    };

    // Allocate buffer for positions
    
    t_zdf_part_info info = {
        .name = (char *) spec -> name,
        .nquants = 5,
        .quants = (char **) quants,
        .qlabels = (char **) qlabels,
        .qunits = (char **) qunits,
        .np = spec ->np
    };

    // Create file and add description
    char path[1024];
    snprintf(path, 1024, "PARTICLES/%s", spec -> name );
    zdf_open_part_file( &part_file, &info, &iter, path );

    // Add positions and generalized velocities
    size_t size = ( spec -> np ) * sizeof( float );
    float* data = malloc( size );

    // x1
    for( i = 0; i < spec ->np; i++ )
        data[i] = ( spec -> part[i].ix + (spec -> part[i].x +0.5f ) ) * spec -> dx[0];
    zdf_add_quant_part_file( &part_file, quants[0], data, spec ->np );

    // x2
    for( i = 0; i < spec ->np; i++ )
        data[i] = ( spec -> part[i].iy + (spec -> part[i].y + 0.5f ) ) * spec -> dx[1];
    zdf_add_quant_part_file( &part_file, quants[1], data, spec ->np );

    // ux
    for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].ux;
    zdf_add_quant_part_file( &part_file, quants[2], data, spec ->np );

    // uy
    for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uy;
    zdf_add_quant_part_file( &part_file, quants[3], data, spec ->np );

    // uz
    for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uz;
    zdf_add_quant_part_file( &part_file, quants[4], data, spec ->np );

    free( data );

    zdf_close_file( &part_file );
}	

/**
 * @brief Saves particle species charge density information to disk
 * 
 * @param spec      Particle species
 */
void spec_rep_charge( const t_species *spec )
{
    
    // Add 1 guard cell to the upper boundary
    size_t size = ( spec -> nx[0] + 1 ) * ( spec -> nx[1] + 1 ) * sizeof( float );
    float* charge = malloc( size );
    memset( charge, 0, size );
    
    // Deposit the charge
    spec_deposit_charge( spec, charge );
    
    // Compact the data to save the file (throw away guard cells)
    size = ( spec -> nx[0] ) * ( spec -> nx[1] );
    float * buf = malloc( size * sizeof( float ) );
    
    float * b = buf;
    float * c = charge;
    for( int j = 0; j < spec->nx[1]; j++) {
        for ( int i = 0; i < spec->nx[0]; i++ ) {
            b[i] = c[i];
        }
        b += spec->nx[0];
        c += spec->nx[0] + 1;
    }
    
    free( charge );

    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
        .min = 0.0,
        .max = spec->box[0],
        .name = "x",
        .label = "x",
        .units = "c/\\omega_p"
    };

    axis[1] = (t_zdf_grid_axis) {
        .min = 0.0,
        .max = spec->box[1],
        .name = "y",
        .label = "y",
        .units = "c/\\omega_p"
    };

    char name[128], label[128];
    snprintf(name, 128, "%s-charge", spec -> name);
    snprintf(label, 128, "%s \\rho", spec -> name);

    t_zdf_grid_info info = {
        .ndims = 2,
        .name = name,
        .label = label,
        .units = "n_e",
        .axis  = axis
    };

    info.count[0] = spec->nx[0];
    info.count[1] = spec->nx[1];

    t_zdf_iteration iter = {
        .name = "ITERATION",
        .n = spec->iter,
        .t = spec -> iter * spec -> dt,
        .time_units = "1/\\omega_p"
    };

    char path[1024];
    snprintf(path, 1024, "CHARGE/%s", spec -> name );
    zdf_save_grid( (void *) charge, zdf_float32, &info, &iter, path );	

    free( buf );
}	

/**
 * @brief Gets axis quantity for phasespace density calculation.
 * 
 * Depending on the selected quantity, the appropriate data will be copied
 * or calculated and stored in the output array.
 * 
 * @param spec      Particle species
 * @param i0        Start particle index
 * @param np        Number of particles to process
 * @param quant     Quantity for axis {X1, U1, U2, U3}
 * @param axis      Axis data
 */
void spec_pha_axis( const t_species *spec, int i0, int np, int quant, float *axis )
{
    switch (quant) {
        case X1:
            for (int i = 0; i < np; i++) 
                axis[i] = ( (spec -> part[i0+i].x + 0.5f) + spec -> part[i0+i].ix ) * spec -> dx[0];
            break;
        case X2:
            for (int i = 0; i < np; i++) 
                axis[i] = ( (spec -> part[i0+i].y + 0.5f) + spec -> part[i0+i].iy ) * spec -> dx[1];
            break;
        case U1:
            for (int i = 0; i < np; i++) 
                axis[i] = spec -> part[i0+i].ux;
            break;
        case U2:
            for (int i = 0; i < np; i++) 
                axis[i] = spec -> part[i0+i].uy;
            break;
        case U3:
            for (int i = 0; i < np; i++) 
                axis[i] = spec -> part[i0+i].uz;
            break;
    }
}

/**
 * @brief Phasespace density axis units
 * 
 * @param quant     Quantity for axis
 * @return          Units for phasespace density axis
 */
const char * spec_pha_axis_units( int quant ) {
    switch (quant) {
        case X1:
            return("c/\\omega_p");
            break;
        case U1:
        case U2:
        case U3:
            return("m_e c");
    }
    return("");
}

/**
 * @brief Deposit 2D phasespace density.
 * 
 * @param spec      Particle species
 * @param rep_type  Type of phasespace, use the PHASESPACE macro to define
 * @param pha_nx    Number of grid points in the phasespace density grid
 * @param pha_range Physical range of each of the phasespace axis
 * @param buf       Phasespace density grid
 */
void spec_deposit_pha( const t_species *spec, const int rep_type,
              const int pha_nx[], const float pha_range[][2], float* restrict buf )
{
    const int BUF_SIZE = 1024;
    float pha_x1[BUF_SIZE], pha_x2[BUF_SIZE];


    const int nrow = pha_nx[0];

    const int quant1 = rep_type & 0x000F;
    const int quant2 = (rep_type & 0x00F0)>>4;

    const float x1min = pha_range[0][0];
    const float x2min = pha_range[1][0];

    const float rdx1 = pha_nx[0] / ( pha_range[0][1] - pha_range[0][0] );
    const float rdx2 = pha_nx[1] / ( pha_range[1][1] - pha_range[1][0] );

    for ( int i = 0; i<spec->np; i+=BUF_SIZE ) {
        int np = ( i + BUF_SIZE > spec->np )? spec->np - i : BUF_SIZE;

        spec_pha_axis( spec, i, np, quant1, pha_x1 );
        spec_pha_axis( spec, i, np, quant2, pha_x2 );

        for ( int k = 0; k < np; k++ ) {

            float nx1 = ( pha_x1[k] - x1min ) * rdx1;
            float nx2 = ( pha_x2[k] - x2min ) * rdx2;

            int i1 = (int)(nx1 + 0.5f);
            int i2 = (int)(nx2 + 0.5f);

            float w1 = nx1 - i1 + 0.5f;
            float w2 = nx2 - i2 + 0.5f;

            int idx = i1 + nrow*i2;

            if ( i2 >= 0 && i2 < pha_nx[1] ) {

                if (i1 >= 0 && i1 < pha_nx[0]) {
                    buf[ idx ] += (1.0f-w1)*(1.0f-w2)*spec->q;
                }

                if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
                    buf[ idx + 1 ] += w1*(1.0f-w2)*spec->q;
                }
            }

            idx += nrow;
            if ( i2+1 >= 0 && i2+1 < pha_nx[1] ) {

                if (i1 >= 0 && i1 < pha_nx[0]) {
                    buf[ idx ] += (1.0f-w1)*w2*spec->q;
                }

                if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
                    buf[ idx + 1 ] += w1*w2*spec->q;
                }
            }

        }

    }
}

/**
 * @brief Saves particle species phasespace density information to disk
 * 
 * @param spec      Particle species
 * @param rep_type  Type of phasespace, use the PHASESPACE macro to define
 * @param pha_nx    Number of grid points in the phasespace density grid
 * @param pha_range Physical range of each of the phasespace axis
 */
void spec_rep_pha( const t_species *spec, const int rep_type,
              const int pha_nx[], const float pha_range[][2] )
{
    // Allocate phasespace buffer
    float* restrict buf = malloc( pha_nx[0] * pha_nx[1] * sizeof( float ));
    memset( buf, 0, pha_nx[0] * pha_nx[1] * sizeof( float ));

    // Deposit the phasespace
    spec_deposit_pha( spec, rep_type, pha_nx, pha_range, buf );

    // save the data in hdf5 format
    int quant1 = rep_type & 0x000F;
    int quant2 = (rep_type & 0x00F0)>>4;

    const char * pha_ax1_units = spec_pha_axis_units(quant1);
    const char * pha_ax2_units = spec_pha_axis_units(quant2);

    char const * const pha_ax_name[] = {"x1","x2","x3","u1","u2","u3"};
    char const * const pha_ax_label[] = {"x","y","z","u_x","u_y","u_z"};

    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
        .min = pha_range[0][0],
        .max = pha_range[0][1],
        .name  = (char *) pha_ax_name[ quant1 - 1 ],
        .label = (char *) pha_ax_name[ quant1 - 1 ],
        .units = (char *) pha_ax1_units
    };

    axis[1] = (t_zdf_grid_axis) {
        .min = pha_range[1][0],
        .max = pha_range[1][1],
        .name  = (char *) pha_ax_name[ quant2 - 1 ],
        .label = (char *) pha_ax_name[ quant2 - 1 ],
        .units = (char *) pha_ax2_units
    };

    char pha_name[64], pha_label[64];
    
    snprintf( pha_name, 64,"%s-%s%s", spec -> name, 
        pha_ax_name[quant1-1], pha_ax_name[quant2-1] );
    
    snprintf( pha_label, 64,"%s %s-%s", spec -> name,
        pha_ax_label[quant1-1], pha_ax_label[quant2-1] );

    t_zdf_grid_info info = {
        .ndims = 2,
        .name = pha_name,
        .label = pha_label,
        .units = "a.u.",
        .axis  = axis
    };

    info.count[0] = pha_nx[0];
    info.count[1] = pha_nx[1];

    t_zdf_iteration iter = {
        .name = "ITERATION",
        .n = spec->iter,
        .t = spec -> iter * spec -> dt,
        .time_units = "1/\\omega_p"
    };

    char path[1024];
    snprintf(path, 1024, "PHASESPACE/%s", spec -> name );
    zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, path );

    // Free temp. buffer
    free( buf );

}

/**
 * @brief Saves particle species diagnostic information to disk
 * 
 * @param spec      Particle species
 * @param rep_type  Type of diagnostic information {CHARGE, PHASESPACE(a,b), PARTICLES}
 * @param pha_nx    Number of grid points in the phasespace density grid, set to NULL for
 *                  diagnostics other then phasespace density
 * @param pha_range Physical range of each of the phasespace axis, set to NULL for
 *                  diagnostics other then phasespace density
 */
void spec_report( const t_species *spec, const int rep_type, 
                  const int pha_nx[], const float pha_range[][2] )
{
    
    switch (rep_type & 0xF000) {
        case CHARGE:
            spec_rep_charge( spec );
            break;

        case PHA:
            spec_rep_pha( spec, rep_type, pha_nx, pha_range );
            break;

        case PARTICLES:
            spec_rep_particles( spec );
            break;
    }
}
