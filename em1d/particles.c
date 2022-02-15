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

#include "particles.h"

#include "random.h"
#include "emf.h"
#include "current.h"

#include "zdf.h"
#include "timer.h"

static double _spec_time = 0.0;
static uint64_t _spec_npush = 0;

void spec_sort( t_species *spec );
void spec_move_window( t_species *spec );


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
    return (_spec_npush > 0 )? _spec_time / _spec_npush: -1.0;
}

/*********************************************************************************************
 Initialization
 *********************************************************************************************/

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
    /**
     * Version 1 momentum initialization
     */

    // Initialize thermal component
    for (int i = start; i <= end; i++) {
        spec->part[i].ux = spec -> uth[0] * rand_norm();
        spec->part[i].uy = spec -> uth[1] * rand_norm();
        spec->part[i].uz = spec -> uth[2] * rand_norm();
    }

    // Calculate net momentum in each cell
    float3 * restrict net_u = (float3 *) malloc( spec->nx * sizeof(float3));
    int * restrict    npc   = (int *) malloc( spec->nx * sizeof(int));

    // Zero momentum grids
    memset(net_u, 0, spec->nx * sizeof(float3) );
    memset(npc, 0, (spec->nx) * sizeof(int) );

    // Accumulate momentum in each cell
    for (int i = start; i <= end; i++) {
        const int idx  = spec -> part[i].ix;

        net_u[ idx ].x += spec->part[i].ux;
        net_u[ idx ].y += spec->part[i].uy;
        net_u[ idx ].z += spec->part[i].uz;

        npc[ idx ] += 1;
    }

    // Normalize to the number of particles in each cell to get the
    // average momentum in each cell
    for(int i =0; i< spec->nx; i++ ) {
        const float norm = (npc[ i ] > 0) ? 1.0f/npc[i] : 0;

        net_u[ i ].x *= norm;
        net_u[ i ].y *= norm;
        net_u[ i ].z *= norm;
    }

    // Subtract average momentum and add fluid component
    for (int i = start; i <= end; i++) {
        const int idx  = spec -> part[i].ix;

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
 * @brief Gets number of particles to be injected.
 *
 * Calculates the number of particles to be injected in the specified range according
 * to the specified density profile. The returned value is not exact but it is
 * guaranteed to be larger than the actual number of particles to be injected
 *
 * @param spec      Particle species
 * @param range     Range of cells in which to inject
 * @return          Number of particles to be injected
 */
int spec_np_inj( t_species* spec, const int range[] )
{
    int np_inj;

    switch ( spec -> density.type ) {
    case STEP: // Step like density profile
        {
            int i0 = spec -> density.start / spec -> dx - spec -> n_move;

            if ( i0 > range[1] ) {
                np_inj = 0;
            } else {
                if ( i0 < range[0] ) i0 = range[0];
                np_inj = ( range[1] - i0 + 1 ) * spec -> ppc;
            }
        }
        break;

    case SLAB: // Slab like density profile
        {
            int i0 = spec -> density.start / spec -> dx - spec -> n_move;
            int i1 = spec -> density.end / spec -> dx - spec -> n_move;

            if ( (i0 > range[1]) || (i1 < range[0]) ) {
                np_inj = 0;
            } else {
                if ( i0 < range[0] ) i0 = range[0];
                if ( i1 > range[1] ) i1 = range[1];
                np_inj = ( i1 - i0 + 1 ) * spec -> ppc;
            }
        }
        break;

    case RAMP: // ramp density profile
        {
            // Ramp start / finish
            float x0 = spec -> density.start;
            float x1 = spec -> density.end;

            // Injection start / finish positions (in simulation units)
            float a = (range[0] + spec -> n_move) * spec->dx;
            float b = (range[1] + 1 + spec -> n_move) * spec->dx;

            // If outside of ramp (or invalid ramp) return 0
            if ( (x1 <= x0) || (a > x1) || (b < x0) ) {
                np_inj = 0;
            } else {
                // limit integration boundaries to ramp start/end
                if ( a < x0 ) a = x0;
                if ( b > x1 ) b = x1;

                // Get total injected charge
                float n0 = spec -> density.ramp[0];
                float n1 = spec -> density.ramp[1];
                float q = (b-a)*( n0 + 0.5 * (a+b-2*x0)*(n1-n0)/(x1-x0));

                // Get corresponding number of simulation particles
                np_inj = q * spec -> ppc / spec -> dx;

            }
        }
        break;

    case CUSTOM: // custom density profile
        {
            // Integrate total charge
            double q = 0.5 * ( (*spec -> density.custom)((range[0] + spec -> n_move) * spec->dx,
                                                         spec -> density.custom_data) +
                               (*spec -> density.custom)((range[1] + 1 + spec -> n_move) * spec->dx,
                                                            spec -> density.custom_data) );

            for( int i = range[0]+1; i <= range[1]; i++) {
                q += (*spec -> density.custom)((i + spec -> n_move) * spec->dx,
                                               spec -> density.custom_data);
            }

            // Get corresponding number of simulation particles, rounding up
            np_inj = ceil(q * spec -> ppc);
        }
        break;
    
    case EMPTY: // Empty profile
        np_inj = 0;
        break;

    default: // Uniform density
        np_inj = ( range[1] - range[0] + 1 ) * spec -> ppc;
    }

    // printf("Predicts injecting %d particles\n", np_inj);
    return np_inj;

}

/**
 * @brief Sets initial position of particles according to density profile
 * 
 * @param spec      Particle species
 * @param range     Range of cells in which to inject
 */
void spec_set_x( t_species* spec, const int range[] )
{

    int i, k, ip;
    float start, end;

    // Calculate particle positions inside the cell
    const int npc = spec->ppc;

    float poscell[npc];

    for (i=0; i<spec->ppc; i++) {
        poscell[i]   = ( i + 0.5 ) / npc;
    }

    ip = spec -> np;

    // Set position of particles in the specified grid range according to the density profile
    switch ( spec -> density.type ) {
    case STEP: // Step like density profile

        // Get edge position normalized to cell size;
        start = spec -> density.start / spec -> dx - spec -> n_move;

        for (i = range[0]; i <= range[1]; i++) {

            for (k=0; k<npc; k++) {
                if ( i + poscell[k] > start ) {
                    spec->part[ip].ix = i;
                    spec->part[ip].x = poscell[k];
                    ip++;
                }
            }
        }

        // printf("Injected %d particles with step injection \n", ip - spec -> np );
        break;

    case SLAB: // Slab like density profile

        // Get edge position normalized to cell size;
        start = spec -> density.start / spec -> dx - spec -> n_move;
        end   = spec -> density.end / spec -> dx - spec -> n_move;

        for (i = range[0]; i <= range[1]; i++) {

            for (k=0; k<npc; k++) {
                if ( i + poscell[k] > start &&  i + poscell[k] < end ) {
                    spec->part[ip].ix = i;
                    spec->part[ip].x = poscell[k];
                    ip++;
                }
            }
        }

        // printf("Injected %d particles with slab injection \n", ip - spec -> np );
        break;

    case RAMP: // ramp like density profile

        {
            // Ramp start/finish in cell units
            double r0 = spec -> density.start / spec -> dx;
            double r1 = spec -> density.end / spec -> dx;

            // If outside ramp return
            if (((range[0] + spec -> n_move) > r1 ) ||
                ((range[1] + spec -> n_move) < r0 )) break;

            double n0 = spec -> density.ramp[0];
            double n1 = spec -> density.ramp[1];

            // Only consider the ramp for x > 0
            if ( r0 < 0 ) {
                n0 += - r0 * (n1-n0) / (r1-r0);
                r0 = 0;
            }

            // Charge per simulation particle
            double cpp = 1.0 / spec->ppc;

            for( k = spec -> density.total_np_inj; ; k++ ) {
                // Desired cumulative density, normalized to the [0,1] interval
                double Rs = (k+0.5) * cpp / (r1 - r0);

                // Position normalized to the [0,1] interval
                // double pos = (-a + sqrt( a*a + 2 * b * Rs ))/b;
                double pos = 2 * Rs / (sqrt( n0*n0 + 2 * (n1-n0) * Rs ) + n0);

                // If outside of ramp interval we are done
                if ( pos > 1 ) break;

                // Position in simulation cell units
                pos = r0 + (r1-r0) * pos;

                // Injection cell
                int ix = pos;

                // (*debug*) This must never happen
                if ( ix - spec -> n_move < range[0] ) {
                    fprintf(stderr, "(*error*) attempting to inject outside of valid range.\n");
                    break;
                }

                // If outside injection range we are done
                if ( ix - spec -> n_move > range[1] ) break;

                // Inject particle
                spec->part[ip].ix = ix - spec -> n_move;
                spec->part[ip].x = pos - ix;
                ip++;

            }
        }

        // printf("Injected %d particles with ramp injection \n", ip - spec -> np );
        break;

    case CUSTOM: // custom density profile

        {

            const double dx = spec -> dx;

            // Charge per simulation particle
            const double cpp = 1.0 / spec->ppc;

            // Injected particles
            k = spec -> density.total_np_inj;

            int ix = range[0];

            // Density on cell edges
            double n0;
            double n1 = (*spec -> density.custom)((ix + spec -> n_move) * dx, spec -> density.custom_data);

            // Accumulated density on cell edges
            double d0;
            double d1 = spec -> density.custom_q_inj;

            double Rs;

            while( ix <= range[1] ){

                // Get density on the edges of current cell
                n0 = n1;
                n1 = (*spec -> density.custom)((ix + 1 + spec -> n_move)*dx, spec -> density.custom_data);

                // Get cumulative density on the edges of current cell
                d0 = d1;
                d1 += 0.5 * (n0+n1);

                while( ( Rs =  (k+0.5) * cpp ) < d1 ) {

                    // Quadratic formula
                    // double pos = (-n0 + sqrt( n0*n0 + 2 * (n1-n0) * (Rs-d0) ))/(n1-n0);

                    // This version avoids a division by 0 if n1 = n0
                    double pos = 2 * (Rs-d0) /( sqrt( n0*n0 + 2 * (n1-n0) * (Rs-d0) ) + n0 );

                    spec->part[ip].ix = ix;
                    spec->part[ip].x = pos;
                    ip++;

                    k++;
                }

                // Move to next cell
                ix++;
            }

            spec -> density.custom_q_inj = d1;
        }

        // printf("Injected %d particles with custom injection \n", ip - spec -> np );
        break;
    
    case EMPTY: // Empty profile
        break;

    default: // Uniform density
        for (i = range[0]; i <= range[1]; i++) {

            for (k=0; k<npc; k++) {
                spec->part[ip].ix = i;
                spec->part[ip].x = poscell[k];
                ip++;
            }
        }
        // printf("Injected %d particles with uniform injection \n", ip - spec -> np );

    }

    // Update total number of injected particles
    spec -> density.total_np_inj += ip - spec -> np;

    // Update number of particles in buffer
    spec -> np = ip;

}

/**
 * @brief Grows particle buffer to specified size.
 * 
 * If the new size is smaller than the previous size the buffer size is not changed
 * and the function returns silently.
 * 
 * @param spec  Particle species
 * @param size  New buffer size (will be rounded up to next multiple of 1024)
 **/
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
 * @param   range   Grid range [ix0, ix1] where to inject particles
 **/
void spec_inject_particles( t_species* spec, const int range[] )
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
 * @param ppc       Reference number of particles per cell
 * @param ufl       Initial fluid momentum of particles, may be set to NULL 
 * @param uth       Initial thermal momentum of particles, may be set to NULL
 * @param nx        Number of grid points
 * @param box       Simulation box size in simulation units
 * @param dt        Simulation time step, in simulation units
 * @param density   Density profile for particle injection, may be set to NULL
 **/
void spec_new( t_species* spec, char name[], const float m_q, const int ppc,
              const float *ufl, const float * uth,
              const int nx, float box, const float dt, t_density* density )
{

    int i, npc;

    // Species name
    strncpy( spec -> name, name, MAX_SPNAME_LEN );

    spec->nx = nx;
    spec->ppc = ppc;
    npc = ppc;

    spec->box = box;
    spec->dx = box / nx;

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
    spec -> density.total_np_inj = 0;
    spec -> density.custom_q_inj = 0.;

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

    // Reset moving window information
    spec -> moving_window = 0;
    spec -> n_move = 0;

    // Inject initial particle distribution
    spec -> np = 0;

    const int range[2] = {0, nx-1};

    spec_inject_particles( spec, range );

    // Set default sorting frequency
    spec -> n_sort = 16;

    // Default to periodic boundary condtions
    spec -> bc_type = PART_BC_PERIODIC;

}

/**
 * @brief Move simulation window
 * 
 * When using a moving simulation window checks if a window move is due
 * at the current iteration and if so shifts left the particle cell indices
 * 
 * @param spec      Particle species
 */
void spec_move_window( t_species *spec ){

    if ((spec->iter * spec->dt ) > (spec->dx * (spec->n_move + 1)))  {

        // shift all particles left
        // particles leaving the box will be removed later
        int i;
        for( i = 0; i < spec->np; i++ ) {
            spec->part[i].ix--;
        }

        // Increase moving window counter
        spec -> n_move++;

        // Inject particles in the right edge of the simulation box
        const int range[2] = {spec->nx-1,spec->nx-1};
        spec_inject_particles( spec, range );

    }

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

 Cuurent deposition

 *********************************************************************************************/

/**
 * @brief Deposit single particle current using Esirkepov method
 * 
 * @param ix0       Initial cell index of particle
 * @param di        Number of cells moved {-1,0,1}
 * @param x0        Initial position of particle inside cell
 * @param x1        Final position of particle inside cell
 * @param qnx       Normalization for x current (q * cell size / dt)
 * @param qvy       Y current ( q * vy )
 * @param qvz       Z current ( q * vz )
 * @param current   Electric current density
 */
void dep_current_esk( int ix0, int di,
                        float x0, float x1,
                        float qnx, float qvy, float qvz,
                        t_current *current )
{

    float S0x[4], S1x[4], DSx[4];
    float Wx[4], Wy[4], Wz[4];

    S0x[0] = 0.0f;
    S0x[1] = 1.0f - x0;
    S0x[2] = x0;
    S0x[3] = 0.0f;

    for (int i=0; i<4; i++) {
        S1x[i] = 0.0f;
    }

    S1x[ 1 + di ] = 1.0f - x1;
    S1x[ 2 + di ] = x1;

    for (int i=0; i<4; i++) {
        DSx[i] = S1x[i] - S0x[i];
    }

    for (int i=0; i<4; i++) {
        Wx[i] = qnx * DSx[i];
        Wy[i] = qvy * (S0x[i] + DSx[i]/2.0f);
        Wz[i] = qvz * (S0x[i] + DSx[i]/2.0f);
    }

    float3* restrict const J = current -> J;
    // jx
    float c;

    c = - Wx[0];
    J[ ix0 - 1 ].x += c;
    for (int i=1; i<4; i++) {
        c -=  Wx[i];
        J[ ix0 + i ].x += c;
    }

    // jy, jz
    for (int i=0; i<4; i++) {
        J[ ix0 + i - 1 ].y += Wy[ i ];
        J[ ix0 + i - 1 ].z += Wz[ i ];
    }

}

/**
 * @brief Deposit single particle current using zamb method
 * 
 * @param ix0       Initial cell index of particle
 * @param di        Number of cells moved {-1,0,1}
 * @param x0        Initial position of particle inside cell
 * @param dx        Particle motion normalized to cell size
 * @param qnx       Normalization for x current (q * cell size / dt)
 * @param qvy       Y current ( q * vy )
 * @param qvz       Z current ( q * vz )
 * @param current   Electric current density
 */
void dep_current_zamb( int ix0, int di,
                        float x0, float dx,
                        float qnx, float qvy, float qvz,
                        t_current *current )
{
    // Split the particle trajectory
    typedef struct {
        float x0, x1, dx, qvy, qvz;
        int ix;
    } t_vp;

    t_vp vp[3];
    int vnp = 1;

    // split
    vp[0].x0 = x0;
    vp[0].dx = dx;

    vp[0].x1 = x0+dx;

    vp[0].qvy = qvy/2.0;
    vp[0].qvz = qvz/2.0;

    vp[0].ix = ix0;

    // x split
    if ( di != 0 ) {

        //int ib = ( di+1 )>>1;
        int ib = ( di == 1 );

        float delta = (x0+dx-ib)/dx;

        // Add new particle
        vp[1].x0 = 1-ib;
        vp[1].x1 = (x0 + dx) - di;
        vp[1].dx = dx*delta;
        vp[1].ix = ix0 + di;

        vp[1].qvy = vp[0].qvy*delta;
        vp[1].qvz = vp[0].qvz*delta;

        // Correct previous particle
        vp[0].x1 = ib;
        vp[0].dx *= (1.0f-delta);

        vp[0].qvy *= (1.0f-delta);
        vp[0].qvz *= (1.0f-delta);

        vnp++;

    }

    // Deposit virtual particle currents
    float3* restrict const J = current -> J;

    for (int k = 0; k < vnp; k++) {
        float S0x[2], S1x[2];

        S0x[0] = 1.0f - vp[k].x0;
        S0x[1] = vp[k].x0;

        S1x[0] = 1.0f - vp[k].x1;
        S1x[1] = vp[k].x1;

        J[ vp[k].ix     ].x += qnx * vp[k].dx;
        J[ vp[k].ix     ].y += vp[k].qvy * (S0x[0]+S1x[0]+(S0x[0]-S1x[0])/2.0f);
        J[ vp[k].ix + 1 ].y += vp[k].qvy * (S0x[1]+S1x[1]+(S0x[1]-S1x[1])/2.0f);
        J[ vp[k].ix     ].z += vp[k].qvz * (S0x[0]+S1x[0]+(S0x[0]-S1x[0])/2.0f);
        J[ vp[k].ix  +1 ].z += vp[k].qvz * (S0x[1]+S1x[1]+(S0x[1]-S1x[1])/2.0f);
    }

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

    const int ncell = spec->nx;

    // Allocate index memory
    int * restrict idx  = (int *) malloc(spec->np*sizeof(int));

    // Allocate temp. array with number of particles in cell
    int * restrict npic = (int *) malloc( ncell * sizeof(int));
    memset( npic, 0, ncell * sizeof(int));

    // Generate sorted index
    for (int i=0; i<spec->np; i++) {
        idx[i] = spec->part[i].ix;
        npic[idx[i]]++;
    }

    int isum = 0;
    for (int i=0; i<ncell; i++) {
        int j = npic[i];
        npic[i] = isum;
        isum += j;
    }

    for (int i=0; i< spec->np; i++) {
        int j = idx[i];
        idx[i] = npic[j]++;
    }

    // free temp. array
    free(npic);

    // low mem
    for (int i=0; i < spec->np; i++) {
        int k = idx[i];
        while ( k > i ) {
            t_part tmp = spec->part[k];
            spec->part[k] = spec->part[i];
            spec->part[i] = tmp;

            int t = idx[k];
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
 * Routine uses linear interpolation and accounts for a staggered (Yee)
 * mesh, with the charge at the lower corner of the cell
 * 
 * @param E     Electric field grid
 * @param B     Magnetic field grid
 * @param part  Particle data
 * @param Ep    E-field interpolated at particle position
 * @param Bp    B-field interpolated at particle position
 */
void interpolate_fld( const float3* restrict const E, const float3* restrict const B,
              const t_part* restrict const part, float3* restrict const Ep, float3* restrict const Bp )
{
    int i, ih;
    float w1, w1h;

    i = part->ix;

    w1 = part->x;
    ih = (w1 <0.5f)? -1 : 0;
    w1h = w1 + ((w1 <0.5f)?0.5f:-0.5f);

    ih += i;

    Ep->x = E[ih].x * (1.0f - w1h) + E[ih+1].x * w1h;
    Ep->y = E[i ].y * (1.0f -  w1) + E[i+1 ].y * w1;
    Ep->z = E[i ].z * (1.0f -  w1) + E[i+1 ].z * w1;

    Bp->x = B[i ].x * (1.0f  - w1) + B[i+1 ].x * w1;
    Bp->y = B[ih].y * (1.0f - w1h) + B[ih+1].y * w1h;
    Bp->z = B[ih].z * (1.0f - w1h) + B[ih+1].z * w1h;

}

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
    return ( x >= 1.0f ) - ( x < 0.0f );
}

/**
 * @brief Advance Particle species 1 timestep
 * 
 * Particles are advanced in time using a leap-frog method; the velocity
 * advance is done using a relativistic Boris pusher. Particle motion is
 * used to deposit electric current on the grid using an exact charge
 * conservation method.
 * 
 * The routine will also:
 * 1. Calculate total time-centered kinetic energy for the Species
 * 2. Apply boundary conditions
 * 3. Move simulation window
 * 4. Sort particle buffer
 * 
 * @param spec      Particle species
 * @param emf       EM fields
 * @param current   Current density
 */
void spec_advance( t_species* spec, t_emf* emf, t_current* current )
{

    uint64_t t0;
    t0 = timer_ticks();

    const float tem   = 0.5 * spec->dt/spec -> m_q;
    const float dt_dx = spec->dt / spec->dx;

    // Auxiliary values for current deposition
    const float qnx = spec -> q *  spec->dx / spec->dt;

    const int nx0 = spec -> nx;

    double energy = 0;

    // Advance particles
    for (int i=0; i<spec->np; i++) {

        float3 Ep, Bp;
        float utx, uty, utz;
        float ux, uy, uz, u2;
        float gamma, rg, gtem, otsq;

        float x1;

        int di;
        float dx;

        // Load particle momenta
        ux = spec -> part[i].ux;
        uy = spec -> part[i].uy;
        uz = spec -> part[i].uz;

        // interpolate fields
        interpolate_fld( emf -> E_part, emf -> B_part, &spec -> part[i], &Ep, &Bp );
        // Ep.x = Ep.y = Ep.z = Bp.x = Bp.y = Bp.z = 0;

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

        x1 = spec -> part[i].x + dx;

        di = ltrim(x1);

        x1 -= di;

        float qvy = spec->q * uy * rg;
        float qvz = spec->q * uz * rg;

        // deposit current using Eskirepov method
        // dep_current_esk( spec -> part[i].ix, di,
        // 				 spec -> part[i].x, x1,
        // 				 qnx, qvy, qvz,
        // 				 current );

        dep_current_zamb( spec -> part[i].ix, di,
                         spec -> part[i].x, dx,
                         qnx, qvy, qvz,
                         current );

        // Store results
        spec -> part[i].x = x1;
        spec -> part[i].ix += di;

    }

    // Store energy
    spec -> energy = spec-> q * spec -> m_q * energy * spec -> dx;

    // Advance internal iteration number
    spec -> iter += 1;

    // Check for particles leaving the box
    if ( spec -> moving_window || spec -> bc_type == PART_BC_OPEN ){

        // Move simulation window if needed
        if (spec -> moving_window )	spec_move_window( spec );

        // Use absorbing boundaries along x
        int i = 0;
        while ( i < spec -> np ) {
            if (( spec -> part[i].ix < 0 ) || ( spec -> part[i].ix >= nx0 )) {
                spec -> part[i] = spec -> part[ -- spec -> np ];
                continue;
            }
            i++;
        }

    } else {
        // Use periodic boundaries in x
        for (int i=0; i<spec->np; i++) {
            spec -> part[i].ix += (( spec -> part[i].ix < 0 ) ? nx0 : 0 ) - (( spec -> part[i].ix >= nx0 ) ? nx0 : 0);
        }
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
 * Deposition is done using linear interpolation, charge grid is expected
 * to have 1 guard cell at the upper boundary.
 * 
 * Used for diagnostics purposes only.
 * 
 * @param spec      Particle species
 * @param charge    Electric charge density
 */
void spec_deposit_charge( const t_species* spec, float* charge )
{
    const float q = spec -> q;

    // Charge array is expected to have 1 guard cell at the upper boundary
    for (int i=0; i<spec->np; i++) {
        int idx = spec->part[i].ix;
        float w1 = spec->part[i].x;

        charge[ idx            ] += ( 1.0f - w1 ) * q;
        charge[ idx + 1        ] += (        w1 ) * q;
    }

    // Correct boundary values

    // x
    if ( ! spec -> moving_window ){
        charge[ 0 ] += charge[ spec -> nx ];
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
 * "PARTICLES/sp_name" directory.
 *
 * @param spec 		Particle species
 */
void spec_rep_particles( const t_species *spec )
{

    t_zdf_file part_file;

    int i;

    const char * quants[] = {
        "x",
        "ux","uy","uz"
    };

    const char * qlabels[] = {
        "x",
        "u_x","u_y","u_z"
    };

    const char * qunits[] = {
        "c/\\omega_p",
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
        .label = (char *) spec -> name,
        .nquants = 4,
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

    // x
    for( i = 0; i < spec ->np; i++ )
        data[i] = (spec -> n_move + spec -> part[i].ix + spec -> part[i].x ) * spec -> dx;
    zdf_add_quant_part_file( &part_file, quants[0], data, spec ->np );

    // ux
    for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].ux;
    zdf_add_quant_part_file( &part_file, quants[1], data, spec ->np );

    // uy
    for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uy;
    zdf_add_quant_part_file( &part_file, quants[2], data, spec ->np );

    // uz
    for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uz;
    zdf_add_quant_part_file( &part_file, quants[3], data, spec ->np );

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
    size_t size = ( spec -> nx + 1 ) * sizeof( float );
    float *charge = malloc( size );
    memset( charge, 0, size );

    // Deposit the charge
    spec_deposit_charge( spec, charge );

    // Compact the data to save the file (throw away guard cells)
    float buffer[ spec -> nx ];

    for ( int i = 0; i < spec->nx; i++ ) {
        buffer[i] = charge[i];
    }

    free( charge );

    // Set grid boundaries accounting for moving window
    t_zdf_grid_axis axis = {
        .min = spec -> n_move * spec -> dx,
        .max = spec->box + spec -> n_move * spec -> dx,
        .name = "x",
        .label = "x",
        .units = "c/\\omega_p"
    };

    char name[128], label[128];
    snprintf(name, 128, "%s-charge", spec -> name);
    snprintf(label, 128, "%s \\rho", spec -> name);

    t_zdf_grid_info info = {
        .ndims = 1,
        .name = name,
        .label = label,
        .units = "n_e",
        .axis  = &axis
    };

    info.count[0] = spec->nx;

    t_zdf_iteration iter = {
        .name = "ITERATION",
        .n = spec->iter,
        .t = spec -> iter * spec -> dt,
        .time_units = "1/\\omega_p"
    };

    char path[1024];
    snprintf(path, 1024, "CHARGE/%s", spec -> name );
    zdf_save_grid( (void *) buffer, zdf_float32, &info, &iter, path );
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
void spec_pha_axis( const t_species *spec, int i0, int np, int quant,
    float * restrict axis )
{
    switch (quant) {
        case X1:
            for (int i = 0; i < np; i++)
                axis[i] = ( spec -> part[i0+i].x + spec -> part[i0+i].ix ) * spec -> dx;
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
        .label = (char *) pha_ax_label[ quant1 - 1 ],
        .units = (char *) pha_ax1_units
    };

    axis[1] = (t_zdf_grid_axis) {
        .min = pha_range[1][0],
        .max = pha_range[1][1],
        .name  = (char *) pha_ax_name[ quant2 - 1 ],
        .label = (char *) pha_ax_label[ quant2 - 1 ],
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
