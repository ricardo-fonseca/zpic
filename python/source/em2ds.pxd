#cython: language_level=3

from libc.stdint cimport uint32_t

cdef extern from "../../em2ds/random.h":
	void set_rand_seed( uint32_t m_z_, uint32_t m_w_ )


#########################################################################################
# Grids
#
cdef extern from "../../em2ds/grid2d.h":
	ctypedef struct t_scalar_grid2d:
		float *s
		float *buffer
		int nx[2]
		int nrow
		int gc[2][2]

	ctypedef struct t_cscalar_grid2d:
		float complex *s
		float complex *buffer
		int nx[2]
		int nrow
		int gc[2][2]

	ctypedef struct float3:
		float x
		float y
		float z

	ctypedef struct t_float3_grid2d:
		float *x
		float *y
		float *z
		float *buffer
		int nx[2]
		int nrow
		int gc[2][2]

	ctypedef struct t_cfloat3_grid2d:
		float complex *x
		float complex *y
		float complex *z
		float complex *buffer
		int nx[2]
		int nrow
		int gc[2][2]

#########################################################################################
# FFT
#
cdef extern from "../../em2ds/fft.h":
	cdef enum:
		MAX_FACTORS

	cdef enum fft_direction:
		FFT_FORWARD, FFT_BACKWARD

	ctypedef struct t_fft_factors:
		unsigned int p
		unsigned int n

	ctypedef struct t_fft_cfg:
		unsigned int n
		fft_direction direction
		t_fft_factors factors[MAX_FACTORS]
		float complex *phase

	cdef enum fftr_storage:
		FFTR_CCS, FFTR_PERM

	ctypedef struct t_fftr_cfg:
		fftr_storage storage
		t_fft_cfg cfg
		float complex *phase

	ctypedef struct t_fftr2d_cfg:
		int nx
		int ny
		int stridey
		t_fftr_cfg cfgx
		t_fft_cfg cfgy
		float complex *phase

#########################################################################################
# Particles
#
cdef extern from "../../em2ds/particles.h":
	cdef enum:
		MAX_SPNAME_LEN

	ctypedef struct t_part:
		int ix
		int iy
		float x
		float y
		float ux
		float uy
		float uz

	cdef enum density_type:
		UNIFORM, EMPTY, STEP, SLAB

	ctypedef struct t_density:
		float n
		density_type type
		float start
		float end

	ctypedef struct t_species:
		char name[MAX_SPNAME_LEN]
		t_part *part
		int np
		int np_max
		float m_q
		float q
		double energy
		int ppc[2]
		t_density density
		float ufl[3]
		float uth[3]
		int nx[2]
		float dx[2]
		float box[2]
		float dt
		int iter
		int n_sort

	void spec_new( t_species* spec, char name[], const float m_q, const int ppc[],
				  const float ufl[], const float uth[],
				  const int nx[], float box[], const float dt, t_density* density )

	void spec_grow_buffer( t_species* spec, const int size )

	cdef int CHARGE
	cdef int PHA
	cdef int PARTICLES
	cdef int X1
	cdef int X2
	cdef int U1
	cdef int U2
	cdef int U3

	int PHASESPACE( int a, int b )

	void spec_report( const t_species *spec, const int rep_type,
				  const int pha_nx[], const float pha_range[][2] )

	void spec_deposit_charge( const t_species* spec, float* charge )

	void spec_deposit_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2], float* buf )

#########################################################################################
# Charge
#

cdef extern from "../../em2ds/charge.h":

	ctypedef struct t_charge:
		t_scalar_grid2d rho
		t_scalar_grid2d neutral
		t_cscalar_grid2d frho
		float box[2]
		float dx[2];
		float dt
		int iter;
		t_fftr2d_cfg *fft_forward

	void charge_report( const t_charge *charge )

#########################################################################################
# Current
#

cdef extern from "../../em2ds/current.h":

	ctypedef struct t_current:
		t_float3_grid2d J
		t_cfloat3_grid2d fJt
		float box[2]
		float dx[2]
		float dt
		int iter
		t_fftr2d_cfg *fft_forward
		t_fftr2d_cfg *fft_backward

	void current_report( const t_current *current, const char jc )

#########################################################################################
# EMF
#
cdef extern from "../../em2ds/emf.h":

	cdef enum emf_diag:
		EFLD, BFLD

	ctypedef struct t_emf:
		t_float3_grid2d E
		t_float3_grid2d B
		t_cfloat3_grid2d fEl
		t_cfloat3_grid2d fEt
		t_cfloat3_grid2d fB
		float box[2]
		float dx[2]
		float dt
		int iter
		t_fftr2d_cfg *fft_forward
		t_fftr2d_cfg *fft_backward

	cdef enum emf_laser_type:
		PLANE, GAUSSIAN

	ctypedef struct t_emf_laser:
		emf_laser_type type
		float start
		float fwhm
		float rise
		float flat
		float fall
		float a0
		float omega0
		float polarization
		float W0
		float focus
		float axis

	void emf_report( const t_emf *emf, const char field, const char fc )
	void emf_get_energy( const t_emf *emf, double energy[] )

#########################################################################################
# Simulation
#

cdef extern from "../../em2ds/simulation.h":
	ctypedef struct t_simulation:
		float dt
		float tmax
		int ndump
		int n_species
		t_species* species
		t_emf emf
		t_current current
		t_charge charge

	void sim_new( t_simulation* sim, int nx[], float box[], float dt, float tmax, int ndump, t_species* species, int n_species )

	void sim_add_laser( t_simulation* sim,  t_emf_laser* laser )
	void sim_add_neutral_bkg( t_simulation* sim )

	void sim_iter( t_simulation* sim )
	void sim_report_energy( t_simulation* sim )

	void sim_delete( t_simulation* sim )


