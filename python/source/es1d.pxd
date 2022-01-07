#cython: language_level=3

from libc.stdint cimport uint32_t

cdef extern from "../../es1d/random.h":
	void set_rand_seed( uint32_t m_z_, uint32_t m_w_ )

#########################################################################################
# Grids
#
cdef extern from "../../es1d/grid.h":
	ctypedef struct t_scalar_grid:
		float *s
		float *buffer
		int nx
		int gc[2]

	ctypedef struct t_cscalar_grid:
		float complex *s
		float complex *buffer
		int nx
		int gc[2]

#########################################################################################
# FFT
#
cdef extern from "../../es1d/fft.h":
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

#########################################################################################
# Particles
#
cdef extern from "../../es1d/particles.h":
	cdef enum:
		MAX_SPNAME_LEN

	ctypedef struct t_part:
		int ix
		float x
		float vx

	cdef enum density_type:
		UNIFORM, EMPTY, STEP, SLAB, RAMP, CUSTOM

	ctypedef struct t_density:
		float n
		density_type type
		float start
		float end
		float ramp[2]
		float (*custom)(float, void*)
		void *custom_data
		unsigned long total_np_inj
		double custom_q_inj


	ctypedef struct t_species:
		char name[MAX_SPNAME_LEN]
		t_part *part
		int np
		int np_max
		float m_q
		double energy
		float q
		int ppc
		t_density density
		float vfl
		float vth
		int nx
		float dx
		float box
		float dt
		int iter
		int n_sort

	void spec_new( t_species* spec, char name[], const float m_q, const int ppc,
				  const float* vfl, const float* vth,
				  const int nx, float box, const float dt, t_density* density )

	void spec_grow_buffer( t_species* spec, const int size )

	cdef int CHARGE
	cdef int PHA
	cdef int PARTICLES
	cdef int X1
	cdef int V1

	int PHASESPACE( int a, int b )

	void spec_report( const t_species *spec, const int rep_type,
				  const int pha_nx[], const float pha_range[][2] )

	void spec_deposit_charge( const t_species* spec, float* charge )

	void spec_deposit_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2], float* buf )


#########################################################################################
# Field
#
cdef extern from "../../es1d/field.h":

	ctypedef struct t_field:
		t_scalar_grid E
		t_cscalar_grid fE
		float box
		float dx
		float dt
		int iter
		t_fftr_cfg *fft_backward

	void field_report( const t_field *field )

#########################################################################################
# Charge
#

cdef extern from "../../es1d/charge.h":

	ctypedef struct t_charge:
		t_scalar_grid rho
		t_scalar_grid neutral
		t_cscalar_grid frho
		float box
		float dx;
		float dt
		int iter;
		t_fftr_cfg *fft_forward

	void charge_report( const t_charge *charge )


#########################################################################################
# Simulation
#

cdef extern from "../../es1d/simulation.h":
	ctypedef struct t_simulation:
		float dt
		float tmax
		int ndump
		int n_species
		t_species* species
		t_field field
		t_charge charge

		t_fftr_cfg fft_forward
		t_fftr_cfg fft_backward

	void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species )

	void sim_iter( t_simulation* sim )
	void sim_report_energy( t_simulation* sim )

	void sim_delete( t_simulation* sim )

	void sim_add_neutral_bkg( t_simulation* sim )


