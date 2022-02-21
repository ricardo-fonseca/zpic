#cython: language_level=3

from libc.stdint cimport uint32_t

cdef extern from "../../em1ds/random.h":
	void set_rand_seed( uint32_t m_z_, uint32_t m_w_ )


#########################################################################################
# Grids
#
cdef extern from "../../em1ds/grid.h":
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

	ctypedef struct float3:
		float x
		float y
		float z

	ctypedef struct t_float3_grid:
		float *x
		float *y
		float *z
		float *buffer
		int nx
		int gc[2]

	ctypedef struct t_cfloat3_grid:
		float complex *x
		float complex *y
		float complex *z
		float complex *buffer
		int nx
		int gc[2]

#########################################################################################
# FFT
#
cdef extern from "../../em1ds/fft.h":
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
# Spectral Filter
#
cdef extern from "../../em1ds/filter.h":
	cdef enum filter_type:
		FILTER_NONE, FILTER_GAUSS, FILTER_SHARP

	ctypedef struct t_filter:
		int type
		int nk
		float *Sk

#########################################################################################
# Particles
#
cdef extern from "../../em1ds/particles.h":
	cdef enum:
		MAX_SPNAME_LEN

	ctypedef struct t_part:
		int ix
		float x
		float ux
		float uy
		float uz

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
		float q
		double energy
		int ppc
		t_density density
		float ufl[3]
		float uth[3]
		int nx
		float dx
		float box
		int iter
		float dt
		int n_sort

	void spec_new( t_species* spec, char name[], const float m_q, const int ppc,
				  const float ufl[], const float uth[],
				  const int nx, float box, const float dt, t_density* density )

	void spec_grow_buffer( t_species* spec, const int size )

	cdef int CHARGE
	cdef int PHA
	cdef int PARTICLES
	cdef int X1
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

cdef extern from "../../em1ds/charge.h":

	ctypedef struct t_charge:
		t_scalar_grid rho
		t_scalar_grid neutral
		t_cscalar_grid frho
		float box
		float dx
		float dt
		int iter;
		t_fftr_cfg *fft_forward
		t_filter *filter

	void charge_report( const t_charge *charge )

#########################################################################################
# Current
#

cdef extern from "../../em1ds/current.h":

	ctypedef struct t_current:
		t_float3_grid J
		t_cfloat3_grid fJt
		float box
		float dx
		float dt
		int iter
		t_fftr_cfg *fft_forward
		t_filter *filter

	void current_report( const t_current *current, const char jc )

#########################################################################################
# EMF
#
cdef extern from "../../em1ds/emf.h":

	cdef enum emf_fld_type:
		EMF_FLD_TYPE_NONE, EMF_FLD_TYPE_UNIFORM, EMF_FLD_TYPE_CUSTOM

	ctypedef struct t_emf_ext_fld:
		emf_fld_type E_type
		emf_fld_type B_type
		float3 E_0
		float3 B_0
		float3 (*E_custom)(int, float, void*)
		float3 (*B_custom)(int, float, void*)
		void *E_custom_data
		void *B_custom_data
		t_float3_grid E_part_buf
		t_float3_grid B_part_buf

	ctypedef struct t_emf_init_fld:
		emf_fld_type E_type
		emf_fld_type B_type
		float3 E_0
		float3 B_0
		float3 (*E_custom)(int, float, void*)
		float3 (*B_custom)(int, float, void*)
		void *E_custom_data
		void *B_custom_data

	cdef enum emf_diag:
		EFLD, BFLD

	cdef enum emf_solver:
		EMF_SOLVER_PSTD, EMF_SOLVER_PSATD

	ctypedef struct t_emf:
		t_float3_grid E
		t_float3_grid B
		t_cfloat3_grid fEl
		t_cfloat3_grid fEt
		t_cfloat3_grid fB
		float box
		float dx
		float dt
		int iter
		t_fftr_cfg *fft_forward
		t_fftr_cfg *fft_backward
		t_filter *filter
		t_float3_grid *E_part
		t_float3_grid *B_part
		t_emf_ext_fld ext_fld
		int solver_type

	ctypedef struct t_emf_laser:
		float start
		float fwhm
		float rise
		float flat
		float fall
		float a0
		float omega0
		float polarization

	void emf_report( const t_emf *emf, const char field, const char fc )
	void emf_get_energy( const t_emf *emf, double energy[] )
	void emf_init_fld( t_emf* const emf, t_emf_init_fld* init_fld )
	void emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* ext_fld )

#########################################################################################
# Simulation
#

cdef extern from "../../em1ds/simulation.h":
	ctypedef struct t_simulation:
		float dt
		float tmax
		int ndump
		int n_species
		t_species* species
		t_emf emf
		t_current current
		t_charge charge
		t_fftr_cfg *fft_forward
		t_fftr_cfg *fft_backward
		t_filter *filter

	void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species )

	int  sim_filter_set( t_simulation* sim, int type, float ck )
	void sim_add_laser( t_simulation* sim,  t_emf_laser* laser )
	void sim_add_neutral_bkg( t_simulation* sim )

	void sim_iter( t_simulation* sim )
	void sim_report_energy( t_simulation* sim )

	void sim_delete( t_simulation* sim )


