#cython: language_level=3

from libc.stdint cimport uint32_t

cdef extern from "../../em1d/random.h":
	void set_rand_seed( uint32_t m_z_, uint32_t m_w_ )

cdef extern from "../../em1d/zpic.h":
	ctypedef struct float3:
		float x
		float y
		float z

#########################################################################################
# Particles
#
cdef extern from "../../em1d/particles.h":
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

	cdef enum part_boundary:
		PART_BC_NONE, PART_BC_PERIODIC, PART_BC_OPEN

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
		float ufl[3]
		float uth[3]
		int nx
		float dx
		float box
		float dt
		int iter
		int moving_window
		int n_move
		int bc_type
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
# EMF
#
cdef extern from "../../em1d/emf.h":

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
		float3 *E_part_buf
		float3 *B_part_buf

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
		EFLD, BFLD, EPART, BPART

	cdef enum emf_boundary:
		EMF_BC_NONE, EMF_BC_PERIODIC, EMF_BC_OPEN

	ctypedef struct t_emf:
		float3 *E
		float3 *B
		float3 *E_buf
		float3 *B_buf
		float3 *E_part
		float3 *B_part
		int nx
		int gc[2]
		float box
		float dx
		float dt
		int iter
		int moving_window
		int n_move
		int bc_type
		float3 mur_fld[2]
		float3 mur_tmp[2]
		t_emf_ext_fld ext_fld

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
# Current
#

cdef extern from "../../em1d/current.h":
	cdef enum smooth_type:
		NONE, BINOMIAL, COMPENSATED

	cdef enum current_boundary:
		CURRENT_BC_NONE, CURRENT_BC_PERIODIC

	ctypedef struct t_smooth:
		smooth_type xtype
		int xlevel

	ctypedef struct t_current:
		float3 *J
		float3 *J_buf
		int nx
		int gc[2]
		float box
		float dx;
		t_smooth smooth
		float dt
		int iter;
		int bc_type;

	void current_report( const t_current *current, const char jc )


#########################################################################################
# Simulation
#

cdef extern from "../../em1d/simulation.h":
	ctypedef struct t_simulation:
		int moving_window
		float dt
		float tmax
		int ndump
		int n_species
		t_species* species
		t_emf emf
		t_current current

	void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species )

	void sim_add_laser( t_simulation* sim,  t_emf_laser* laser )
	void sim_set_moving_window( t_simulation* sim )
	void sim_set_smooth( t_simulation* sim,  t_smooth* smooth )

	void sim_iter( t_simulation* sim )
	void sim_report_energy( t_simulation* sim )

	void sim_delete( t_simulation* sim )


