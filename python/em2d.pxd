#cython: language_level=3

from libc.stdint cimport uint32_t

cdef extern from "../em2d/random.h":
	void set_rand_seed( uint32_t m_z_, uint32_t m_w_ )

cdef extern from "../em2d/zpic.h":
	ctypedef struct t_vfld:
		float x
		float y
		float z

#########################################################################################
# Particles
#
cdef extern from "../em2d/particles.h":
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
		UNIFORM, STEP, SLAB, CUSTOM

	ctypedef struct t_density:
		float n
		density_type type
		float start
		float end
		float (*custom_x)(float, void*)
		void *custom_data_x
		float (*custom_y)(float, void*)
		void *custom_data_y
		unsigned long custom_x_total_part
		double custom_x_total_q

	ctypedef struct t_species:
		char name[MAX_SPNAME_LEN]
		t_part *part
		int np
		int np_max
		float m_q
		double energy
		float q
		int ppc[2]
		t_density density
		float ufl[3]
		float uth[3]
		int nx[2]
		float dx[2]
		float box[2]
		float dt
		int iter
		int moving_window
		int n_move
		int sort_t

	void spec_new( t_species* spec, char name[], const float m_q, const int ppc[],
				  const float ufl[], const float uth[],
				  const int nx[], float box[], const float dt, t_density* density ,const int tsort)

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
# EMF
#
cdef extern from "../em2d/emf.h":

	cdef enum emf_ext_fld:
		EMF_EXT_FLD_NONE, EMF_EXT_FLD_UNIFORM, EMF_EXT_FLD_CUSTOM

	ctypedef struct t_emf_ext_fld:
		emf_ext_fld E_type
		emf_ext_fld B_type
		t_vfld E_0
		t_vfld B_0
		t_vfld (*E_custom)(int, float, int, float, void*)
		t_vfld (*B_custom)(int, float, int, float, void*)
		void *E_custom_data
		void *B_custom_data
		t_vfld *E_part_buf
		t_vfld *B_part_buf

	cdef enum emf_diag:
		EFLD, BFLD, EPART, BPART

	ctypedef struct t_emf:
		t_vfld *E
		t_vfld *B
		t_vfld *E_buf
		t_vfld *B_buf
		t_vfld *E_part
		t_vfld *B_part
		int nx[2]
		int nrow
		int gc[2][2]
		float box[2]
		float dx[2]
		float dt
		int iter
		int moving_window
		int n_move
		t_emf_ext_fld ext_fld

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
	void emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* ext_fld )

#########################################################################################
# Current
#

cdef extern from "../em2d/current.h":
	cdef enum smooth_type:
		NONE, BINOMIAL, COMPENSATED

	ctypedef struct t_smooth:
		smooth_type xtype
		smooth_type ytype
		int xlevel
		int ylevel

	ctypedef struct t_current:
		t_vfld *J
		t_vfld *J_buf
		int nx[2]
		int nrow
		int gc[2][2]
		float box[2]
		float dx[2];
		t_smooth smooth
		float dt
		int iter;
		int moving_window;

	void current_report( const t_current *current, const char jc )


#########################################################################################
# Simulation
#

cdef extern from "../em2d/simulation.h":
	ctypedef struct t_simulation:
		int moving_window
		float dt
		float tmax
		int ndump
		int n_species
		t_species* species
		t_emf emf
		t_current current

	void sim_new( t_simulation* sim, int nx[], float box[], float dt, float tmax, int ndump, t_species* species, int n_species )

	void sim_add_laser( t_simulation* sim,  t_emf_laser* laser )
	void sim_set_moving_window( t_simulation* sim )
	void sim_set_smooth( t_simulation* sim,  t_smooth* smooth )

	void sim_iter( t_simulation* sim )
	void sim_report_energy( t_simulation* sim )

	void sim_delete( t_simulation* sim )
