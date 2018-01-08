#cython: language_level=3

cimport em1d
from libc.stdlib cimport calloc, free

import numpy as np

cdef float custom_density( float x, void *f ):
	cdef Density d = <object> f
	return d.custom_func(x)

cdef class Density:
	"""Extension type to wrap t_density objects"""
	cdef t_density *_thisptr

	cdef object custom_func

	uniform = UNIFORM
	step = STEP
	slab = SLAB
	# Because of a namespace conflict we cannot call this "ramp"
	linear_ramp = RAMP
	custom = CUSTOM

	def __cinit__( self, *, type = UNIFORM, float start = 0.0, float end = 0.0,
		           list ramp = [0.,0.], custom = None):
		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_density *> calloc(1, sizeof(t_density))

		self._thisptr.type = type
		self._thisptr.start = start
		self._thisptr.end = end
		self._thisptr.ramp = np.array(ramp, dtype=np.float32)
		if ( custom ):
			self.custom_func = custom
			self._thisptr.custom = custom_density
			self._thisptr.custom_data = <void *> self

	def __dealloc__(self):
		free( self._thisptr )

	def copy(self):
		new = Density()
		new.n	 = self.n
		new.type  = self.type
		new.start = self.start
		new.end   = self.end
		new.ramp  = self.ramp

		new.custom_func = self.custom_func
		new._thisptr.custom = self._thisptr.custom
		new._thisptr.custom_data = self._thisptr.custom_data

		return new

	@property
	def n(self):
		return self._thisptr.n

	@n.setter
	def n(self,value):
		self._thisptr.n = value

	@property
	def type(self):
		return self._thisptr.type

	@type.setter
	def type(self,value):
		self._thisptr.type = value

	@property
	def start(self):
		return self._thisptr.start

	@start.setter
	def start(self,value):
		self._thisptr.start = value

	@property
	def end(self):
		return self._thisptr.end

	@end.setter
	def end(self,value):
		self._thisptr.end = value

	@property
	def ramp(self):
		return self._thisptr.ramp

	@ramp.setter
	def ramp(self,value):
		self._thisptr.ramp = value

cdef class Species:
	"""Extension type to wrap t_species objects"""

	cdef t_species _this
	cdef t_species* _thisptr
	cdef Density _density
	cdef str _name

	# Diagnostic types
	charge	= CHARGE
	pha	   = PHA
	particles = PARTICLES
	x1		= X1
	u1		= U1
	u2		= U2
	u3		= U3

	# Boundary condition types
	bc_none     = PART_BC_NONE
	bc_periodic = PART_BC_PERIODIC
	bc_open     = PART_BC_OPEN

	def __cinit__( self, str name, const float m_q, const int ppc, *,
				  list ufl = [0.,0.,0.], list uth = [0.,0.,0.], Density density = None):

		self._thisptr = &self._this
		self._name = name
		self._this.m_q = m_q
		self._this.ppc = ppc
		self._this.ufl = np.array(ufl, dtype=np.float32)
		self._this.uth = np.array(uth, dtype=np.float32)

		if ( density ):
			self._density = density.copy()
		else:
			self._density._thisptr = NULL

	cdef new( self, t_species* ptr, int nx, float box, float dt ):
		self._thisptr = ptr
		spec_new( self._thisptr, self._name.encode(), self._this.m_q, self._this.ppc,
			self._this.ufl, self._this.uth,
			nx, box, dt, self._density._thisptr )

	def report( self, type, *, pha_nx = 0, pha_range = 0 ):
		cdef int _nx[2]
		cdef float _range[2][2]

		if ( type == PARTICLES or type == CHARGE ):
			spec_report( self._thisptr, type, NULL, NULL )
		else:
			# Phasespace diagnostic
			_nx = np.array( pha_nx, dtype = np.int32)
			_range = np.array( pha_range, dtype = np.float32)
			spec_report( self._thisptr, type, _nx, _range )


def phasespace( int a, int b ):
	"""Returns the type of the requested phasespace"""
	return PHASESPACE(a,b)


cdef class EMF:
	"""Extension type to wrap t_emf objects"""

	cdef t_emf* _thisptr

	# Diagnostic types
	efld = EFLD
	bfld = BFLD

	# Boundary condition types
	bc_none     = EMF_BC_NONE
	bc_periodic = EMF_BC_PERIODIC
	bc_open     = EMF_BC_OPEN

	cdef associate( self, t_emf* ptr ):
		self._thisptr = ptr

	def report( self, char field, char fc ):
		emf_report( self._thisptr, field, fc )

	@property
	def bc_type( self ):
		return self._thisptr.bc_type

	@bc_type.setter
	def bc_type( self, value ):
		self._thisptr.bc_type = value

	@property
	def E( self ):
		cdef float *buf = <float *> self._thisptr.E
		return np.asarray( <float [:self._thisptr.nx,:3]> buf, dtype = np.float32 )

	@property
	def B( self ):
		cdef float *buf = <float *> self._thisptr.B
		return np.asarray( <float [:self._thisptr.nx,:3]> buf, dtype = np.float32 )

cdef class Laser:
	"""Extension type to wrap t_emf_laser objects"""

	cdef t_emf_laser * _thisptr

	def __cinit__( self, *, float start = 0.0, float fwhm = 0.0,
		           float rise = 0.0, float flat = 0.0, float fall = 0.0,
	               float a0 = 0.0, float omega0 = 0.0, float polarization = 0.0 ):
		self._thisptr = <t_emf_laser *> calloc(1, sizeof(t_emf_laser))

		self._thisptr.start = start
		self._thisptr.fwhm = fwhm
		self._thisptr.rise = rise
		self._thisptr.flat = flat
		self._thisptr.fall = fall
		self._thisptr.a0 = a0
		self._thisptr.omega0 = omega0
		self._thisptr.polarization = polarization

	def __dealloc__(self):
		free( self._thisptr )

	@property
	def start(self):
		return self._thisptr.start

	@start.setter
	def start(self,value):
		self._thisptr.start = value

	@property
	def fwhm(self):
		return self._thisptr.fwhm

	@fwhm.setter
	def fwhm(self,value):
		self._thisptr.fwhm = value

	@property
	def rise(self):
		return self._thisptr.rise

	@rise.setter
	def rise(self,value):
		self._thisptr.rise = value

	@property
	def flat(self):
		return self._thisptr.flat

	@flat.setter
	def flat(self,value):
		self._thisptr.flat = value

	@property
	def fall(self):
		return self._thisptr.fall

	@fall.setter
	def fall(self,value):
		self._thisptr.fall = value

	@property
	def a0(self):
		return self._thisptr.a0

	@a0.setter
	def a0(self,value):
		self._thisptr.a0 = value

	@property
	def omega0(self):
		return self._thisptr.omega0

	@omega0.setter
	def omega0(self,value):
		self._thisptr.omega0 = value

	@property
	def polarization(self):
		return self._thisptr.polarization

	@polarization.setter
	def polarization(self,value):
		self._thisptr.polarization = value


cdef class Current:
	"""Extension type to wrap t_current objects"""

	cdef t_current* _thisptr

	cdef associate( self, t_current* ptr ):
		self._thisptr = ptr

	def report( self, char jc ):
		current_report( self._thisptr, jc )

	@property
	def J( self ):
		cdef float *buf = <float *> self._thisptr.J
		return np.asarray( <float [:self._thisptr.nx,:3]> buf, dtype = np.float32 )


cdef class Smooth:
	"""Extension type to wrap t_smooth objects"""

	cdef t_smooth* _thisptr

	none        = NONE
	binomial    = BINOMIAL
	compensated = COMPENSATED

	def __cinit__( self, *, int xtype = NONE, int xlevel = 0 ):
		self._thisptr = <t_smooth *> calloc(1, sizeof(t_smooth))

		self._thisptr.xtype = <smooth_type>  xtype
		self._thisptr.xlevel = xlevel

	def __dealloc__(self):
		free( self._thisptr )

	@property
	def xtype(self):
		return self._thisptr.xtype

	@xtype.setter
	def xtype( self, value ):
		self._thisptr.xtype = value

	@property
	def xlevel(self):
		return self._thisptr.xlevel

	@xlevel.setter
	def xlevel( self, value ):
		self._thisptr.xlevel = value


cdef class Simulation:
	"""Extension type to wrap t_simulation objects"""

	cdef t_simulation *_thisptr
	cdef int n
	cdef float t

	cdef EMF emf
	cdef Current current
	cdef list species

	cdef object report

	def __cinit__( self, int nx, float box, float dt, species, *,
	               report = None ):

		# Allocate the simulation object
		self._thisptr = <t_simulation *> calloc(1, sizeof(t_simulation))

		# Initialize the random number generator
		# These are the value set when launching a new C simulation
		set_rand_seed( 12345, 67890 )

		# Initialize particle species data
		self.species = []
		cdef Species s

		cdef int n_species
		cdef t_species* species_

		if ( isinstance( species, Species )):
			n_species = 1
			species_ = <t_species *> calloc( 1, sizeof(t_species))
			s = species
			s.new( &species_[0], nx, box, dt )
			self.species.append( s )

		elif ( isinstance( species, (list,tuple) ) ):
			n_species = len( species )
			species_ = <t_species *> calloc( n_species, sizeof(t_species))
			for i in range(n_species):
				s = species[i]
				s.new( &species_[i], nx, box, dt )
				self.species.append( s )
		else:
			n_species = 0
			species_ = NULL

		# Diagnostics
		self.report = report

		# Initialize simulation
		sim_new( self._thisptr, nx, box, dt, 0.0, 0, species_, n_species )

		self.n = 0
		self.t = 0.0

		self.emf = EMF()
		self.emf.associate( &self._thisptr.emf )

		self.current = Current()
		self.current.associate( &self._thisptr.current )

	def __dealloc__(self):
		sim_delete( self._thisptr )
		free(self._thisptr)

	def set_moving_window(self):
		sim_set_moving_window( self._thisptr )

	def set_smooth(self, Smooth smooth):
		sim_set_smooth( self._thisptr, smooth._thisptr )

	def add_laser(self, Laser laser):
		sim_add_laser( self._thisptr, laser._thisptr )

	def run( self, float tmax ):

		if ( tmax < self.t ):
			print("Simulation is already at t = {:g}".format(self.t))
			return

		print("Running simulation up to t = {:g} ...".format(tmax))

		if ( self.report ):

			# Run simulation with diagnostics
			while self.t <= tmax:
				print('n = {:d}, t = {:g}'.format(self.n,self.t), end = '\r')
				self.report( self )
				sim_iter( self._thisptr )
				self.n = self.n+1
				self.t = self.n * self._thisptr.dt
		else:
			# Run simulation without diagnostics
			while self.t <= tmax:
				print('n = {:d}, t = {:g}'.format(self.n,self.t), end = '\r')
				sim_iter( self._thisptr )
				self.n = self.n+1
				self.t = self.n * self._thisptr.dt

		print('n = {:d}, t = {:g}'.format(self.n,self.t), end = '\r')
		print("\nDone.")

	@property
	def emf(self):
		return self.emf

	@property
	def current(self):
		return self.current

	@property
	def n(self):
		return self.n

	@property
	def t(self):
		return self.t

	@property
	def report(self):
		return self.report

	@report.setter
	def report( self, f ):
		self.report = f





