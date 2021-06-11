#cython: language_level=3

cimport em2ds
from libc.stdlib cimport calloc, free

import numpy as np

cdef class Density:
	"""Extension type to wrap t_density objects"""
	cdef t_density *_thisptr

	_density_types = {'uniform':UNIFORM,
					  'empty':EMPTY,
	                  'step':STEP,
	                  'slab':SLAB}

	def __cinit__( self, *, str type = 'uniform', float start = 0.0, float end = 0.0,
		           float n = 1.0):

		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_density *> calloc(1, sizeof(t_density))

		self._thisptr.type = <density_type> self._density_types[type]
		self._thisptr.n = n
		self._thisptr.start = start
		self._thisptr.end = end

	def __dealloc__(self):
		free( self._thisptr )

	def copy(self):
		new = Density()
		new.n	 = self.n
		new.type  = self.type
		new.start = self.start
		new.end   = self.end

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


cdef class Species:
	"""Extension type to wrap t_species objects"""

	cdef t_species _this
	cdef t_species* _thisptr
	cdef Density _density
	cdef str _name

	# Diagnostic types
	_diag_types  = { 'charge':CHARGE, 'pha':PHA, 'particles':PARTICLES }
	_pha_quants = { 'x1':X1, 'x2':X2, 'u1':U1, 'u2':U2, 'u3':U3 }

	def __cinit__( self, str name, const float m_q, list ppc = [1,1], *,
				  list ufl = [0.,0.,0.], list uth = [0.,0.,0.], Density density = None):

		self._thisptr = &self._this
		self._name = name
		self._this.m_q = m_q
		self._this.ppc = np.array(ppc, dtype=np.int32)
		self._this.ufl = np.array(ufl, dtype=np.float32)
		self._this.uth = np.array(uth, dtype=np.float32)

		if ( density ):
			self._density = density.copy()
		else:
			# Use default uniform density
			self._density = Density()

	cdef new( self, t_species* ptr, int nx[], float box[], float dt ):
		self._thisptr = ptr
		spec_new( self._thisptr, self._name.encode(), self._this.m_q, self._this.ppc,
			self._this.ufl, self._this.uth,
			nx, box, dt, self._density._thisptr )

	def add( self, int[:] ix, float[:] x, float[:] u):
		# insure we have enough room for new particle
		spec_grow_buffer( self._thisptr, self._thisptr.np + 1 )
		
		cdef t_part particle
		particle.ix = ix[0]
		particle.iy = ix[1]
		particle.x = x[0]
		particle.y = x[1]
		particle.ux = u[0]
		particle.uy = u[1]
		particle.uz = u[2]

	def report( self, str type, *, list quants = [], list pha_nx = [], list pha_range = [] ):
		cdef int _nx[2]
		cdef float _range[2][2]

		cdef int rep_type = self._diag_types[type]

		if ( rep_type == PHA ):
			# Phasespace diagnostics get special treatment
			_nx = np.array( pha_nx, dtype = np.int32)
			_range = np.array( pha_range, dtype = np.float32)
			rep_type = PHASESPACE( self._pha_quants[quants[0]],
				                   self._pha_quants[quants[1]])
			spec_report( self._thisptr, rep_type, _nx, _range )
		else:
			# Other diagnostic
			spec_report( self._thisptr, rep_type, NULL, NULL )

	@property
	def dx(self):
		return self._thisptr.dx

	@property
	def particles(self):
		cdef t_part[::1] buf = <t_part[:self._thisptr.np]>self._thisptr.part
		return np.asarray( buf )

	def charge(self):
		charge = np.zeros( shape = (self._thisptr.nx[1]+1,self._thisptr.nx[0]+1),
		                   dtype = np.float32 )
		cdef float [:,:] buf = charge
		spec_deposit_charge( self._thisptr, &buf[0,0] )

		# Throw away guard cells
		return charge[ 0 : self._thisptr.nx[1], 0 : self._thisptr.nx[0] ]

	def phasespace( self, list quants, pha_nx, pha_range ):

		cdef int _nx[2]
		cdef float _range[2][2]
		cdef int rep_type = PHASESPACE( self._pha_quants[quants[0]],
				                        self._pha_quants[quants[1]])

		_nx = np.array( pha_nx, dtype = np.int32)
		_range = np.array( pha_range, dtype = np.float32)

		pha = np.zeros( shape = (_nx[1],_nx[0]), dtype = np.float32 )
		cdef float [:,:] buf = pha

		spec_deposit_pha( self._thisptr, rep_type, _nx, _range, &buf[0,0] )

		return pha


def phasespace( int a, int b ):
	"""Returns the type of the requested phasespace"""
	return PHASESPACE(a,b)

cdef class EMF:
	"""Extension type to wrap t_emf objects"""

	cdef t_emf* _thisptr

	# Diagnostic types
	_diag_types = { 'E' : EFLD,	'B' : BFLD }

	cdef associate( self, t_emf* ptr ):
		self._thisptr = ptr

	def report( self, str type, char fc ):
		cdef int rep_type = self._diag_types[type];
		emf_report( self._thisptr, rep_type, fc )

	def get_energy( self ):
		cdef double energy[6]
		emf_get_energy( self._thisptr, energy )
		return np.array( energy, dtype = np.float64 )

	@property
	def nx(self):
		return self._thisptr.E.nx

	@property
	def dx(self):
		return self._thisptr.dx

	@property
	def box(self):
		return self._thisptr.box

	@property
	def Ex( self ):
		cdef int nx = self._thisptr.E.gc[0][0] + self._thisptr.E.nx[0] + self._thisptr.E.gc[0][1]
		cdef int ny = self._thisptr.E.gc[1][0] + self._thisptr.E.nx[1] + self._thisptr.E.gc[1][1]

		cdef float *buf = <float *> self._thisptr.E.x - self._thisptr.E.gc[0][0] - self._thisptr.E.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.E.gc[1][0]:self._thisptr.E.gc[1][0] + self._thisptr.E.nx[1], \
		            self._thisptr.E.gc[0][0]:self._thisptr.E.gc[0][0] + self._thisptr.E.nx[0] ]

	@property
	def Ey( self ):
		cdef int nx = self._thisptr.E.gc[0][0] + self._thisptr.E.nx[0] + self._thisptr.E.gc[0][1]
		cdef int ny = self._thisptr.E.gc[1][0] + self._thisptr.E.nx[1] + self._thisptr.E.gc[1][1]

		cdef float *buf = <float *> self._thisptr.E.y - self._thisptr.E.gc[0][0] - self._thisptr.E.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.E.gc[1][0]:self._thisptr.E.gc[1][0] + self._thisptr.E.nx[1], \
		            self._thisptr.E.gc[0][0]:self._thisptr.E.gc[0][0] + self._thisptr.E.nx[0] ]

	@property
	def Ez( self ):
		cdef int nx = self._thisptr.E.gc[0][0] + self._thisptr.E.nx[0] + self._thisptr.E.gc[0][1]
		cdef int ny = self._thisptr.E.gc[1][0] + self._thisptr.E.nx[1] + self._thisptr.E.gc[1][1]

		cdef float *buf = <float *> self._thisptr.E.z - self._thisptr.E.gc[0][0] - self._thisptr.E.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.E.gc[1][0]:self._thisptr.E.gc[1][0] + self._thisptr.E.nx[1], \
		            self._thisptr.E.gc[0][0]:self._thisptr.E.gc[0][0] + self._thisptr.E.nx[0] ]


	@property
	def Bx( self ):
		cdef int nx = self._thisptr.B.gc[0][0] + self._thisptr.B.nx[0] + self._thisptr.B.gc[0][1]
		cdef int ny = self._thisptr.B.gc[1][0] + self._thisptr.B.nx[1] + self._thisptr.B.gc[1][1]

		cdef float *buf = <float *> self._thisptr.B.x - self._thisptr.B.gc[0][0] - self._thisptr.B.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.B.gc[1][0]:self._thisptr.B.gc[1][0] + self._thisptr.B.nx[1], \
		            self._thisptr.B.gc[0][0]:self._thisptr.B.gc[0][0] + self._thisptr.B.nx[0] ]

	@property
	def By( self ):
		cdef int nx = self._thisptr.B.gc[0][0] + self._thisptr.B.nx[0] + self._thisptr.B.gc[0][1]
		cdef int ny = self._thisptr.B.gc[1][0] + self._thisptr.B.nx[1] + self._thisptr.B.gc[1][1]

		cdef float *buf = <float *> self._thisptr.B.y - self._thisptr.B.gc[0][0] - self._thisptr.B.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.B.gc[1][0]:self._thisptr.B.gc[1][0] + self._thisptr.B.nx[1], \
		            self._thisptr.B.gc[0][0]:self._thisptr.B.gc[0][0] + self._thisptr.B.nx[0] ]

	@property
	def Bz( self ):
		cdef int nx = self._thisptr.B.gc[0][0] + self._thisptr.B.nx[0] + self._thisptr.B.gc[0][1]
		cdef int ny = self._thisptr.B.gc[1][0] + self._thisptr.B.nx[1] + self._thisptr.B.gc[1][1]

		cdef float *buf = <float *> self._thisptr.B.z - self._thisptr.B.gc[0][0] - self._thisptr.B.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.B.gc[1][0]:self._thisptr.B.gc[1][0] + self._thisptr.B.nx[1], \
		            self._thisptr.B.gc[0][0]:self._thisptr.B.gc[0][0] + self._thisptr.B.nx[0] ]


cdef class Laser:
	"""Extension type to wrap t_emf_laser objects"""

	cdef t_emf_laser * _thisptr

	# Laser types
	_laser_types = {'plane':PLANE,'gaussian':GAUSSIAN }

	def __cinit__( self, *, type = 'plane', float start = 0.0, float fwhm = 0.0,
		           float rise = 0.0, float flat = 0.0, float fall = 0.0,
	               float a0 = 0.0, float omega0 = 0.0, float polarization = 0.0,
	               float W0 = 0.0, float focus = 0.0, float axis = 0.0 ):
		self._thisptr = <t_emf_laser *> calloc(1, sizeof(t_emf_laser))

		self._thisptr.type = self._laser_types[type]

		self._thisptr.start = start
		self._thisptr.fwhm = fwhm
		self._thisptr.rise = rise
		self._thisptr.flat = flat
		self._thisptr.fall = fall
		self._thisptr.a0 = a0
		self._thisptr.omega0 = omega0
		self._thisptr.polarization = polarization
		self._thisptr.W0 = W0
		self._thisptr.focus = focus
		self._thisptr.axis = axis

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

	@property
	def W0(self):
		return self._thisptr.W0

	@W0.setter
	def W0(self,value):
		self._thisptr.W0 = value

	@property
	def focus(self):
		return self._thisptr.focus

	@focus.setter
	def focus(self,value):
		self._thisptr.focus = value

	@property
	def axis(self):
		return self._thisptr.axis

	@axis.setter
	def axis(self,value):
		self._thisptr.axis = value


cdef class Current:
	"""Extension type to wrap t_current objects"""

	cdef t_current* _thisptr

	cdef associate( self, t_current* ptr ):
		self._thisptr = ptr

	def report( self, char jc ):
		current_report( self._thisptr, jc )

	@property
	def Jx( self ):
		cdef int nx = self._thisptr.J.gc[0][0] + self._thisptr.J.nx[0] + self._thisptr.J.gc[0][1]
		cdef int ny = self._thisptr.J.gc[1][0] + self._thisptr.J.nx[1] + self._thisptr.J.gc[1][1]

		cdef float *buf = <float *> self._thisptr.J.x - self._thisptr.J.gc[0][0] - self._thisptr.J.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.J.gc[1][0]:self._thisptr.J.gc[1][0] + self._thisptr.J.nx[1], \
		            self._thisptr.J.gc[0][0]:self._thisptr.J.gc[0][0] + self._thisptr.J.nx[0] ]

	@property
	def Jy( self ):
		cdef int nx = self._thisptr.J.gc[0][0] + self._thisptr.J.nx[0] + self._thisptr.J.gc[0][1]
		cdef int ny = self._thisptr.J.gc[1][0] + self._thisptr.J.nx[1] + self._thisptr.J.gc[1][1]

		cdef float *buf = <float *> self._thisptr.J.y - self._thisptr.J.gc[0][0] - self._thisptr.J.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.J.gc[1][0]:self._thisptr.J.gc[1][0] + self._thisptr.J.nx[1], \
		            self._thisptr.J.gc[0][0]:self._thisptr.J.gc[0][0] + self._thisptr.J.nx[0] ]

	@property
	def Jz( self ):
		cdef int nx = self._thisptr.J.gc[0][0] + self._thisptr.J.nx[0] + self._thisptr.J.gc[0][1]
		cdef int ny = self._thisptr.J.gc[1][0] + self._thisptr.J.nx[1] + self._thisptr.J.gc[1][1]

		cdef float *buf = <float *> self._thisptr.J.z - self._thisptr.J.gc[0][0] - self._thisptr.J.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.J.gc[1][0]:self._thisptr.J.gc[1][0] + self._thisptr.J.nx[1], \
		            self._thisptr.J.gc[0][0]:self._thisptr.J.gc[0][0] + self._thisptr.J.nx[0] ]

cdef class Charge:
	"""Extension type to wrap t_charge objects"""

	cdef t_charge* _thisptr

	cdef associate( self, t_charge* ptr ):
		self._thisptr = ptr

	def report( self ):
		charge_report( self._thisptr )

	@property
	def rho( self ):
		cdef int nx = self._thisptr.rho.gc[0][0] + self._thisptr.rho.nx[0] + self._thisptr.rho.gc[0][1]
		cdef int ny = self._thisptr.rho.gc[1][0] + self._thisptr.rho.nx[1] + self._thisptr.rho.gc[1][1]

		cdef float *buf = <float *> self._thisptr.rho.s - self._thisptr.rho.gc[0][0] - self._thisptr.rho.gc[1][0] * nx

		tmp = np.asarray( <float [:ny,:nx]> buf )
		return tmp[ self._thisptr.rho.gc[1][0]:self._thisptr.rho.gc[1][0] + self._thisptr.rho.nx[1], \
		            self._thisptr.rho.gc[0][0]:self._thisptr.rho.gc[0][0] + self._thisptr.rho.nx[0] ]


cdef class Simulation:
	"""Extension type to wrap t_simulation objects"""

	cdef t_simulation *_thisptr
	cdef int n
	cdef float t

	cdef EMF emf
	cdef Current current
	cdef Charge charge
	cdef list species

	cdef object report

	def __cinit__( self, list nx, list box, float dt, *, species = None,
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

		cdef int[2] _nx = np.array( nx, dtype = np.int32)
		cdef float[2] _box = np.array( box, dtype = np.float32)

		if ( isinstance( species, Species )):
			n_species = 1
			species_ = <t_species *> calloc( 1, sizeof(t_species))
			s = species
			s.new( &species_[0], _nx, _box, dt )
			self.species.append( s )

		elif ( isinstance( species, (list,tuple) ) ):
			n_species = len( species )
			species_ = <t_species *> calloc( n_species, sizeof(t_species))
			for i in range(n_species):
				s = species[i]
				s.new( &species_[i], _nx, _box, dt )
				self.species.append( s )
		else:
			n_species = 0
			species_ = NULL

		# Diagnostics
		self.report = report

		# Initialize simulation
		sim_new( self._thisptr, _nx, _box, dt, 0.0, 0, species_, n_species )

		self.n = 0
		self.t = 0.0

		self.emf = EMF()
		self.emf.associate( &self._thisptr.emf )

		self.current = Current()
		self.current.associate( &self._thisptr.current )

		self.charge = Charge()
		self.charge.associate( &self._thisptr.charge )

	def __dealloc__(self):
		sim_delete( self._thisptr )
		free(self._thisptr)

	def add_laser(self, Laser laser):
		sim_add_laser( self._thisptr, laser._thisptr )

	def add_neutral_bkg(self):
		sim_add_neutral_bkg( self._thisptr )

	def iter( self ):
		sim_iter( self._thisptr )
		self.n = self.n+1
		self.t = self.n * self._thisptr.dt

	def run( self, float tmax ):

		if ( tmax < self.t ):
			print("Simulation is already at t = {:g}".format(self.t))
			return

		print("\nRunning simulation up to t = {:g} ...".format(tmax))

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
	def charge(self):
		return self.charge

	@property
	def n(self):
		return self.n

	@property
	def t(self):
		return self.t

	@property
	def dx(self):
		return self.emf.dx

	@property
	def nx(self):
		return self.emf.nx

	@property
	def box(self):
		return self.emf.box

	@property
	def report(self):
		return self.report

	@report.setter
	def report( self, f ):
		self.report = f





