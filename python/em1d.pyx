cimport em1d
from libc.stdlib cimport calloc, free

import numpy as np

cdef class Density:
	cdef t_density *_thisptr

	def __cinit__( self ):
		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_density *> calloc(1, sizeof(t_density))

	def __dealloc__(self):
		free( self._thisptr )

	def copy(self):
		new = Density()
		new.n	 = self.n
		new.type  = self.type
		new.start = self.start
		new.end   = self.end
		new.ramp  = self.ramp
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
	cdef t_species _this
	cdef t_species* _thisptr
	cdef Density _density
	cdef str _name

	charge	= CHARGE
	pha	   = PHA
	particles = PARTICLES
	x1		= X1
	u1		= U1
	u2		= U2
	u3		= U3

	def __cinit__( self, str name, const float m_q, const int ppc,
				  list ufl, list uth, Density density = None):

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
		if ( type == PARTICLES ):
			spec_report( self._thisptr, PARTICLES, NULL, NULL )
		else:
			print("Unsupported species diagnostic type")

cdef class EMF:
	cdef t_emf* _thisptr

	efld = EFLD
	bfld = BFLD

	cdef associate( self, t_emf* ptr ):
		self._thisptr = ptr

	def report( self, char field, char fc ):
		emf_report( self._thisptr, field, fc )

cdef class Simulation:
	cdef t_simulation *_thisptr
	cdef int n
	cdef float t

	def __cinit__( self, int nx, float box, float dt, float tmax, list species ):
		# Allocate the simulation object
		self._thisptr = <t_simulation *> calloc(1, sizeof(t_simulation))

		# Allocate the species array
		cdef int n_species = len(species)
		cdef t_species* species_ = <t_species *> calloc( n_species, sizeof(t_species))

		# Initialize each species
		cdef Species s
		for i in range(n_species):
			s = species[i]
			s.new( &species_[i], nx, box, dt )

		# Initialize simulation
		sim_new( self._thisptr, nx, box, dt, tmax, 0, species_, n_species )

		self.n = 0
		self.t = 0.0

	def __dealloc__(self):
		sim_delete( self._thisptr )
		free(self._thisptr)

	def run( self, report ):

		print("Starting simulation...")

		while self.t <= self._thisptr.tmax:
			print("n = {:d}, t = {:g}".format(self.n,self.t))
			report( self )
			sim_iter( self._thisptr )
			self.n = self.n+1
			self.t = self.n * self._thisptr.dt

		print("\nSimulation completed.")

	@property
	def emf(self):
		obj = EMF()
		obj.associate( &self._thisptr.emf )
		return obj







