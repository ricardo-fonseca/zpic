"""# 1D electrostatic PIC code

1D, electrostatic, Particle-in-Cell code, using a spectral EM field solver
"""
#cython: language_level=3

cimport es1d
from libc.stdlib cimport calloc, free

import numpy as np
import sys

cdef float custom_density( float x, void *f ):
	"""Internal function for custom density profiles"""
	cdef Density d = <object> f
	return d.custom_func(x)

cdef class Density:
	"""Density(type='uniform', start=0.0, end=0.0, n=1.0, ramp=[0.,0.], custom=None)
	
	Class representing charge density profiles for particle species
	initialization

	Parameters
	----------
	type : str, optional
		Density profile type, one of "uniform", "empty", "step", "slab"
		or "custom", defaults to "uniform"
	start : float, optional
		Position of the plasma start position for "step", "slab" or "ramp"
		profiles, defaults to 0.0
	end : float, optional
		Position of the plasma end position for the "slab" or "ramp"
		profiles, defaults to 0.0
	n : float, optional
		Reference density to use, multiplies density profile value,
		defaults to 1.0
	ramp : list of float
		2 element list specifying the required density at the start and
		end positions for the "ramp" density profile, defaults to [0.,0.]
	custom : function, optional
		Custom density function, defaults to None
	
	See Also
	--------
	es1d.Species
	"""
	cdef t_density *_thisptr

	cdef object custom_func

	_density_types = {'uniform':UNIFORM,
					  'empty':EMPTY,
					  'step':STEP,
					  'slab':SLAB,
					  'ramp':RAMP,
					  'custom':CUSTOM}

	def __cinit__( self, *, str type = 'uniform', float n = 1.0, float start = 0.0, float end = 0.0,
		           list ramp = [0.,0.], custom = None):
		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_density *> calloc(1, sizeof(t_density))

		self._thisptr.type = <density_type> self._density_types[type]
		self._thisptr.n = n
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
		"""copy()

		Object copy.
		"""
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
		"""Reference density to use
		
		Returns
		-------
		n : float
			Reference density to use
		"""
		return self._thisptr.n

	@n.setter
	def n(self,value):
		self._thisptr.n = value

	@property
	def type(self):
		"""Density profile type

		Returns
		-------
		type : {'uniform', 'empty', 'step', 'slab', 'ramp','custom'}
			Density profile type
		"""
		tmp = {UNIFORM:'uniform', EMPTY:'empty', STEP:'step', SLAB:'slab',
			RAMP:'ramp', CUSTOM:'custom'}
		return tmp[self._thisptr.type]

	@type.setter
	def type(self,value):
		self._thisptr.type = self._density_types[value]

	@property
	def start(self):
		"""
		Position of the plasma start position for "step", "slab" or "ramp"
		profiles
		
		Returns
		-------
		start : float
			Position of the plasma start position
		"""
		return self._thisptr.start

	@start.setter
	def start(self,value):
		self._thisptr.start = value

	@property
	def end(self):
		"""
		Position of the plasma end position for "slab" or "ramp" profiles
		
		Returns
		-------
		end : float
			Position of the plasma end position
		"""
		return self._thisptr.end

	@end.setter
	def end(self,value):
		self._thisptr.end = value

	@property
	def ramp(self):
		"""
		Initial and final density values for the "ramp" profile
		
		Returns
		-------
		ramp : list of float (2)
			[start,end] density values
		"""
		return self._thisptr.ramp

	@ramp.setter
	def ramp(self,value):
		self._thisptr.ramp = value

cdef class Species:
	"""Species(name, m_q, ppc, vfl=0.0, vth=0.0, density=None, n_sort=16)
	
	Class representing particle species. Particle data can be accessed
	(read/write) using the `particles` property.
	
	Parameters
	----------
	name : str
		Name used to identify the particle species
	m_q : float
		Mass over charge ration in for particles in the species in
		simulaition units (e.g. for electrons use -1)
	ppc : list
		Number of particles per cell in the form [nx,ny]
	vfl : float, optional
		Initial fluid (generalized) velocity for the particles, 
		defaults to 0 (no fluid velocity)
	vth : float, optional
		Initial thermal velocity for the particles, defaults to 0
		(no thermal velocity).
	density : `Density`, optional
		Density profile for the particle species specified as using a
		`Density` object. Defaults to `None` which corresponds to a
		uniform density of value 1
	n_sort : int, optional
		Number of iterations between particle sort, defaults to 16
	
	See also
	--------
	es1d.Density
	"""

	cdef t_species _this
	cdef t_species* _thisptr
	cdef Density _density
	cdef str _name

	# Diagnostic types
	_diag_types  = { 'charge':CHARGE, 'pha':PHA, 'particles':PARTICLES }
	_pha_quants = { 'x1':X1, 'v1':V1 }

	def __cinit__( self, str name, const float m_q, const int ppc, *,
				  const float vfl = 0, const float vth = 0, Density density = None,
				  int n_sort = 16 ):

		self._thisptr = &self._this
		self._name = name
		self._this.m_q = m_q
		self._this.ppc = ppc
		self._this.vfl = vfl
		self._this.vth = vth

		if ( density ):
			self._density = density.copy()
		else:
			# Use default uniform density
			self._density = Density()

	def add( self, int ix, float x, float vx):
		"""add(ix, x, vx)

		Adds a new particle to the particle buffer

		Parameters
		----------
		ix : int
			New particle cell index
		x : float
			New particle position inside the cell
		vx : float
			New particle velocity
		"""
		# insure we have enough room for new particle
		spec_grow_buffer( self._thisptr, self._thisptr.np + 1 )
		
		cdef t_part particle
		particle.ix = ix
		particle.x = x
		particle.vx = vx

		self._thisptr.part[ self._thisptr.np ] = particle
		self._thisptr.np = self._thisptr.np + 1


	cdef new( self, t_species* ptr, int nx, float box, float dt ):
		self._thisptr = ptr
		spec_new( self._thisptr, self._name.encode(), self._this.m_q, self._this.ppc,
			&self._this.vfl, &self._this.vth,
			nx, box, dt, self._density._thisptr )

	def report( self, str type, *, list quants = [], list pha_nx = [], list pha_range = [] ):
		"""report(type, quants=[], pha_nx=[], pha_range=[])
		
		Saves diagnostic information to disk

		Parameters
		----------
		type : {"charge", "pha", "particles"}
			Type of information to save, must be one of "charge" (charge
			density), "pha" (phasespace density), or "particles" (raw
			particle data)
		quants : list, optional
			2 element list of quantities to use for "pha" diagnostics.
			Each quantity must be one of 'x1' or 'v1'
		pha_nx : list of int, optional
			2 element list specifying the size of the phasespace grid
		pha_range : list of float, optional
			2x2 element list specifying the physical limits of the
			phasespace grid in the form [[xmin,xmax],[ymin,ymax]]
		"""
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
	def particles(self):
		"""ndarray of particle data
		
		Allows full read/write access to particle data
		"""
		cdef t_part[::1] buf = <t_part[:self._thisptr.np]>self._thisptr.part
		return np.asarray( buf )

	def charge(self):
		"""charge()
		
		Calculate charge density of particle species

		Returns
		-------
		n : numpy.array, (nx)
			Charge density of particle species. Array will jhave the same
			shape as the simulation grid
		"""
		charge = np.zeros( shape = self._thisptr.nx+1, dtype = np.float32 )
		cdef float [::1] buf = charge
		spec_deposit_charge( self._thisptr, &buf[0] )

		# Throw away guard cell
		return charge[ 0 : self._thisptr.nx ]

	def phasespace( self, int type, pha_nx, pha_range ):
		"""phasespace(quants, pha_nx, pha_range)
		
		Calculate phasespace density of particle species

		Parameters
		----------
		quants : list of str
			2 element list of quantities to use for "pha" diagnostics.
			Each quantity must be one of 'x1' or 'v1'
		pha_nx : list of int, optional
			2 element list specifying the size of the phasespace grid
		pha_range : list of float, optional
			2x2 element list specifying the physical limits of the
			phasespace grid in the form [[xmin,xmax],[ymin,ymax]]

		Returns
		-------
		pha : numpy.ndarray, (pha_nx[0],pha_nx[1])
			Phasespace density of particle species.
		"""
		cdef int _nx[2]
		cdef float _range[2][2]

		_nx = np.array( pha_nx, dtype = np.int32)
		_range = np.array( pha_range, dtype = np.float32)

		pha = np.zeros( shape = (_nx[1],_nx[0]), dtype = np.float32 )
		cdef float [:,:] buf = pha

		spec_deposit_pha( self._thisptr, type, _nx, _range, &buf[0,0] )

		return pha

	@property
	def energy(self):
		"""Total kinetic energy of particle species

		Returns
		-------
		ene : float
			Time-centered total kinetic energy of particle species

		Note
		----
		To ensure the correct time-centering the particle kinetic energy
		is calculated during the particle advance, so it will be 0 before
		the first iteration is completed
		"""
		return self._thisptr.energy

	@property
	def dx(self):
		"""Cell size used for the species

		Returns
		-------
		dx :float
			Cell size
		"""
		return self._thisptr.dx

	@property
	def dt(self):
		"""Time step used for advancing the species

		Returns
		-------
		dt : float
			Time step in simulation units
		"""
		return self._thisptr.dt

	@property
	def iter(self):
		"""Last iteration completed by the particle species
		
		Returns
		-------
		iter : int
			Iteration number
		"""
		return self._thisptr.iter

	@property
	def ppc(self):
		"""Number of particles per cell used for initializing new particles
		from density profile.
		
		Returns
		-------
		ppc : float
			Number of particles per cell
		"""
		return self._thisptr.ppc

	@property
	def n_sort(self):
		"""Number of iterations between sorting particle data buffer

		Returns
		-------
		n_sort : int
			Number of iteration between sorts
		"""
		return self._thisptr.n_sort

	@n_sort.setter
	def n_sort(self, int value):
		if ( value < 0 ):
			print("(*error*) Invalid value for n_sort, must be >= 0.", file = sys.stderr)
			return
		self._thisptr.n_sort = value

cdef class Field:
	"""Field()
	
	Electric fields class

	This class allows access to the electric field data structures in the
	simulation. An object of this class is created automatically when
	creating a `es1d.Simulation` object.

	See Also
	--------
	es1d.Simulation
	"""
	cdef t_field* _thisptr

	cdef associate( self, t_field* ptr ):
		self._thisptr = ptr

	def report( self ):
		"""report( )

		Save diagnostic field information to disk. Files will be saved in the
		"field" directory below the current working directory.
		"""
		field_report( self._thisptr )

	@property
	def nx(self):
		"""Grid size used for the Field object

		Returns
		-------
		nx :int
			Number of grid cells for the simulation
		"""
		return self._thisptr.E.nx

	@property
	def dx(self):
		"""Cell size used for the Field object

		Returns
		-------
		dx : int
			Cell size
		"""
		return self._thisptr.dx

	@property
	def box(self):
		"""Simulation box physical size

		Returns
		-------
		box : float
			Simulation box size
		"""
		return self._thisptr.box

	@property
	def E( self ):
		"""Electric field

		Grid of (scalar) Electric field values excluding guard cells

		Returns
		-------
		E : numpy.array, (nx)
			Electric field values.
		"""
		cdef float *buf = <float *> self._thisptr.E.s
		return np.asarray( <float [:self._thisptr.E.nx]> buf )


cdef class Charge:
	"""Charge()
	
	Electric charge density class

	This class allows access to the electric chrage density data
	structures in the simulation. An object of this class is created
	automatically when creating a `es1d.Simulation` object.

	See Also
	--------
	es1d.Simulation
	"""

	cdef t_charge* _thisptr

	cdef associate( self, t_charge* ptr ):
		self._thisptr = ptr

	def report( self ):
		"""report( )

		Save charge density information to disk. Files will be saved in
		the CHARGE directory below the current working directory.
		"""
		charge_report( self._thisptr )

	@property
	def rho( self ):
		"""Charge density

		Grid of charge density values excluding guard cells

		Returns
		-------
		rho : numpy.array, (nx)
			Charge density values
		"""
		cdef float *buf = <float *> self._thisptr.rho.s
		return np.asarray( <float [:self._thisptr.rho.nx]> buf )


cdef class Simulation:
	"""Simulation(nx, box, dt, species=None, report=None, init_fld=None, ext_fld=None, neutral_bkg=False)
	
	ZPIC ES1D Simulation class

	Parameters
	----------
	nx : int
		Number of grid cells for the simulation
	box : float
		Simulation box (phyiscal) size, in simulation units
	dt : float
		Simulation time step, in simulation units
	species : `Species` or list of `Species`, optional
		Particle species to use in the simulation, defaults to None
		(no particles)
	report : function, optional
		Python function used for simulation reporting, defaults to None
	neutral_bkg : `bool`, optional
		Controls adding a neutralizing charge background to the
		simulation at initialization, defaults to False.
	
	See also
	--------
	es1d.Species
	"""

	cdef t_simulation *_thisptr
	cdef int n
	cdef float t

	cdef Field field
	cdef Charge charge
	cdef list species

	cdef object report

	def __cinit__( self, int nx, float box, float dt, *, species = None,
	               report = None, bint neutral_bkg = False ):

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

		self.field = Field()
		self.field.associate( &self._thisptr.field )

		self.charge = Charge()
		self.charge.associate( &self._thisptr.charge )

		# Neutralizing background
		if ( neutral_bkg ):
			sim_add_neutral_bkg( self._thisptr )

	def __dealloc__(self):
		sim_delete( self._thisptr )
		free(self._thisptr)

	def add_neutral_bkg( self ):
		"""add_neutral_bkg()
		
		Adds a (initial) neutralizing charge background to the simulation

		Note
		----
		Use of this function has been deprecated and will be removed soon.
		Use the `neutral_bkg` parameter of the `Simulation` class instead.
		"""
		sim_add_neutral_bkg( self._thisptr )

	def iter( self ):
		"""iter()

		Advance simulation 1 iteration.
		"""
		sim_iter( self._thisptr )
		self.n = self.n+1
		self.t = self.n * self._thisptr.dt

	def run( self, float tmax ):
		"""run(tmax)

		Advance simulation up to time `tmax`. If specified earlier, the
		`report` function will be called before each iteration.

		Parameters
		----------
		tmax : float
			Intended final simulation time. If smaller than current
			simulation time, a warning message will be displayed
		"""
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
	def species(self):
		"""Simulation particle species list

		Returns
		-------
		species : list of `es1d.Species`
			Simulation particle species list
		"""
		return self.species

	@property
	def field(self):
		"""Simulation Electric field object

		Returns
		-------
		field : `es1d.Field`
			Simulation Electric field object
		"""
		return self.field

	@property
	def charge(self):
		"""Simulation electric charge density object

		Returns
		-------
		current : `es1d.Charge`
			Simulation electric charge density object
		"""
		return self.charge

	@property
	def n(self):
		"""Current simulation iteration number

		This number is advanced automatically by calls to the `iter()`
		and `run()` methods

		Returns
		-------
		n : int
			Current simulation iteration number
		"""
		return self.n

	@property
	def t(self):
		"""Current simulation time value

		This value is advanced automatically by calls to the `iter()`
		and `run()` methods

		Returns
		-------
		t : float
			Current simulation simulation time
		"""
		return self.t

	@property
	def dx(self):
		"""Cell size used for the simulation

		Returns
		-------
		dx : float
			Cell size
		"""
		return self.field.dx

	@property
	def nx(self):
		"""Grid size used for the simulation

		Returns
		-------
		nx : int
			Number of grid cells for the simulation
		"""
		return self.field.nx

	@property
	def box(self):
		"""Simulation box physical size

		Returns
		-------
		box : float
			Simulation box size
		"""
		return self.field.box

	@property
	def report(self):
		"""Report function for the simulation

		This function will be called once before each time step when using
		the `run()` method.

		Returns
		-------
		report : function
			Report function for the simulation
		"""
		return self.report

	@report.setter
	def report( self, f ):
		self.report = f





