"""# 1D spectral EM-PIC code

1D, electro-magnetic, fully relativistic, Particle-in-Cell code, using
a spectral EM field solver
"""
#cython: language_level=3

cimport em1ds
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
	em1ds.Species
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
	"""Species(name, m_q, ppc, ufl=[0.,0.,0.], uth=[0.,0.,0.], density=None, n_sort=16)
	
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
	ufl : list, optional
		Initial fluid (generalized) velocity for the particles, 
		defaults to [0,0,0] (no fluid velocity)
	uth : list, optional
		Initial thermal velocity for the particles, defaults to [0,0,0]
		(no thermal velocity).
	density : `Density`, optional
		Density profile for the particle species specified as using a
		`Density` object. Defaults to `None` which corresponds to a
		uniform density of value 1
	n_sort : int, optional
		Number of iterations between particle sort, defaults to 16
	
	See also
	--------
	em1ds.Density
	"""

	cdef t_species _this
	cdef t_species* _thisptr
	cdef Density _density
	cdef str _name

	# Diagnostic types
	_diag_types  = { 'charge':CHARGE, 'pha':PHA, 'particles':PARTICLES }
	_pha_quants = { 'x1':X1, 'u1':U1, 'u2':U2, 'u3':U3 }

	def __cinit__( self, str name, const float m_q, const int ppc, *,
				  list ufl = [0.,0.,0.], list uth = [0.,0.,0.], Density density = None,
				  int n_sort = 16):

		self._thisptr = &self._this
		self._name = name
		self._this.m_q = m_q
		self._this.ppc = ppc
		self._this.ufl = np.array(ufl, dtype=np.float32)
		self._this.uth = np.array(uth, dtype=np.float32)
		self._this.n_sort = n_sort

		if ( density ):
			self._density = density.copy()
		else:
			# Use default uniform density
			self._density = Density()

	cdef new( self, t_species* ptr, int nx, float box, float dt ):
		self._thisptr = ptr
		n_sort = self._this.n_sort
		spec_new( self._thisptr, self._name.encode(), self._this.m_q, self._this.ppc,
			self._this.ufl, self._this.uth,
			nx, box, dt, self._density._thisptr )
		self._this.n_sort = n_sort

	def add( self, int ix, float x, float[:] u):
		"""add(ix, x, u)

		Adds a new particle to the particle buffer

		Parameters
		----------
		ix : int
			New particle cell index
		x : float
			New particle position inside the cell
		u : float[3]
			New particle (generalized velocity)
		"""
		# insure we have enough room for new particle
		spec_grow_buffer( self._thisptr, self._thisptr.np + 1 )
		
		cdef t_part particle
		particle.ix = ix
		particle.x = x
		particle.ux = u[0]
		particle.uy = u[1]
		particle.uz = u[2]

		self._thisptr.part[ self._thisptr.np ] = particle
		self._thisptr.np = self._thisptr.np + 1


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
			Each quantity must be one of 'x1', 'u1', 'u2' or 'u3'
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

	def phasespace( self, list quants, pha_nx, pha_range ):
		"""phasespace(quants, pha_nx, pha_range)
		
		Calculate phasespace density of particle species

		Parameters
		----------
		quants : list of str
			2 element list of quantities to use for "pha" diagnostics.
			Each quantity must be one of 'x1', 'u1', 'u2' or 'u3'
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
		cdef int rep_type = PHASESPACE( self._pha_quants[quants[0]],
				                        self._pha_quants[quants[1]])

		_nx = np.array( pha_nx, dtype = np.int32)
		_range = np.array( pha_range, dtype = np.float32)

		pha = np.zeros( shape = (_nx[1],_nx[0]), dtype = np.float32 )
		cdef float [:,:] buf = pha

		spec_deposit_pha( self._thisptr, rep_type, _nx, _range, &buf[0,0] )

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

cdef float3 custom_ext_E( int ix, float dx, void *f ):
	"""Internal auxiliary function for using external E fields"""
	cdef ExternalField ext = <object> f
	val = ext.custom_func_E(ix,dx)
	cdef float3 e
	e.x = val[0]
	e.y = val[1]
	e.z = val[2]
	return e

cdef float3 custom_ext_B( int ix, float dx, void *f ):
	"""Internal auxiliary function for using external B fields"""
	cdef ExternalField ext = <object> f
	val = ext.custom_func_B(ix,dx)
	cdef float3 b
	b.x = val[0]
	b.y = val[1]
	b.z = val[2]
	return b

cdef class ExternalField:
	"""ExternalField( E_type='none', B_type='none', E_0=[0.,0.,0.], B_0 = [0.,0.,0.], E_custom=None, B_custom=None)
	
	Used for defining external EM fields in the simulation

	Parameters
	----------
	E_type : {'none','uniform','custom'}, optional
		Type of external electric field to use, must be one of 'none'
		(no external field), 'uniform' (uniform external field) or
		'custom' (custom external field defined by a function),
		defaults to 'none'
	B_type : {'none','uniform','custom'}, optional
		Type of magnetic electric field to use, must be one of 'none'
		(no external field), 'uniform' (uniform external field) or
		'custom' (custom external field defined by a function), 
		defaults to 'none'
	E_0 : list of float, optional
		3 element list specifying electric field value for 'uniform'
		external field type, defaults to [0.,0.,0.]
	B_0 : list of float, optional
		3 element list specifying magnetic field value for 'uniform'
		external field type, defaults to [0.,0.,0.]
	E_custom : function`
		Python function for calculating all 3 components of external
		electric field at every cell
	B_custom : function`
		Python function for calculating all 3 components of external
		magnetic field at every cell
	"""
	
	cdef t_emf_ext_fld *_thisptr

	cdef object custom_func_E
	cdef object custom_func_B

	_ext_types = {'none':EMF_FLD_TYPE_NONE, 'uniform':EMF_FLD_TYPE_UNIFORM, 'custom':EMF_FLD_TYPE_CUSTOM}

	def __cinit__( self, *, str E_type = 'none', str B_type = 'none', 
				list E_0 = [0.,0.,0.], list B_0 = [0.,0.,0.],
				E_custom = None, B_custom = None ):

		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_emf_ext_fld *> calloc(1, sizeof(t_emf_ext_fld))

		self._thisptr.E_type = <emf_fld_type> self._ext_types[E_type]
		self._thisptr.B_type = <emf_fld_type> self._ext_types[B_type]

		buf = np.array( E_0, dtype=np.float32)
		self._thisptr.E_0.x = buf[0]
		self._thisptr.E_0.y = buf[1]
		self._thisptr.E_0.z = buf[2]

		buf = np.array( B_0, dtype=np.float32)
		self._thisptr.B_0.x = buf[0]
		self._thisptr.B_0.y = buf[1]
		self._thisptr.B_0.z = buf[2]

		if ( E_custom ):
			self.custom_func_E = E_custom
			self._thisptr.E_custom = custom_ext_E
			self._thisptr.E_custom_data = <void *> self
		if ( B_custom ):
			self.custom_func_B = B_custom
			self._thisptr.B_custom = custom_ext_B
			self._thisptr.B_custom_data = <void *> self


	def __dealloc__(self):
		free( self._thisptr )

	def copy(self):
		"""copy()

		Object copy.
		"""
		new = ExternalField()
		new.E_type  = self.E_type
		new.B_type  = self.B_type
		new.E_0	    = self.E_0
		new.B_0	    = self.B_0
		new.custom_func_E = self.custom_func_E
		new.custom_func_B = self.custom_func_B
		new._thisptr.E_custom = self._thisptr.E_custom
		new._thisptr.B_custom = self._thisptr.B_custom
		new._thisptr.E_custom_data = self._thisptr.E_custom_data
		new._thisptr.B_custom_data = self._thisptr.B_custom_data

		return new

	@property
	def E_type(self):
		"""Type of external E field

		Returns
		-------
		type : {'none','uniform','custom'}
			Type of external field
		"""
		tmp = {EMF_FLD_TYPE_NONE:'none', EMF_FLD_TYPE_UNIFORM:'uniform', EMF_FLD_TYPE_CUSTOM:'custom'}
		return tmp[self._thisptr.E_type]

	@E_type.setter
	def E_type(self,value):
		self._thisptr.E_type = self._init_types[value]

	@property
	def B_type(self):
		"""Type of external B field

		Returns
		-------
		type : {'none','uniform','custom'}
			Type of external field
		"""
		tmp = {EMF_FLD_TYPE_NONE:'none', EMF_FLD_TYPE_UNIFORM:'uniform', EMF_FLD_TYPE_CUSTOM:'custom'}
		return tmp[self._thisptr.B_type]

	@B_type.setter
	def B_type(self,value):
		self._thisptr.B_type = self._init_types[value]

	@property
	def E_0(self):
		"""Electric field value for 'uniform' external field type

		Returns
		-------
		E0 : numpy.ndarray, (3)
			Vector field value for 'uniform' external field type
		"""
		return self._thisptr.E_0

	@E_0.setter
	def E_0(self,value):
		self._thisptr.E_0 = value

	@property
	def B_0(self):
		"""Magnetic field value for 'uniform' external field type

		Returns
		-------
		B0 : numpy.ndarray, (3)
			Vector field value for 'uniform' external field type
		"""
		return self._thisptr.E_0

	@B_0.setter
	def B_0(self,value):
		self._thisptr.E_0 = value

	@property
	def E_custom(self):
		"""Python function used for 'custom' external E field type

		Returns
		-------
		func : function
			Python function for calculating all 3 components of external
			electric field at every cell
		"""
		return self.custom_func_E

	@E_custom.setter
	def E_custom(self,value):
		self.custom_func_E = value
		self._thisptr.E_custom = custom_ext_E
		self._thisptr.E_custom_data = <void *> self

	@property
	def B_custom(self):
		"""Python function used for 'custom' external B field type

		Returns
		-------
		func : function
			Python function for calculating all 3 components of external
			magnetic field at every cell
		"""
		return self.custom_func_B

	@B_custom.setter
	def B_custom(self,value):
		self.custom_func_B = value
		self._thisptr.B_custom = custom_ext_B
		self._thisptr.B_custom_data = <void *> self


cdef float3 custom_init_E( int ix, float dx, void *f ):
	"""Internal auxiliary function for using initial E fields"""
	cdef InitialField init = <object> f
	val = init.custom_func_E(ix,dx)
	cdef float3 e
	e.x = val[0]
	e.y = val[1]
	e.z = val[2]
	return e

cdef float3 custom_init_B( int ix, float dx, void *f ):
	"""Internal auxiliary function for using initial B fields"""
	cdef InitialField init = <object> f
	val = init.custom_func_B(ix,dx)
	cdef float3 b
	b.x = val[0]
	b.y = val[1]
	b.z = val[2]
	return b


cdef class InitialField:
	"""InitialField( E_type='none', B_type='none', E_0=[0.,0.,0.], B_0 = [0.,0.,0.], E_custom=None, B_custom=None)
	
	Used for defining initial EM fields in the simulation

	Parameters
	----------
	E_type : {'none','uniform','custom'}, optional
		Type of initial electric field to use, must be one of 'none'
		(no initial field), 'uniform' (uniform initial field) or
		'custom' (custom initial field defined by a function),
		defaults to 'none'
	B_type : {'none','uniform','custom'}, optional
		Type of initial electric field to use, must be one of 'none'
		(no initial field), 'uniform' (uniform initial field) or
		'custom' (custom initial field defined by a function), 
		defaults to 'none'
	E_0 : list of float, optional
		3 element list specifying electric field value for 'uniform'
		initial field type, defaults to [0.,0.,0.]
	B_0 : list of float, optional
		3 element list specifying magnetic field value for 'uniform'
		initial field type, defaults to [0.,0.,0.]
	E_custom : function`
		Python function for calculating all 3 components of initial
		electric field at every cell
	B_custom : function`
		Python function for calculating all 3 components of initial
		magnetic field at every cell
	"""
	cdef t_emf_init_fld *_thisptr

	cdef object custom_func_E
	cdef object custom_func_B

	_init_types = {'none':EMF_FLD_TYPE_NONE, 'uniform':EMF_FLD_TYPE_UNIFORM, 'custom':EMF_FLD_TYPE_CUSTOM}

	def __cinit__( self, *, str E_type = 'none', str B_type = 'none', 
				list E_0 = [0.,0.,0.], list B_0 = [0.,0.,0.],
				E_custom = None, B_custom = None ):

		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_emf_init_fld *> calloc(1, sizeof(t_emf_init_fld))

		self._thisptr.E_type = <emf_fld_type> self._init_types[E_type]
		self._thisptr.B_type = <emf_fld_type> self._init_types[B_type]

		buf = np.array( E_0, dtype=np.float32)
		self._thisptr.E_0.x = buf[0]
		self._thisptr.E_0.y = buf[1]
		self._thisptr.E_0.z = buf[2]

		buf = np.array( B_0, dtype=np.float32)
		self._thisptr.B_0.x = buf[0]
		self._thisptr.B_0.y = buf[1]
		self._thisptr.B_0.z = buf[2]

		if ( E_custom ):
			self.custom_func_E = E_custom
			self._thisptr.E_custom = custom_init_E
			self._thisptr.E_custom_data = <void *> self
		if ( B_custom ):
			self.custom_func_B = B_custom
			self._thisptr.B_custom = custom_init_B
			self._thisptr.B_custom_data = <void *> self


	def __dealloc__(self):
		free( self._thisptr )

	def copy(self):
		"""copy()

		Object copy.
		"""
		new = InitialField()
		new.E_type  = self.E_type
		new.B_type  = self.B_type
		new.E_0	    = self.E_0
		new.B_0	    = self.B_0
		new.custom_func_E = self.custom_func_E
		new.custom_func_B = self.custom_func_B
		new._thisptr.E_custom = self._thisptr.E_custom
		new._thisptr.B_custom = self._thisptr.B_custom
		new._thisptr.E_custom_data = self._thisptr.E_custom_data
		new._thisptr.B_custom_data = self._thisptr.B_custom_data

		return new

	@property
	def E_type(self):
		"""Type of initial E field

		Returns
		-------
		type : {'none','uniform','custom'}
			Type of initial field
		"""
		tmp = {EMF_FLD_TYPE_NONE:'none', EMF_FLD_TYPE_UNIFORM:'uniform', EMF_FLD_TYPE_CUSTOM:'custom'}
		return tmp[self._thisptr.E_type]

	@E_type.setter
	def E_type(self,value):
		self._thisptr.E_type = self._init_types[value]

	@property
	def B_type(self):
		"""Type of initial B field

		Returns
		-------
		type : {'none','uniform','custom'}
			Type of initial field
		"""
		tmp = {EMF_FLD_TYPE_NONE:'none', EMF_FLD_TYPE_UNIFORM:'uniform', EMF_FLD_TYPE_CUSTOM:'custom'}
		return tmp[self._thisptr.B_type]

	@B_type.setter
	def B_type(self,value):
		self._thisptr.B_type = self._init_types[value]

	@property
	def E_0(self):
		"""Electric field value for 'uniform' initial field type

		Returns
		-------
		E0 : numpy.ndarray, (3)
			Vector field value for 'uniform' initial field type
		"""
		return self._thisptr.E_0

	@E_0.setter
	def E_0(self,value):
		self._thisptr.E_0 = value

	@property
	def B_0(self):
		"""Magnetic field value for 'uniform' initial field type

		Returns
		-------
		E0 : numpy.ndarray, (3)
			Vector field value for 'uniform' initial field type
		"""
		return self._thisptr.E_0

	@B_0.setter
	def B_0(self,value):
		self._thisptr.E_0 = value

	@property
	def E_custom(self):
		return self.custom_func_E

	@E_custom.setter
	def E_custom(self,value):
		"""Python function used for 'custom' initial E field type

		Returns
		-------
		func : function
			Python function for calculating all 3 components of initial
			electric field at every cell
		"""
		self.custom_func_E = value
		self._thisptr.E_custom = custom_init_E
		self._thisptr.E_custom_data = <void *> self

	@property
	def B_custom(self):
		"""Python function used for 'custom' initial B field type

		Returns
		-------
		func : function
			Python function for calculating all 3 components of initial
			magnetic field at every cell
		"""
		return self.custom_func_B

	@B_custom.setter
	def B_custom(self,value):
		self.custom_func_B = value
		self._thisptr.B_custom = custom_init_B
		self._thisptr.B_custom_data = <void *> self


cdef class EMF:
	"""EMF()
	
	Electro-Magnetic fields class

	This class allows access to the EM field data structures in the
	simulation. An object of this class is created automatically when
	creating a `em1ds.Simulation` object.

	See Also
	--------
	em1ds.Simulation
	"""

	cdef t_emf* _thisptr

	# Diagnostic types
	_diag_types = { 'E' : EFLD,	'B' : BFLD }

	# Field solver types
	_solver_types = {'PSTD' : EMF_SOLVER_PSTD,
                     'PSATD': EMF_SOLVER_PSATD}

	cdef associate( self, t_emf* ptr ):
		self._thisptr = ptr

	def report( self, str type, char fc ):
		"""report( type, fc )

		Save diagnostic information to disk. Files will be saved in the
		EMF directory below the current working directory.

		Parameters
		----------
		type : {'E','B'}
			Type of data to save, must be one of 'E' (electric field),
			'B' (magnetic field)
		fc : {0,1,2}
			Field component to save, must be one of 0 (x), 1 (y) or 2 (z)
		"""
		cdef int rep_type = self._diag_types[type];
		emf_report( self._thisptr, rep_type, fc )

	def set_ext_fld(self, ExternalField ext):
		"""set_ext_fld( ext )

		Sets external EM field values. This method can only be called
		before the simulation starts.

		Parameters
		----------
		ext : `ExternalField`
			External field parameters

		See Also
		--------
		em2d.ExternalField

		Note
		----
		Use of this function has been deprecated and will be removed soon.
		Use the `ext_fld` parameter of the `Simulation` class instead.
		"""
		if (self._thisptr.iter == 0):
			emf_set_ext_fld( self._thisptr, ext._thisptr )
		else:
			print("set_ext_fld can only be called before the simulation starts")

	def init_fld(self, InitialField init):
		"""init_fld( init )

		Sets initial EM field values. This method can only be called
		before the simulation starts.

		Parameters
		----------
		init : `InitialField`
			Initial field parameters

		See Also
		--------
		em2d.InitialField

		Note
		----
		Use of this function has been deprecated and will be removed soon.
		Use the `init_fld` parameter of the `Simulation` class instead.
		"""
		if (self._thisptr.iter == 0):
			emf_init_fld( self._thisptr, init._thisptr )
		else:
			print("init_fld can only be called before the simulation starts")

	@property
	def energy( self ):
		"""EM field energy per field component

		Returns
		-------
		energy : numpy.array, (6)
			Total energy in each field component. Electric field energy is
			in the first 3 values and magnetic field energy is in the last
			3 values.
		
		Note
		----
		These values are recalculated each time this function is called.
		"""
		cdef double energy[6]
		emf_get_energy( self._thisptr, energy )
		return np.array( energy, dtype = np.float64 )

	@property
	def nx(self):
		"""Grid size used for the EMF object

		Returns
		-------
		nx :int
			Number of grid cells for the simulation
		"""
		return self._thisptr.E.nx

	@property
	def dx(self):
		"""Cell size used for the EMF object

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
	def solver_type(self):
		"""Field solver algorithm used

		Returns
		-------
		solver : 'PSTD', 'PSATD'
			Field solver in use, either 'PSTD' (pseudo-spectral time
			domain) or 'PSATD' (pseudo-spectral analytical time domain)
		"""
		for key, value in self._solver_types.items():
			if ( value == self._thisptr.solver_type ):
				return key
		return 'unknown'

	@solver_type.setter
	def solver_type( self, str solver ):
		self._thisptr.solver_type = self._solver_types[solver]

	@property
	def Ex( self ):
		"""Ex field component

		Grid of (scalar) Ex field values excluding guard cells

		Returns
		-------
		Ex : numpy.array, (nx)
			X component values of the electric field.
		"""
		cdef float *buf = <float *> self._thisptr.E.x - self._thisptr.E.gc[0]
		cdef int size = self._thisptr.E.gc[0] + self._thisptr.E.nx + self._thisptr.E.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.E.gc[0] : self._thisptr.E.gc[0] + self._thisptr.E.nx ]

	@property
	def Ey( self ):
		"""Ey field component

		Grid of (scalar) Ey field values excluding guard cells

		Returns
		-------
		Ey : numpy.array, (nx)
			Y component values of the electric field.
		"""
		cdef float *buf = <float *> self._thisptr.E.y - self._thisptr.E.gc[0]
		cdef int size = self._thisptr.E.gc[0] + self._thisptr.E.nx + self._thisptr.E.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.E.gc[0] : self._thisptr.E.gc[0] + self._thisptr.E.nx ]

	@property
	def Ez( self ):
		"""Ez field component

		Grid of (scalar) Ey field values excluding guard cells

		Returns
		-------
		Ez : numpy.array, (nx)
			Z component values of the electric field.
		"""
		cdef float *buf = <float *> self._thisptr.E.z - self._thisptr.E.gc[0]
		cdef int size = self._thisptr.E.gc[0] + self._thisptr.E.nx + self._thisptr.E.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.E.gc[0] : self._thisptr.E.gc[0] + self._thisptr.E.nx ]


	@property
	def Bx( self ):
		"""Bx field component

		Grid of (scalar) Bx field values excluding guard cells

		Returns
		-------
		Bx : numpy.array, (nx)
			X component values of the magnetic field.
		"""
		cdef float *buf = <float *> self._thisptr.B.x - self._thisptr.B.gc[0]
		cdef int size = self._thisptr.B.gc[0] + self._thisptr.B.nx + self._thisptr.B.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.B.gc[0] : self._thisptr.B.gc[0] + self._thisptr.B.nx ]

	@property
	def By( self ):
		"""By field component

		Grid of (scalar) By field values excluding guard cells

		Returns
		-------
		By : numpy.array, (nx)
			Y component values of the magnetic field.
		"""
		cdef float *buf = <float *> self._thisptr.B.y - self._thisptr.B.gc[0]
		cdef int size = self._thisptr.B.gc[0] + self._thisptr.B.nx + self._thisptr.B.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.B.gc[0] : self._thisptr.B.gc[0] + self._thisptr.B.nx ]

	@property
	def Bz( self ):
		"""Bz field component

		Grid of (scalar) Bz field values excluding guard cells

		Returns
		-------
		Bz : numpy.array, (nx)
			Z component values of the magnetic field.
		"""
		cdef float *buf = <float *> self._thisptr.B.z - self._thisptr.B.gc[0]
		cdef int size = self._thisptr.B.gc[0] + self._thisptr.B.nx + self._thisptr.B.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.B.gc[0] : self._thisptr.B.gc[0] + self._thisptr.B.nx ]


cdef class Laser:
	"""Laser(start=0.0, fwhm=0.0, rise=0.0, flat=0.0, fall=0.0, a0=0.0, omega0=0.0, polarization=0.0)
	
	Class representing laser pulses. Laser pulses are added to the
	simulation using the `Simulation.add_laser()` method.

	Parameters
	----------
	start : float, optional
		Position of the starting (leading) point for the laser envelope,
		defaults to 0.
	fwhm : float, optional
		Full width at half-max of the laser pulse. If set it overrides the
		`rise`, `flat`, and `fall` parameters, defaults to 0
	rise, flat, fall : float, optional
		Rise time (`rise`), flat time (`flat`) and fall time (`fall`) of
		the temporal envelope, default to 0
	a0 : float, optional
		Normalized vector potential value at peak intensity of the laser
		pulse, default to 0
	omega0 : float, optional
		Laser frequency in simulation units, defaults to 0
	polarization : float, optional
		Laser polarization in radians measured in reference to the y
		direction, defaults to 0

	See Also
	--------
	em1ds.Simulation
	"""
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
		"""Start position of laser envelope

		Returns
		-------
		start : float
			Position of the starting (leading) point for the laser envelope
		"""
		return self._thisptr.start

	@start.setter
	def start(self,value):
		self._thisptr.start = value

	@property
	def fwhm(self):
		return self._thisptr.fwhm

	@fwhm.setter
	def fwhm(self,value):
		"""FWHM of laser envelope
		
		Returns
		-------
		fwhm : float
			Full width at half-max of the laser pulse
		"""
		self._thisptr.fwhm = value

	@property
	def rise(self):
		"""Rise length of laser envelope
		
		Returns
		-------
		rise : float
			Rise length of laser envelope
		"""
		return self._thisptr.rise

	@rise.setter
	def rise(self,value):
		self._thisptr.rise = value

	@property
	def flat(self):
		"""Flat length of laser envelope
		
		Returns
		-------
		flat : float
			flat length of laser envelope
		"""
		return self._thisptr.flat

	@flat.setter
	def flat(self,value):
		self._thisptr.flat = value

	@property
	def fall(self):
		"""Fall length of laser envelope
		
		Returns
		-------
		fall : float
			fall length of laser envelope
		"""
		return self._thisptr.fall

	@fall.setter
	def fall(self,value):
		self._thisptr.fall = value

	@property
	def a0(self):
		"""Normalized vector potential value at peak intensity of the
		laser pulse
		
		Returns
		-------
		a0 : float
			Normalized vector potential value
		"""
		return self._thisptr.a0

	@a0.setter
	def a0(self,value):
		self._thisptr.a0 = value

	@property
	def omega0(self):
		"""Laser frequency in simulation units
		
		Returns
		-------
		omega0 : float
			Laser frequency
		"""
		return self._thisptr.omega0

	@omega0.setter
	def omega0(self,value):
		self._thisptr.omega0 = value

	@property
	def polarization(self):
		"""Laser polarization in radians measured in reference to the y
		direction
		
		Returns
		-------
		pol : float
			Polarization angle
		"""
		return self._thisptr.polarization

	@polarization.setter
	def polarization(self,value):
		self._thisptr.polarization = value


cdef class Current:
	"""Current()
	
	Electric current density class

	This class allows access to the electric current density data
	structures in the simulation. An object of this class is created
	automatically when creating a `em1ds.Simulation` object.

	See Also
	--------
	em1ds.Simulation
	"""

	cdef t_current* _thisptr

	cdef associate( self, t_current* ptr ):
		self._thisptr = ptr

	def report( self, char jc ):
		"""report( jc )

		Save diagnostic information to disk. Files will be saved in the
		CURRENT directory below the current working directory.

		Parameters
		----------
		jc : {0,1,2}
			Current density component to save, must be one of 0 (x), 1 (y)
			or 2 (z)
		"""
		current_report( self._thisptr, jc )

	@property
	def Jx( self ):
		"""Jx current density component

		Grid of (scalar) Jx field values excluding guard cells

		Returns
		-------
		Jx : numpy.array, (nx)
			X component values of the current density
		"""
		cdef float *buf = <float *> self._thisptr.J.x - self._thisptr.J.gc[0]
		cdef int size = self._thisptr.J.gc[0] + self._thisptr.J.nx + self._thisptr.J.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.J.gc[0] : self._thisptr.J.gc[0] + self._thisptr.J.nx ]

	@property
	def Jy( self ):
		"""Jy current density component

		Grid of (scalar) Jy field values excluding guard cells

		Returns
		-------
		Jy : numpy.array, (nx)
			Y component values of the current density
		"""
		cdef float *buf = <float *> self._thisptr.J.y - self._thisptr.J.gc[0]
		cdef int size = self._thisptr.J.gc[0] + self._thisptr.J.nx + self._thisptr.J.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.J.gc[0] : self._thisptr.J.gc[0] + self._thisptr.J.nx ]

	@property
	def Jz( self ):
		"""Jz current density component

		Grid of (scalar) Jz field values excluding guard cells

		Returns
		-------
		Jz : numpy.array, (nx)
			Z component values of the current density
		"""
		cdef float *buf = <float *> self._thisptr.J.z - self._thisptr.J.gc[0]
		cdef int size = self._thisptr.J.gc[0] + self._thisptr.J.nx + self._thisptr.J.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.J.gc[0] : self._thisptr.J.gc[0] + self._thisptr.J.nx ]

cdef class Charge:
	"""Charge()
	
	Electric charge density class

	This class allows access to the electric chrage density data
	structures in the simulation. An object of this class is created
	automatically when creating a `em1ds.Simulation` object.

	See Also
	--------
	em1ds.Simulation
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
		cdef float *buf = <float *> self._thisptr.rho.buffer
		cdef int size = self._thisptr.rho.gc[0] + self._thisptr.rho.nx + self._thisptr.rho.gc[1]
		tmp = np.asarray( <float [:size]> buf )
		return tmp[ self._thisptr.rho.gc[0] : self._thisptr.rho.gc[0] + self._thisptr.rho.nx ]


cdef class Simulation:
	"""Simulation(nx, box, dt, species=None, report=None, init_fld=None, ext_fld=None, neutral_bkg=False)
	
	ZPIC EM1DS Simulation class

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
	init_fld : `InitialField`, optional
		Initial EM fields for the simulation, defaults to None (0 initial
		fields)
	ext_fld : `ExternalField`, optional
		External EM fields for the simulation, defaults to None
	neutral_bkg : `bool`, optional
		Controls adding a neutralizing charge background to the
		simulation at initialization, defaults to False.
	
	See also
	--------
	em1ds.Species
	em1ds.InitialField
	em1ds.ExternalField
	"""

	cdef t_simulation *_thisptr
	cdef int n
	cdef float t

	cdef EMF emf
	cdef Current current
	cdef Charge charge
	cdef list species

	cdef object report

    # Filter types
	_filter_types = {'none' :     FILTER_NONE,
                     'gaussian' : FILTER_GAUSS,
                     'sharp' :    FILTER_SHARP}

	def __cinit__( self, int nx, float box, float dt, *, species = None, report = None,
				   InitialField init_fld = None, ExternalField ext_fld = None,
				   bint neutral_bkg = False ):

		# Sanity checks
		if ( nx < 2 ):
			print("Invalid number of cells", file = sys.stderr)
			return

		if ( box <= 0 ):
			print("Invalid box size, must be > 0", file = sys.stderr)
			return

		if ( dt < 0 ):
			print("Invalid time-step, must be > 0", file = sys.stderr)
			return

		if ( dt >= box/nx ):
			print("Invalid timestep (courant condition violation), dt must be < {:g}".format( box/nx ) , file = sys.stderr)
			return

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

		self.charge = Charge()
		self.charge.associate( &self._thisptr.charge )

		# Set initial fields
		if ( init_fld ):
			emf_init_fld( self.emf._thisptr, init_fld._thisptr )
		
		# Set external fields
		if ( ext_fld ):
			emf_set_ext_fld( self.emf._thisptr, ext_fld._thisptr )
		
		# Neutralizing background
		if ( neutral_bkg ):
			sim_add_neutral_bkg( self._thisptr )

	def __dealloc__(self):
		sim_delete( self._thisptr )
		free(self._thisptr)

	def add_laser(self, Laser laser):
		"""add_laser(laser)
		
		Adds laser pulse to the simulation

		Parameters
		----------
		laser : `Laser`
			Laser pulse to be added to the simulation
		
		See Also
		--------
		em1ds.Laser
		"""
		sim_add_laser( self._thisptr, laser._thisptr )

	def add_neutral_bkg(self):
		"""add_neutral_bkg()
		
		Adds a (initial) neutralizing charge background to the simulation

		Note
		----
		Use of this function has been deprecated and will be removed soon.
		Use the `neutral_bkg` parameter of the `Simulation` class instead.
		"""
		sim_add_neutral_bkg( self._thisptr )

	def filter_set( self, str type, *, float ck = 0.0 ):
		"""filter_set(type, ck = 0.0)

		Sets spectral filtering parameters

		Parameters
		----------
		type : {'none', 'gaussian', 'sharp'}
			Type of spectral filtering to use, must be one of 'none' (no
			filtering), 'gaussian' (gaussian shaped transfer function), or
			'sharp' (perfect low pass filter)
		ck : float
			Filter parameter for 'gaussian' and 'sharp' filters, for
			'sharp' filtering this value must be in the ]0.0,1.0[ range
			(cutoff frequecny in units of the Nyquist frequency) and for
			'gaussian' filtering, this value must be > 0.0. 
		"""
		cdef int filter_type;

		filter_type = self._filter_types[type]

		if ( filter_type == FILTER_SHARP ):
			if ( ck <= 0.0 or ck >= 1.0 ):
				print("For sharp filter ck must be in the ]0.0,1.0[ range", file = sys.stderr)
				return
		elif ( filter_type == FILTER_GAUSS ):
			if ( ck <= 0.0 ):
				print("For gaussian filter ck must be > 0.0", file = sys.stderr)
				return

		sim_filter_set( self._thisptr, filter_type, ck )

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
		species : list of `em1ds.Species`
			Simulation particle species list
		"""
		return self.species

	@property
	def emf(self):
		"""Simulation EM fields object

		Returns
		-------
		emf : `em1ds.EMF`
			Simulation EM fields object
		"""
		return self.emf

	@property
	def current(self):
		"""Simulation electric current density object

		Returns
		-------
		current : `em1ds.Current`
			Simulation electric current density object
		"""
		return self.current

	@property
	def charge(self):
		"""Simulation electric charge density object

		Returns
		-------
		current : `em1ds.Charge`
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
	def dt(self):
		"""Time step used for the simulation

		Returns
		-------
		dt :float
			Time step
		"""
		return self._thisptr.dt

	@property
	def dx(self):
		"""Cell size used for the simulation

		Returns
		-------
		dx : float
			Cell size
		"""
		return self.emf.dx

	@property
	def nx(self):
		"""Grid size used for the simulation

		Returns
		-------
		nx : int
			Number of grid cells for the simulation
		"""
		return self.emf.nx

	@property
	def box(self):
		"""Simulation box physical size

		Returns
		-------
		box : float
			Simulation box size
		"""
		return self.emf.box

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





