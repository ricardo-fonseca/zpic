"""# 2D EM-PIC code

2D, electro-magnetic, fully relativistic, Particle-in-Cell code, using
a finite-difference EM field solver
"""
#cython: language_level=3

cimport em2d
from libc.stdlib cimport calloc, free

import numpy as np
import sys

cdef float custom_density_x( float x, void *f ):
	"""Internal function for custom density profiles"""
	cdef Density d = <object> f
	return d.custom_func_x(x)

cdef float custom_density_y( float y, void *f ):
	"""Internal function for custom density profiles"""
	cdef Density d = <object> f
	return d.custom_func_y(y)


cdef class Density:
	"""Density(type='uniform', start=0.0, end=0.0, n=1.0, custom_x=None, custom_y=None)
	
	Class representing charge density profiles for particle species
	initialization

	Parameters
	----------
	type : str, optional
		Density profile type, one of "uniform", "empty", "step", "slab"
		or "custom", defaults to "uniform"
	start : float, optional
		Position of the plasma start position for "step" or "slab"
		profiles, defaults to 0.0
	end : float, optional
		Position of the plasma end position for the "slab" profiles,
		defaults to 0.0
	n : float, optional
		Reference density to use, multiplies density profile value,
		defaults to 1.0
	custom_x : function, optional
		Custom density function along x, defaults to None
	custom_y : function, optional
		Custom density function along y, defaults to None
	
	See Also
	--------
	em2d.Species
	"""

	cdef t_density *_thisptr

	cdef object custom_func_x
	cdef object custom_func_y

	_density_types = {'uniform':UNIFORM,
					  'empty':EMPTY,
						'step':STEP,
						'slab':SLAB,
						'custom':CUSTOM}

	def __cinit__( self, *, str type = 'uniform', float start = 0.0, float end = 0.0,
				   float n = 1.0, custom_x = None, custom_y = None ):

		# Allocates the structure and initializes all elements to 0
		self._thisptr = <t_density *> calloc(1, sizeof(t_density))

		self._thisptr.type = <density_type> self._density_types[type]
		self._thisptr.n = n
		self._thisptr.start = start
		self._thisptr.end = end
		if ( custom_x ):
			self.custom_func_x = custom_x
			self._thisptr.custom_x = custom_density_x
			self._thisptr.custom_data_x = <void *> self
		if ( custom_y ):
			self.custom_func_y = custom_y
			self._thisptr.custom_y = custom_density_y
			self._thisptr.custom_data_y = <void *> self


	def __dealloc__(self):
		free( self._thisptr )

	def copy(self):
		"""copy()

		Object copy.
		"""
		new = Density()
		new.type  = self.type
		new.n	 = self.n
		new.start = self.start
		new.end   = self.end
		new.custom_func_x = self.custom_func_x
		new.custom_func_y = self.custom_func_y
		new._thisptr.custom_x = self._thisptr.custom_x
		new._thisptr.custom_y = self._thisptr.custom_y
		new._thisptr.custom_data_x = self._thisptr.custom_data_x
		new._thisptr.custom_data_y = self._thisptr.custom_data_y

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
		type : {'uniform', 'empty', 'step', 'slab', 'custom'}
			Density profile type
		"""
		tmp = {UNIFORM:'uniform', EMPTY:'empty', STEP:'step', SLAB:'slab', CUSTOM:'custom'}
		return tmp[self._thisptr.type]

	@type.setter
	def type(self,value):
		self._thisptr.type = self._density_types[value]

	@property
	def start(self):
		"""
		Position of the plasma start position for "step" or "slab" profiles
		
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
		Position of the plasma end position for "slab" profiles
		
		Returns
		-------
		end : float
			Position of the plasma end position
		"""
		return self._thisptr.end

	@end.setter
	def end(self,value):
		self._thisptr.end = value

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
		Number of iterations between particle sort, defaults to 16.
	
	See also
	--------
	em2d.Density
	"""

	cdef t_species _this
	cdef t_species* _thisptr
	cdef Density _density
	cdef str _name

	# Diagnostic types
	_diag_types  = { 'charge':CHARGE, 'pha':PHA, 'particles':PARTICLES }
	_pha_quants = { 'x1':X1, 'x2':X2, 'u1':U1, 'u2':U2, 'u3':U3 }

	def __cinit__( self, str name, const float m_q, list ppc = [1,1], *,
				  list ufl = [0.,0.,0.], list uth = [0.,0.,0.], Density density = None, n_sort = 16 ):

		self._thisptr = &self._this
		self._name = name
		self._this.m_q = m_q
		self._this.ppc = np.array(ppc, dtype=np.int32)
		self._this.ufl = np.array(ufl, dtype=np.float32)
		self._this.uth = np.array(uth, dtype=np.float32)
		self._this.n_sort = n_sort

		if ( density ):
			self._density = density.copy()
		else:
			# Use default uniform density
			self._density = Density()

	cdef new( self, t_species* ptr, int nx[], float box[], float dt ):
		self._thisptr = ptr
		n_sort = self._this.n_sort
		spec_new( self._thisptr, self._name.encode(), self._this.m_q, self._this.ppc,
			self._this.ufl, self._this.uth,
			nx, box, dt, self._density._thisptr )
		self._this.n_sort = n_sort

	def add( self, int[:] ix, float[:] x, float[:] u):
		"""add(ix, x, u)

		Adds a new particle to the particle buffer

		Parameters
		----------
		ix : `int[:]`
			New particle cell index
		x : `float[:]`
			New particle position inside the cell
		u : `float[:]`
			New particle (generalized velocity)
		"""
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
			Each quantity must be one of 'x1', 'x2', 'u1', 'u2' or 'u3'
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
		n : numpy.ndarray, (nx,ny)
			Charge density of particle species. Array will jhave the same
			shape as the simulation grid
		"""
		charge = np.zeros( shape = (self._thisptr.nx[1]+1,self._thisptr.nx[0]+1),
						   dtype = np.float32 )
		cdef float [:,:] buf = charge
		spec_deposit_charge( self._thisptr, &buf[0,0] )

		# Throw away guard cells
		return charge[ 0 : self._thisptr.nx[1], 0 : self._thisptr.nx[0] ]

	def phasespace( self, list quants, pha_nx, pha_range ):
		"""phasespace(quants, pha_nx, pha_range)
		
		Calculate phasespace density of particle species

		Parameters
		----------
		quants : list of str
			2 element list of quantities to use for "pha" diagnostics.
			Each quantity must be one of 'x1', 'x2', 'u1', 'u2' or 'u3'
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
		"""Cell sizes used for the species

		Returns
		-------
		dx : numpy.ndarray, float(2)
			Array of cell sizes
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
		ppc : numpy.ndarray, (2)
			Number of particles per cell along x and y
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

	@property
	def n_move(self):
		"""Number of cell moves by the moving simulation window

		Returns
		-------
		n_move : int
			Number of cell moves by the moving simulation window
		"""
		return self._thisptr.n_move

cdef float3 custom_ext_E( int ix, float dx, int iy, float dy, void *f ):
	"""Internal auxiliary function for using external E fields"""
	cdef ExternalField ext = <object> f
	val = ext.custom_func_E(ix,dx,iy,dy)
	cdef float3 e
	e.x = val[0]
	e.y = val[1]
	e.z = val[2]
	return e

cdef float3 custom_ext_B( int ix, float dx, int iy, float dy, void *f ):
	"""Internal auxiliary function for using external B fields"""
	cdef ExternalField ext = <object> f
	val = ext.custom_func_B(ix,dx,iy,dy)
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
		self._thisptr.E_type = self._ext_types[value]

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
		self._thisptr.B_type = self._ext_types[value]

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


cdef float3 custom_init_E( int ix, float dx, int iy, float dy, void *f ):
	"""Internal auxiliary function for using initial E fields"""
	cdef InitialField init = <object> f
	val = init.custom_func_E(ix,dx,iy,dy)
	cdef float3 e
	e.x = val[0]
	e.y = val[1]
	e.z = val[2]
	return e

cdef float3 custom_init_B( int ix, float dx, int iy, float dy, void *f ):
	"""Internal auxiliary function for using initial E fields"""
	cdef InitialField init = <object> f
	val = init.custom_func_B(ix,dx,iy,dy)
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
		"""Python function used for 'custom' initial E field type

		Returns
		-------
		func : function
			Python function for calculating all 3 components of initial
			electric field at every cell
		"""
		return self.custom_func_E

	@E_custom.setter
	def E_custom(self,value):
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
	creating a `em2d.Simulation` object.

	See Also
	--------
	em2d.Simulation
	"""

	cdef t_emf* _thisptr

	# Diagnostic types
	_diag_types = { 'E' : EFLD,	'B' : BFLD, 'Ep' : EPART, 'Bp': BPART }

	cdef associate( self, t_emf* ptr ):
		self._thisptr = ptr

	def report( self, str type, char fc ):
		"""report( type, fc )

		Save diagnostic information to disk. Files will be saved in the
		EMF directory below the current working directory.

		Parameters
		----------
		type : {'E','B','Ep','Bp'}
			Type of data to save, must be one of 'E' (electric field),
			'B' (magnetic field), 'Ep' (self-consistent + external
			electric field), 'Bp' (self-consistent + external magnetic
			field)
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
	def nx(self):
		"""Grid size used for the EMF object

		Returns
		-------
		nx : numpy.ndarray, int(2)
			Number of grid cells [nx,ny] for the simulation
		"""
		return self._thisptr.nx

	@property
	def dx(self):
		"""Cell sizes used for the EMF object

		Returns
		-------
		dx : numpy.ndarray, float(2)
			Array of cell sizes
		"""
		return self._thisptr.dx

	@property
	def box(self):
		"""Simulation box physical size

		Returns
		-------
		box : numpy.ndarray, float(2)
			Simulation box size
		"""
		return self._thisptr.box

	@property
	def n_move(self):
		"""Number of cell moves by the moving simulation window

		Returns
		-------
		n_move : int
			Number of cell moves by the moving simulation window
		"""
		return self._thisptr.n_move

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
	def solver_type(self):
		"""Field solver algorithm used

		Returns
		-------
		solver : 'Yee'
			Field solver, currently only 'Yee' is supported
		"""
		return "Yee"

	@property
	def Ex( self ):
		"""Ex field component

		Grid of (scalar) Ex field values excluding guard cells

		Returns
		-------
		Ex : numpy.ndarray, [nx,ny]
			X component values of the electric field.
		"""
		cdef float *buf = <float *> self._thisptr.E_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 0]

	@property
	def Ey( self ):
		"""Ey field component

		Grid of (scalar) Ey field values excluding guard cells

		Returns
		-------
		Ey : numpy.ndarray, [nx,ny]
			Y component values of the electric field.
		"""
		cdef float *buf = <float *> self._thisptr.E_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 1]

	@property
	def Ez( self ):
		"""Ez field component

		Grid of (scalar) Ey field values excluding guard cells

		Returns
		-------
		Ez : numpy.ndarray, [nx,ny]
			Z component values of the electric field.
		"""
		cdef float *buf = <float *> self._thisptr.E_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 2]

	@property
	def Bx( self ):
		"""Bx field component

		Grid of (scalar) Bx field values excluding guard cells

		Returns
		-------
		Bx : numpy.ndarray, [nx,ny]
			X component values of the magnetic field.
		"""
		cdef float *buf = <float *> self._thisptr.B_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 0]

	@property
	def By( self ):
		"""By field component

		Grid of (scalar) By field values excluding guard cells

		Returns
		-------
		By : numpy.ndarray, [nx,ny]
			Y component values of the magnetic field.
		"""
		cdef float *buf = <float *> self._thisptr.B_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 1]

	@property
	def Bz( self ):
		"""Bz field component

		Grid of (scalar) Bz field values excluding guard cells

		Returns
		-------
		Bz : numpy.ndarray, [nx,ny]
			Z component values of the magnetic field.
		"""
		cdef float *buf = <float *> self._thisptr.B_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 2]

	@property
	def Ex_part( self ):
		"""Ex field component (as seen by particles)

		Grid of (scalar) Ex field values excluding guard cells. If
		external fields are enabled this will be the sum of the self-
		consistent fields and the external fields

		Returns
		-------
		Ex_part : numpy.ndarray, [nx,ny]
			X component values of the total electric field.
		"""
		cdef float *buf
		if (self._thisptr.ext_fld.E_part_buf == NULL):
			buf = <float *> self._thisptr.E_buf
		else:
			buf = <float *> self._thisptr.ext_fld.E_part_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 0]

	@property
	def Ey_part( self ):
		"""Ey field component (as seen by particles)

		Grid of (scalar) Ey field values excluding guard cells. If
		external fields are enabled this will be the sum of the self-
		consistent fields and the external fields

		Returns
		-------
		Ey_part : numpy.ndarray, [nx,ny]
			Y component values of the total electric field.
		"""
		cdef float *buf
		if (self._thisptr.ext_fld.E_part_buf == NULL):
			buf = <float *> self._thisptr.E_buf
		else:
			buf = <float *> self._thisptr.ext_fld.E_part_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 1]

	@property
	def Ez_part( self ):
		"""Ez field component (as seen by particles)

		Grid of (scalar) Ex field values excluding guard cells. If
		external fields are enabled this will be the sum of the self-
		consistent fields and the external fields

		Returns
		-------
		Ez_part : numpy.ndarray, [nx,ny]
			Z component values of the total electric field.
		"""
		cdef float *buf
		if (self._thisptr.ext_fld.E_part_buf == NULL):
			buf = <float *> self._thisptr.E_buf
		else:
			buf = <float *> self._thisptr.ext_fld.E_part_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 2]

	@property
	def Bx_part( self ):
		"""Bx field component (as seen by particles)

		Grid of (scalar) Bx field values excluding guard cells. If
		external fields are enabled this will be the sum of the self-
		consistent fields and the external fields

		Returns
		-------
		Bx_part : numpy.ndarray, [nx,ny]
			X component values of the total magnetic field.
		"""
		cdef float *buf
		if (self._thisptr.ext_fld.B_part_buf == NULL):
			buf = <float *> self._thisptr.B_buf
		else:
			buf = <float *> self._thisptr.ext_fld.B_part_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 0]

	@property
	def By_part( self ):
		"""By field component (as seen by particles)

		Grid of (scalar) Bx field values excluding guard cells. If
		external fields are enabled this will be the sum of the self-
		consistent fields and the external fields

		Returns
		-------
		By_part : numpy.ndarray, [nx,ny]
			Y component values of the total magnetic field.
		"""
		cdef float *buf
		if (self._thisptr.ext_fld.B_part_buf == NULL):
			buf = <float *> self._thisptr.B_buf
		else:
			buf = <float *> self._thisptr.ext_fld.B_part_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 1]

	@property
	def Bz_part( self ):
		"""Bz field component (as seen by particles)

		Grid of (scalar) Bx field values excluding guard cells. If
		external fields are enabled this will be the sum of the self-
		consistent fields and the external fields

		Returns
		-------
		Bz_part : numpy.ndarray, [nx,ny]
			Z component values of the total magnetic field.
		"""
		cdef float *buf
		if (self._thisptr.ext_fld.B_part_buf == NULL):
			buf = <float *> self._thisptr.B_buf
		else:
			buf = <float *> self._thisptr.ext_fld.B_part_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 2]


cdef class Laser:
	"""Laser(type='plane', start=0.0, fwhm=0.0, rise=0.0, flat=0.0, fall=0.0, a0=0.0, omega0=0.0, polarization=0.0, W0=0.0, focus=0.0, axis=0.0)
	
	Class representing laser pulses. Laser pulses are added to the
	simulation using the `Simulation.add_laser()` method.

	Parameters
	----------
	type : {'plane', 'gaussian'}, optional
		Type of laser pulse to launch, defaults to 'plane'
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
	W0 : float, optional
		Gaussian pulse waist size in simulation units, defaults to 0
	focus : float, optional
		Position of focal plane for Gaussian pulse in simulation units,
		may be set outside of the simulation box, defaults to 0
	axis : float, optional
		y position of the propagation axis for Gaussian pulses in
		simulation units, defaults to 0.

	See Also
	--------
	em2d.Simulation
	"""

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
		"""FWHM of laser envelope
		
		Returns
		-------
		fwhm : float
			Full width at half-max of the laser pulse
		"""
		return self._thisptr.fwhm

	@fwhm.setter
	def fwhm(self,value):
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

	@property
	def W0(self):
		"""Gaussian pulse waist size in simulation units
		
		Returns
		-------
		W0 : float
			Gaussian envelope waist size
		"""
		return self._thisptr.W0

	@W0.setter
	def W0(self,value):
		self._thisptr.W0 = value

	@property
	def focus(self):
		"""Position of focal plane for Gaussian pulse in simulation units
		
		Returns
		-------
		focus : float
			Focal plane position
		"""
		return self._thisptr.focus

	@focus.setter
	def focus(self,value):
		self._thisptr.focus = value

	@property
	def axis(self):
		"""Y position of the propagation axis for Gaussian pulses
		
		Returns
		-------
		focus : float
			Y axis position
		"""
		return self._thisptr.axis

	@axis.setter
	def axis(self,value):
		self._thisptr.axis = value


cdef class Current:
	"""Current()
	
	Electric current density class

	This class allows access to the electric current density data
	structures in the simulation. An object of this class is created
	automatically when creating a `em2d.Simulation` object.

	See Also
	--------
	em2d.Simulation
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
		Jx : numpy.ndarray, [nx,ny]
			X component values of the current density
		"""
		cdef float *buf = <float *> self._thisptr.J_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 0]

	@property
	def Jy( self ):
		"""Jy current density component

		Grid of (scalar) Jy field values excluding guard cells

		Returns
		-------
		Jy : numpy.ndarray, [nx,ny]
			Y component values of the current density
		"""
		cdef float *buf = <float *> self._thisptr.J_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 1]

	@property
	def Jz( self ):
		"""Jz current density component

		Grid of (scalar) Jz field values excluding guard cells

		Returns
		-------
		Jz : numpy.ndarray, [nx,ny]
			Z component values of the current density
		"""
		cdef float *buf = <float *> self._thisptr.J_buf
		cdef int nx = self._thisptr.gc[0][0] + self._thisptr.nx[0] + self._thisptr.gc[0][1]
		cdef int ny = self._thisptr.gc[1][0] + self._thisptr.nx[1] + self._thisptr.gc[1][1]
		tmp = np.asarray( <float [:ny,:nx,:3]> buf )
		return tmp[ self._thisptr.gc[1][0] : self._thisptr.gc[1][0] + self._thisptr.nx[1], \
					self._thisptr.gc[0][0] : self._thisptr.gc[1][0] + self._thisptr.nx[0], 2]


cdef class Smooth:
	"""Smooth(xtype='none', ytype = 'none', xlevel=0, ylevel=0)

	Digital smoothing parameter for the Electric current density.

	Parameters
	----------
	xtype : {'none', 'binomial', 'compensated'}, optional
		Type of digital filtering to apply along the x direction, defaults
		to 'none'
	xlevel : int
		Number of filtering passes to apply along the x direction,
		defaults to 0
	ytype : {'none', 'binomial', 'compensated'}, optional
		Type of digital filtering to apply along the y direction, defaults
		to 'none'
	ylevel : int
		Number of filtering passes to apply along the y direction,
		defaults to 0
	"""

	cdef t_smooth* _thisptr

	# Filter types
	_filter_types = {'none'        : NONE,
					 'binomial'    : BINOMIAL,
					 'compensated' : COMPENSATED}

	def __cinit__( self, *, str xtype = 'none', str ytype = 'none',
				   int xlevel = 0, int ylevel = 0 ):

		self._thisptr = <t_smooth *> calloc(1, sizeof(t_smooth))

		self._thisptr.xtype = self._filter_types[ xtype ]
		self._thisptr.ytype = self._filter_types[ ytype ]

		self._thisptr.xlevel = xlevel
		self._thisptr.ylevel = ylevel

	def __dealloc__(self):
		free( self._thisptr )

	@property
	def xtype(self):
		"""Type of digital filtering to use along x

		Returns
		-------
		xtype : {'none', 'binomial', 'compensated'}
			Type of digital filtering to use along x
		"""
		tmp = {NONE:'none', BINOMIAL:'binomial', COMPENSATED:'compensated'}
		return tmp[self._thisptr.xtype]

	@xtype.setter
	def xtype( self, value ):
		self._thisptr.xtype = self._filter_types[value]

	@property
	def xlevel(self):
		"""Level of digital filtering to use along x

		Returns
		-------
		xlevel : int
			Number of filtering passes to apply along x
		"""
		return self._thisptr.xlevel

	@xlevel.setter
	def xlevel( self, value ):
		self._thisptr.xlevel = value

	@property
	def ytype(self):
		"""Type of digital filtering to use along y

		Returns
		-------
		ytype : {'none', 'binomial', 'compensated'}
			Type of digital filtering to use along y
		"""
		tmp = {NONE:'none', BINOMIAL:'binomial', COMPENSATED:'compensated'}
		return tmp[self._thisptr.ytype]

	@ytype.setter
	def ytype( self, value ):
		self._thisptr.ytype = self._filter_types[value]

	@property
	def ylevel(self):
		"""Level of digital filtering to use along y

		Returns
		-------
		ylevel : int
			Number of filtering passes to apply along y
		"""
		return self._thisptr.xlevel

	@ylevel.setter
	def ylevel( self, value ):
		self._thisptr.ylevel = value

cdef class Simulation:
	"""Simulation(nx, box, dt, species=None, report=None, mov_window=False, smooth=None, init_fld=None, ext_fld=None)
	
	ZPIC EM2D Simulation class

	Parameters
	----------
	nx : list of int
		Number of grid cells [nx,ny] for the simulation
	box : list of float
		Simulation box (phyiscal) size, in simulation units
	dt : float
		Simulation time step, in simulation units
	species : `Species` or list of `Species`, optional
		Particle species to use in the simulation, defaults to None
		(no particles)
	report : function, optional
		Python function used for simulation reporting, defaults to None
	mov_window : bool, optional
		Use a moving simulation window, defaults to False
	smooth : `Smooth`, optional
		Electric current smoothing options, defaults to None
	init_fld : `InitialField`, optional
		Initial EM fields for the simulation, defaults to None (0 initial
		fields)
	ext_fld : `ExternalField`, optional
		External EM fields for the simulation, defaults to None

	See also
	--------
	em2d.Species
	em2d.Smooth
	em2d.InitialField
	em2d.ExternalField
	"""

	cdef t_simulation *_thisptr
	cdef int n
	cdef float t

	cdef EMF emf
	cdef Current current
	cdef list species

	cdef object report

	def __cinit__( self, list nx, list box, float dt, *, species = None,
				   report = None, bint mov_window = False, Smooth smooth = None,
				   InitialField init_fld = None, ExternalField ext_fld = None ):

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

		# Set moving window
		if ( mov_window ):
			sim_set_moving_window( self._thisptr )
		
		# Set current smoothing
		if ( smooth ):
			sim_set_smooth( self._thisptr, smooth._thisptr )
		
		# Set initial fields
		if ( init_fld ):
			emf_init_fld( self.emf._thisptr, init_fld._thisptr )
		
		# Set external fields
		if ( ext_fld ):
			emf_set_ext_fld( self.emf._thisptr, ext_fld._thisptr )

	def __dealloc__(self):
		sim_delete( self._thisptr )
		free(self._thisptr)

	def set_moving_window(self):
		"""set_moving_window()
		
		Turns on moving window for the simulation

		Note
		----
		Use of this function has been deprecated and will be removed soon.
		Use the `mov_window` parameter of the `Simulation` class instead.
		"""
		sim_set_moving_window( self._thisptr )

	def set_smooth(self, Smooth smooth):
		"""set_smooth(smooth)
		
		Sets electric current smoothing (digital filtering) parameters

		Parameters
		----------
		smooth : `Smooth`
			`Smooth` object specifying smoothing parameters
		
		See Also
		--------
		em2d.Smooth

		Note
		----
		Use of this function has been deprecated and will be removed soon.
		Use the `smooth` parameter of the `Simulation` class instead.
		"""
		sim_set_smooth( self._thisptr, smooth._thisptr )

	def add_laser(self, Laser laser):
		"""add_laser(laser)
		
		Adds laser pulse to the simulation

		Parameters
		----------
		laser : `Laser`
			Laser pulse to be added to the simulation
		
		See Also
		--------
		em2d.Laser
		"""
		sim_add_laser( self._thisptr, laser._thisptr )

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
		species : list of `em2d.Species`
			Simulation particle species list
		"""
		return self.species
	
	@property
	def emf(self):
		"""Simulation EM fields object

		Returns
		-------
		emf : `em2d.EMF`
			Simulation EM fields object
		"""
		return self.emf

	@property
	def current(self):
		"""Simulation electric current density object

		Returns
		-------
		current : `em2d.Current`
			Simulation electric current density object
		"""
		return self.current

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
		"""Cell sizes used for the simulation

		Returns
		-------
		dx : numpy.ndarray, float(2)
			Array of cell sizes
		"""
		return self.emf.dx

	@property
	def nx(self):
		"""Grid size used for the simulation

		Returns
		-------
		nx : numpy.ndarray, int(2)
			Number of grid cells [nx,ny] for the simulation
		"""
		return self.emf.nx

	@property
	def box(self):
		"""Simulation box physical size

		Returns
		-------
		box : numpy.ndarray, float(2)
			Simulation box size
		"""
		return self.emf.box

	@property
	def n_move(self):
		"""Number of cell moves by the moving simulation window

		Returns
		-------
		n_move : int
			Number of cell moves by the moving simulation window
		"""
		return self.emf.n_move

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
