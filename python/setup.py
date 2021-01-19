from distutils.core import setup, Extension
from Cython.Build import cythonize

from distutils import sysconfig

# Removes the debug options and replaces the optimization flag
cflags = sysconfig.get_config_var('CFLAGS')
cflags = cflags.replace(' -g', '')
cflags = cflags.replace(' -DNDEBUG', '')
cflags = cflags.replace('-O3', '-Ofast')

# Sets the standard to c99
cflags += ' -std=c99'

sysconfig._config_vars['CFLAGS'] = cflags

# Set custom compiler
sysconfig._config_vars['CC'] = 'gcc'

# Finit difference Electromagnetic codes
em1d = Extension("em1d",
                sources=["em1d.pyx",
                "../em1d/current.c",
				"../em1d/emf.c",
				"../em1d/particles.c",
				"../em1d/random.c",
				"../em1d/simulation.c",
				"../em1d/timer.c",
				"../em1d/zdf.c"]
)

em2d = Extension("em2d",
                sources=["em2d.pyx",
                "../em2d/current.c",
				"../em2d/emf.c",
				"../em2d/particles.c",
				"../em2d/random.c",
				"../em2d/simulation.c",
				"../em2d/timer.c",
				"../em2d/zdf.c"]
)

# Electrostatic (spectral) codes
es1d = Extension("es1d",
                sources=["es1d.pyx",
				"../es1d/charge.c",
				"../es1d/fft.c",
				"../es1d/field.c",
				"../es1d/grid.c",
				"../es1d/particles.c",
				"../es1d/random.c",
				"../es1d/simulation.c",
				"../es1d/timer.c",
				"../es1d/zdf.c"]
)

# Spectral Electromagnetic codes
em1ds = Extension("em1ds",
                sources=["em1ds.pyx",
				"../em1ds/filter.c",
				"../em1ds/charge.c",
				"../em1ds/current.c",
				"../em1ds/emf.c",
				"../em1ds/fft.c",
				"../em1ds/grid.c",
				"../em1ds/particles.c",
				"../em1ds/random.c",
				"../em1ds/simulation.c",
				"../em1ds/timer.c",
				"../em1ds/zdf.c"]
)

em2ds = Extension("em2ds",
                sources=["em2ds.pyx",
				"../em2ds/charge.c",
				"../em2ds/current.c",
				"../em2ds/emf.c",
				"../em2ds/fft.c",
				"../em2ds/filter.c",
				"../em2ds/grid2d.c",
				"../em2ds/particles.c",
				"../em2ds/random.c",
				"../em2ds/simulation.c",
				"../em2ds/timer.c",
				"../em2ds/zdf.c"]
)


setup(name="zpic",
      ext_modules = cythonize([em1d, em2d, es1d, em1ds, em2ds]))
