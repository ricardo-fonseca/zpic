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

es1d = Extension("es1d",
                sources=["es1d.pyx",
				"../es1d/charge.c",
				"../es1d/fft.c",
				"../es1d/field.c",
				"../es1d/grid.c",
				"../es1d/main.c",
				"../es1d/particles.c",
				"../es1d/random.c",
				"../es1d/simulation.c",
				"../es1d/timer.c",
				"../es1d/zdf.c"]
)


setup(name="zpic",
      ext_modules = cythonize([em1d, es1d]))
