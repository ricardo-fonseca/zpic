from distutils.core import setup, Extension
from Cython.Build import cythonize

from distutils import sysconfig

# Removes the debug options and replaces the optimization flag
cflags = sysconfig.get_config_var('CFLAGS')
cflags = cflags.replace(' -g', '')
cflags = cflags.replace(' -DNDEBUG', '')
cflags = cflags.replace('-O3', '-Ofast')
sysconfig._config_vars['CFLAGS'] = cflags

ext = Extension("em1d",
                sources=["em1d.pyx",
                "em1d/current.c",
				"em1d/emf.c",
				"em1d/particles.c",
				"em1d/random.c",
				"em1d/simulation.c",
				"em1d/timer.c",
				"em1d/zdf.c"]
)

setup(name="zpic_em1d",
      ext_modules = cythonize([ext]))
