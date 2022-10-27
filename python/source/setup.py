from setuptools import setup, Extension
from Cython.Build import cythonize

# Removes the debug options and replaces the optimization flag
import sysconfig
from setuptools.command.build_ext import build_ext

cflags = sysconfig.get_config_var('CFLAGS')
cflags = cflags.replace(' -g', '')
cflags = cflags.replace('-O3', '-Ofast')
cflags += ' -std=c99 -fPIC'

# Removes unrecognized gcc option if present
cflags = cflags.replace(' -fexceptionsrecord-gcc-switches', '')

cc = sysconfig.get_config_var('CC')

class custom_build_ext(build_ext):
	def build_extensions(self):
		self.compiler.set_executable("compiler_so", cc + " " + cflags)
		build_ext.build_extensions(self)

# Get the path to the C source file root
import pathlib
csource = str(pathlib.Path("../../").resolve()) + "/"

# Finit difference Electromagnetic codes
em1d = Extension("em1d",
                sources=["em1d.pyx",
                csource + "em1d/current.c",
				csource + "em1d/emf.c",
				csource + "em1d/particles.c",
				csource + "em1d/random.c",
				csource + "em1d/simulation.c",
				csource + "em1d/timer.c",
				csource + "em1d/zdf.c"]
)

em2d = Extension("em2d",
                sources=["em2d.pyx",
                csource + "em2d/current.c",
				csource + "em2d/emf.c",
				csource + "em2d/particles.c",
				csource + "em2d/random.c",
				csource + "em2d/simulation.c",
				csource + "em2d/timer.c",
				csource + "em2d/zdf.c"]
)

# Electrostatic (spectral) codes
es1d = Extension("es1d",
                sources=["es1d.pyx",
				csource + "es1d/charge.c",
				csource + "es1d/fft.c",
				csource + "es1d/field.c",
				csource + "es1d/grid.c",
				csource + "es1d/particles.c",
				csource + "es1d/random.c",
				csource + "es1d/simulation.c",
				csource + "es1d/timer.c",
				csource + "es1d/zdf.c"]
)

# Spectral Electromagnetic codes
em1ds = Extension("em1ds",
                sources=["em1ds.pyx",
				csource + "em1ds/filter.c",
				csource + "em1ds/charge.c",
				csource + "em1ds/current.c",
				csource + "em1ds/emf.c",
				csource + "em1ds/fft.c",
				csource + "em1ds/grid.c",
				csource + "em1ds/particles.c",
				csource + "em1ds/random.c",
				csource + "em1ds/simulation.c",
				csource + "em1ds/timer.c",
				csource + "em1ds/zdf.c"]
)

em2ds = Extension("em2ds",
                sources=["em2ds.pyx",
				csource + "em2ds/charge.c",
				csource + "em2ds/current.c",
				csource + "em2ds/emf.c",
				csource + "em2ds/fft.c",
				csource + "em2ds/filter.c",
				csource + "em2ds/grid2d.c",
				csource + "em2ds/particles.c",
				csource + "em2ds/random.c",
				csource + "em2ds/simulation.c",
				csource + "em2ds/timer.c",
				csource + "em2ds/zdf.c"]
)

# Compile extensions
setup(name="zpic",
    ext_modules = cythonize([em1d, em2d, es1d, em1ds, em2ds]), 
	zip_safe=False,
    cmdclass={"build_ext": custom_build_ext}
)
