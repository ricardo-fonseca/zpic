# To compile the extension type:
#	python setup.py build_ext --inplace
#
# To remove all compilation results
#	python setup.py clean
#	rm *.so
#

# Note that from a file named fname.pyx, Cython will create another file named
# fname.c without warning. Long story short, don't name your .pyx file with the
# same name as any other file in the same dir 

cdef extern from "ccode.h":
	ctypedef struct t_struct:
		int a
		float b

	int prints( t_struct *s )

# This can be called with a dictionary like this:
#
# import cwrap
# s = {'a':123, 'b':23.4}
# cwrap.print_s(s)

def print_s( t_struct s ):
	prints(&s)

