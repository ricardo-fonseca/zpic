# The ZPIC educational code suite

Particle-in-Cell (PIC) codes are used in almost all areas of plasma physics, such as fusion energy research, plasma accelerators, space physics, ion propulsion, and plasma processing, and many other areas. Leveraging on our expertise and experience from the development and use of the OSIRIS PIC code, we have developed a suite of 1D/2D fully relativistic electromagnetic PIC codes, as well as 1D electrostatic. These codes are self-contained and require only a standard laptop/desktop computer with a C compiler to be run. The output files are written in a new file format called ZDF that can be easily read using the supplied routines in a number of languages, such as Python, and IDL. The code suite also includes a number of example problems that can be used to illustrate several textbook and advanced plasma mechanisms, including instructions for parameter space exploration. We also invite contributions to this repository of test problems that will be made freely available to the community provided the input files comply with the format defined by the ZPIC team.

Please see the [doc](https://github.com/zambzamb/zpic/tree/master/doc) folder for documentation on using the code.

The directory structure is organized as follows:
* [**em1d**](https://github.com/zambzamb/zpic/tree/master/em1d) - 1D electromagnetic (finite difference)
* [**em1ds**](https://github.com/zambzamb/zpic/tree/master/em1ds) - 1D electromagnetic (spectral)
* [**em2d**](https://github.com/zambzamb/zpic/tree/master/em2d)  - 2D electromagnetic (finite difference)
* [**em2ds**](https://github.com/zambzamb/zpic/tree/master/em2ds) - 2D electromagnetic (spectral)
* [**es1d**](https://github.com/zambzamb/zpic/tree/master/es1d)  - 1D electrostatic
* [**mods**](https://github.com/zambzamb/zpic/tree/master/mods)  - Modified versions of the base codes
