
# ZPIC Python interface

## Introduction

This module allows you to use any of the ZPIC codes from a Python environment, including Jupiter notebooks. All code features that are available in the main C versions are available, and it is possible to use this module to extend the codes in terms of features/diagnostics using only Python.

## Dependencies

The python interface was developed using Cython (0.29.x). The notebooks are based on Jupyter use NumPy and Matplotlib for data and visualization.

## Compilation

To compile the module just navigate to this folder and run `make`. There will be a python module created for each of the ZPIC codes and placed in the `lib` directory. Be sure to add this library to your python path to use these modules:

```python
import sys
sys.path.append("path_to_module/lib")
```

## Getting Started

Once you compiled the code, the best way to start using the code is to open the [ZPIC](notebooks/tutorial/ZPIC.ipynb) notebook (in the `notebooks/tutorial` folder) in a Jupyter session. This notebook covers the basic steps for running ZPIC simulations in a notebook environment.

## Examples

This module includes a set of Jupyter notebooks exemplifying the use of ZPIC in multiple situations. These are found in the `notebooks` directory, and are organized as follows:

* `tutorial` - Notebooks introducing code usage and functionalities
* `classroom` - Examples of ZPIC notebooks that can be used in a classrom to showcase some of the most fundamental plasma physics phenomena.
* `papers` - Notebooks reproducing and extend the work done in seminal plasma physics papers

Check the [notebooks](notebooks) directory for more details.
