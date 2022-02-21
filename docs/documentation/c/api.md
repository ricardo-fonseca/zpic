---
title: C API reference
permalink: /documentation/c/api

layout: single

sidebar:
  nav: "docs"
---

The C source code is extensively documented using _docstrings_, users may find the documentation by opening either header (`*.h`) or source (`*.c`) files. Alternatively, an html version of the documentation may be generated using [Doxygen](https://www.doxygen.nl/). Check their website for instructions on intalling Doxygen on your system.

To generate the documentation just navigate into the appropriate source folder (e.g. `em2d`) and run `doxygen`:

```bash
$ cd em2d
$ doxygen
```

The documentation will be generated in the `docs` subfolder.
