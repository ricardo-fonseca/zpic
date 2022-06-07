---
title: Getting started in C
permalink: /documentation/c/start

layout: single
toc: true

sidebar:
  nav: "docs"
---

## Choosing the code version

To run a ZPIC simulation in C users must first choose which code version is to be used and compile the code in the corresponding directory. Available versions are:

* `em1d`, `em1ds` - 1D electromagnetic code using a finite difference field solver (`em1d`) or a spectral field solver (`em1ds`)
* `em2d`, `em2ds` - 2D electromagnetic code using a finite difference field solver (`em2d`) or a spectral field solver (`em2ds`)
* `es1d` - 1D electrostatic code using a spectral field solver

Each top level directory (`em1d`, `em2d`, etc.) has a self-contained version of the code, with no additional dependencies. Just navigate to the required directory and run `make`. Detailed instruction on compiling the code(s) can be found in the [Compiling](/start/compile) page.

## Initializing the simulation

Simulation parameters are defined using a C file that is included in the `main.c` file of each code. For example, to launch the example Two-Stream instability simulation on the 1D electrostatic code (file `es1d/input/twostream.c`), you need to edit the `es1d/main.c` file and add:

```c
// Include Simulation parameters here
#include "input/twostream.c"
```

The file specifying the simulation parameters must define two functions:

* `sim_init(t_simulation* sim)`   - Initialization of the simulation
* `sim_report(t_simulation* sim)` - Simulation diagnostics

Please check the examples available on the `input` directories for each code.

## Main simulation parameters

A ZPIC simulation requires defining a set of numerical parameters regarding the simulation spatial / temporal domain, the simulation particles. Additionally, the user may optionally specify additional parameters, such as laser pulses.

These parameters are defined in the `sim_init()` function in the input file. This function *must* call the `sim_new()` function to initialize the `sim` (simulation) object:

```c
void sim_new( t_simulation* sim, int nx[], float box[], 
        float dt, float tmax, int ndump,
        t_species* species, int n_species );
```

The parameters are:

* `sim` -  Pointer to simulation object
* `nx` - Grid size for the simulation (scalar `int` in 1D codes)
* `box` - Physical size of the simulation box (scalar `float` in 1D codes) in simulation units
* `dt` - Simulation time step in simulation units
* `species` - Array of `t_species` objects specifying particle species. May be set to `NULL` for no particle species. See the [Particle injection](particles) page for details.
* `n_species` - Number of objects in the `species` array.

The following example sets the minimal parameters for a 2D ZPIC simulation (with no particles):

```c
void sim_init( t_simulation* sim ){

  // Simulation box
  int   nx[2]  = { 100, 100 };
  float box[2] = { 10.0, 10.0 };

  // Time step
  float dt = 0.07;
  float tmax = 70.0;

  // Initialize Simulation data
  sim_new( sim, nx, box, dt, tmax, 0, NULL, 0 );
}
```

Please check the `input` directories under each code for several complete examples.

### Simulation grid size and physical box

The input parameters must define the number of grid points in each direction of the simulation grid, and the physical dimensions (in simulation units) of this box. These dimensions are measured from the lower boundary of the first cell to the upper boundary of the last cell.

```c
  // Simulation box
  int   nx[2]  = { 100, 100 };
  float box[2] = { 10.0, 10.0 };

  // Initialize Simulation data
  sim_new( sim, nx, box, ... );
```

### Time step and total simulation time

The input parameters must define the simulation time step and the total simulation time (in simulation units)

```c
  // Time step
  float dt = 0.07;
  float tmax = 70.0;

  // Initialize Simulation data
  sim_new( sim, nx, box, dt, tmax, ... );
```

The code will check if the time step is compatible with the specified resolution.

### Diagnostic frequency

The input parameters must also define the frequency at which the `sim_report()` function is called from the main loop. For this purpose the user specifies the number of iterations between these function calls:

```c
  // Diagnostic frequency - call sim_report() every 50 iterations
  int ndump = 50;

  // Initialize Simulation data
  sim_new( sim, nx, box, dt, tmax, ndump, ... );
```

Setting this parameter to 0 disables all diagnostics. For details on simulation diagnostics and the `sim_report()` function see the [Simulation Diagnostics](diag) page.

### Particles

The input parameters may define an arbitrary number of particle species, defining the required plasma properties. See the [Particle injection](particles) page for details.

## Running the code

Once you have finished compiling your code (including the simulation parameters as described above), you will have a `zpic` binary in the source code folder. Executing this binary will run your simulation. Since output files will be written in the same directory where you run your code it is recommended that you do not run the code in the same directory where the binary was created.

The following example compiles the `em1d` code () and runs it in a directory named `test`:

```shell
$ cd source/zpic/em1d
$ make
gcc -c -Ofast -std=c99 -pedantic current.c -o current.o
(...)
gcc -Ofast -std=c99 -pedantic current.o emf.o particles.o random.o timer.o main.o simulation.o zdf.o -lm -o zpic
$ mkdir test
$ cd test
$ ../zpic
Starting simulation ...

n = 0, t = 0.000000
(...)
n = 500, t = 50.000000

Simulation ended.
$
```
