#EM1D - Parameter Scan

This is a modification of the EM1D code that implements a parameter scan. The code will run a simulation multiple times, varying the required parameter(s), and storing the results on different directories for later processing.

The input file has two minor modifications when compared to the default files:

1. It must define the constant `n_scan`. This sets the number of parameter sets to try, e.g.

    ```C
    static const int n_scan = 5;
    ```
2. The `sim_init` routine includes an additional integer parameter, `opt_idx` that should be used to select the parameters for the simulation, e.g.:
    
    ```C
    void sim_init( t_simulation* sim, int opt_idx ){
      ...
      laser.a0 = 1.0 + 0.1 * opt_idx;
      ...
    }
    ```

Please see the `input/lwfa.c` file for a complete example.

## Output

For each parameter set _i_ (ranging from 0 to _n\_scan_ - 1), the code will create a directory named "output._i_" and write all diagnostic output for that parameter set in that directory.

## Using MPI (optional)

The code can use MPI to launch multiple simulations in parallel and speed up the parameter scan. To compile the code for MPI use just set the following options on the Makefile (assuming you have MPI wrapper compilers installed:)

```Make
CC = mpicc
CFLAGS = -D_MPI_ -Ofast -std=c99 -pedantic
```

Notice the `-D_MPI_`, this enables the MPI sections of the code.

When launching the code, the number of processes need not match the number of parameters. The algorithm will attempt to divide the parameter sets as evenly as possible across available processes. 