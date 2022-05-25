User Documentation
==================

Getting Started
---------------

### Download
To get the source code, the repository can be cloned to your chosen directory by typing the following command into a terminal (making sure you are in the correct directory):
`$ git clone https://github.com/HetSys/SpinDec2`

### Installation
An installation script is provided with the repository. This can be used with several flags. To view the flags, run `./install.sh -h`, `./install.sh --help` or just `./install.sh` when inside the repository. A help screen will be displayed as shown below:

       %%%%%%\\          %%\\        %%%%%%\\                    %%%%%%\\
      %%  \_\_%%\\         \\\_\_|       %%   %%\\                        %%\\
      %% /  \\\_\_|%%%%%%\\ %%\\%%%%%%\\ %%    %%\\ %%%%%%\\  %%%%%%%\\     %% |
      \\%%%%%%\\ %%  \_\_%%\\%% %%  \_%%\\%%    %% %%    %%\\%%  \_\_\_\_\_%%%%%%  |
       \\\_\_\_\_%%\\%% /  %% %% %% / %% %%    %% %%%%%%%% %% /     %%  \_\_\_/
      %%\\   %% %% |  %% %% %% | %% %%   %% /%%   \_\_\_\_%% |     %% |
      \\%%%%%%  %%%%%%%  %% %% | %% %%%%%% / \\%%%%%%%\\\\%%%%%%%\\%%%%%%%\\
       \\\_\_\_\_\_\_/%%  \_\_\_\_/\\\_\_\\\_\_| \\\_\_\\\_\_\_\_\_/   \\\_\_\_\_\_\_\_|\\\_\_\_\_\_\_\_\\\_\_\_\_\_\_\_|
               %% |
               %% |
               \\\_\_|

      Modelling Spinoidal Decomposition Using a Phase Field Approach
      \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

     usage: spindec \[-h\]
                    \[-c ARGS\]
                    \[-C ARGS\]
                    \[-t ARGS\]
                    \[-T\]
     options:
       -h, --help              show this help message and exit

       -c, --compile ARGS      compile the code with optional debug or profile option
                               optional arguments: \[ none | d/debug \] (default=none)

       -C, --clean ARGS        remove compiled binaries from repository
                               optional arguments: \[ none | c/confirm \] (default=none)

       -t, --test ARGS         run automated unit tests
                               required arguments: \[ c/compile | r/run | b/both \]

       -T, --test-clean        clean test binaries

The code can be compiled with the installation script as follows:

`$ ./install.sh -c` or `$ ./install.sh --compile`

If the `$PATH` to the `spindec` binary isn't already in `~/.bashrc` , a prompt will appear asking if the user wishes to add the `/bin` location to `$PATH` . The shell will then need to be reloaded using `source ~/.bashrc` for `spindec` to appear as executable. This is the recommended method in which to use the binary, and the rest of the tutorial will assume this step has been taken.

### Using the clean flag to remove compiled binaries
To recompile the code, or just to remove binaries from a previous compilation, the clean flag can be used:

`$ ./install.sh -C` or `$ ./install.sh --clean`

To ask for confirmation before removing, the `c/confirm` argument can be specified:

`$ ./install.sh -Cc` or `$ ./install.sh --compile=confirm`

### Running the code
After compilation of the code, the binaries will be stored in the sub-directory `/bin` from the repository home. To execute the code, it is highly recommended to use the approach of creating a directory per calculation. This directory should contain an input file called `input.txt` of which the syntax is discussed further [here](input.html). Firstly the desired number of OpenMP threads needs to be specified. Then the computation can be run in the directory where the input file is located:

`$ export OMP_NUM_THREADS=n`
`$ mpirun -np n spindec`

where `n` is the number of processors to run with. If it is desired to only run with a single MPI thread, it is possible to call `spindec` simply like so (after specifiying OpenMP threads):

`$ spindec`

This uses the parameters provided in the input file to run the simulation and evolve the concentration grid over a the number of time steps requested.

After the computation has finished running, the `vis_spindec.py` script can be used to generate graphs and animations from the output of the computation. This should also have been added to `$PATH` at the same time as `spindec`.

`$ vis_spindec.py`

### Testing 

#### Running

Unit tests can be automatically run from the install script. To use this feature, specify the argument with one of several mandatory arugments

* `$ ./install.sh -tc` / `$ ./install.sh --test=compile`: compile the tests
* `$ ./install.sh -tr` / `$ ./install.sh --test=run`: run the tests
* `$ ./install.sh -tb` / `$ ./install.sh --test=both`: both compile and run the tests

#### Cleaning 

To clean all the test binaries:

`$ ./install.sh -T` or `$ ./install.sh --test-clean`

Examples
--------
Below are some simple examples of how to use SpinDec2 for the various cases of methods and Mobilities.
*   Euler - constant M - Mobility is taken to be constant everywhere at all times, [Serial Example 1](euler_example_1.html).
*   Euler - M = M(c) - Mobility is dependent on the order paramter c over time, [Serial Example 2](euler_example_2.html).
*   Euler - M = M(c,T) The atomic mobilities are defined in terms of local tempeartures, [Serial Example 3](euler_example_3.html).
*   Validation Example, [50/50 Split Example](test_example.html) .
*   Running in parallel using MPI, [Parallel Example](parra_example.html) .
*   Pseudo-Spectral Method, [Pseudo-Spectral Example](PS_example.html) .
*   How to use check-pointing, [Check-Pointing Example](PS_example.html) .
Authors: Anas Siddiqui, Ben Gosling, Dyaln Morgan, Geraldine Anis, Matyas Parrag
*   Generated by [![doxygen](doxygen.svg)](https://www.doxygen.org/index.html) 1.9.3
