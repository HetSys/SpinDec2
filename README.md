
```
                %%%%%%\          %%\        %%%%%%\                    %%%%%%\
               %%  __%%\         \__|       %%   %%\                        %%\
               %% /  \__|%%%%%%\ %%\%%%%%%\ %%    %%\ %%%%%%\  %%%%%%%\     %% |
               \%%%%%%\ %%  __%%\%% %%  _%%\%%    %% %%    %%\%%  _____%%%%%%  |
                \____%%\%% /  %% %% %% / %% %%    %% %%%%%%%% %% /     %%  ___/
               %%\   %% %% |  %% %% %% | %% %%   %% /%%   ____%% |     %% |
               \%%%%%%  %%%%%%%  %% %% | %% %%%%%% / \%%%%%%%\\%%%%%%%\%%%%%%%\
                \______/%%  ____/\__\__| \__\_____/   \_______|\_______\_______|
                        %% |
                        %% |
                        \__|
```

# SpinDec2
This project models spinoidal decomposition of a binary alloy system and has the following features:
* Explicit Eulerian approach with either:
  * A constant diffusive mobility field
  * A non-constant diffusive mobility field using the Darkens equation
  * Temperature dependent atomic mobilities
 
* Pseudo-Spectral solver with constant mobility

For full documentation, please refer to the wiki: https://hetsys.github.io/SpinDec2/

The following is an example end state using the explicit Eulerian solver with an initial concentration range of 0.31 to 0.32:

![constant_end_3132](https://user-images.githubusercontent.com/78127892/170461306-6012d905-60c8-4752-a883-a6a565a6f858.png)

## Dependencies 
Prior to installation, please ensure that the following dependencies are installed: 

* bash
  * util-linux
  * getopt

* python
  * netCDF4
  * matplotlib
  * numpy
  * scipy
  * argparse

* fortran
  * netCDF
  * fftw

This has been tested on HPC systems, and is known to work with the following dependency versions:
* GCC/7.3.0-2.30  OpenMPI/3.1.1 netCDF-Fortran/4.4.4 FFTW/3.3.8
* GCC/10.3.0 OpenMPI/4.1.1 netCDF-Fortran/4.5.3 FFTW/3.3.9

## Download
To download the repository, clone the repository somewhere in your filesystem:

`git clone https://github.com/HetSys/SpinDec2`

## Installation 
An installation script is provided with the repository. All options can be viewed by running the script without any flags, or with `-h/--help`. 

All flags have a short and a long version. If using the short option, and the flag supports an additional argument, the argument is called by a single letter immediately after the flag without leaving whitespace. For example, for compilation with debugging: 

`./install.sh -cd`.

When using the long version, additional arguments can be specified like so (also using the compile flag):

`./install.sh --compile=debug`.

A description of the different flags and their arguments are given below:

### Compile 
To compile the code, run 

`./install.sh -c` or `./install.sh –-compile`.

The compilation argument takes an optional argument `d/debug`, however this is only for developer use. Running with this argument can cause decreases in performance.

After compilation, you will be prompted to add the binary to your default shell `$PATH` in `$HOME/.bashrc`. This feature only works for bash, so if you are using a different shell, you will have to do this manually. 

### Clean 
Using the option will remove all installed binaries, netcdf, and checkpoint files in `bin`. To do so, run 

`./install.sh -C` or `./install.sh –-clean`

The clean argument takes an option argument `c/confirm`. This will prompt for confirmation before removing binaries.

## Testing
### Running 
Unit tests can be automatically run from the install script. To use this feature, specify the argument 

`./install.sh -t` or `./install.sh --test`. 

This takes one of several required arguments:
- `c/compile` to only compile the tests 
- `r/run` to only run the tests 
- `b/both` to both compile and run the tests.

All test binaries are compiled to `test/test_bin/`, and the main executable is saved as `test_spindec`.

### Cleaning 
To clean all unit tests binaries, use 

`./install.sh -T` or `./install.sh --test-clean`.

## Contributors 
Geraldine Anis, Ben Gosling, Dylan Morgan, Matyas Parrag, and Anas Siddiqui

HetSys CDT, University of Warwick

Copyright © 2022, Geraldine Anis, Ben Gosling, Dylan Morgan, Matyas Parrag, and Anas Siddiqui
