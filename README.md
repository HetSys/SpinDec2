
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
Modelling Spinoidal Decomposition Using a Phase Field Approach

[LIST FEATURES HERE]

[SHOW SOME PRETTY PICTURES/DIAGRAMS/GRAPHS]

## Dependencies 
Prior to installation, please ensure that the following dependencies are installed: 
- util-linux
- [WRITE MORE HERE]

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
