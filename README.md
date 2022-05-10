
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
Modelling the phase field of spinoidal decomposition

[LIST FEATURES HERE]

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

A description of the different flags and their arugments are given below:

### Compile 
To compile the code, run 

`./install.sh -c` or `./install.sh –-compile`.

The compilation argument takes an optional argument `d/debug`, however this is only for developer use. Running with this argument can cause decreases in performance.

After compilation, you will be prompted to add the binary to your default shell `$PATH` in `$HOME/.bashrc`. This feature only works for bash, so if you are using a different shell, you will have to do this manually. 

### Clean 
Using the option will remove all installed binaries, netcdf, and checkpoint files in `bin`. To do so, run 

`./install.sh -C/–-clean`

The clean argument takes an option argument `c/confirm`. This will prompt for confirmation before removing binaries.

### Test 
TODO

### Example 
TODO
