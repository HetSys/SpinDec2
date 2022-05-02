#!/bin/bash 
# install.sh
# Compilation/install script for SpinDec

# Exit if something fails
set -e

### Compilation ###
compile () {
    # Option for choosing the compiler
    if [[ $1 == "t" ]] || [[ $1 == "true" ]]; then
        comp_line="gfortran -g -std=f2008 -Wall -fimplicit-none -fcheck=all -Wextra -pedantic -fbacktrace"
    elif [[ -z $1 ]] || [[ $1 == "f" ]] || [[ $1 == "false" ]]; then
        comp_line="gfortran -g"
    else
        echo -e "Not a valid option for -c/--compile\n"
        help_message
        exit 2
    fi

    # Add program files from prog_files.txt
    input_file="./prog_files.txt"
    readarray -t prog_files < "$input_file"

    # Binary name and location
    compd_file="./bin/spindec"

    # NetCDF flags
    flibs=`nf-config --flibs`
    fflags=`nf-config --fflags`

    echo -e "Compile line: $comp_line\n"

    # Compile
    # $comp_line $fflags $prog_files $flibs -o $compd_file

    # Add binary to $PATH (with some checks)
    while true; do
        read -p 'Add binary to $PATH? [Y/n] ' confirmation

        if [[ $confirmation =~ ^[Yy]$ ]] || [[ $confirmation == '' ]]; then
            if grep -Fq "spindec" $HOME/.bashrc; then
                echo 'The binary is already in your $PATH'
                break
            else
                if [[ "$SHELL" == *"bash"* ]]; then
                    working_dir=`pwd`
                    echo "export PATH=`pwd`/bin/spindec:$PATH" >> $HOME/.bashrc
                    echo 'Binary added to $PATH'
                    break
                else
                    echo 'Unable to add to $PATH as bash is not your default shell'
                    break
                fi
            fi

        elif [[ $confirmation =~ ^[Nn]$ ]]; then
            echo 'Not added to $PATH'
            break
        else
            echo 'Not a valid option'
        fi
    done
}

auto_test () {
    # Tests directory
    tests_dir="./test/"
}

example () {
    example_dir="./test/"
}

ascii_art () {
    echo
    echo -E '  %%%%%%\          %%\        %%%%%%\                    %%%%%%\'
    echo -E ' %%  __%%\         \__|       %%   %%\                        %%\'
    echo -E ' %% /  \__|%%%%%%\ %%\%%%%%%\ %%    %%\ %%%%%%\  %%%%%%%\     %% |'
    echo -E ' \%%%%%%\ %%  __%%\%% %%  _%%\%%    %% %%    %%\%%  _____%%%%%%  |'
    echo -E '  \____%%\%% /  %% %% %% / %% %%    %% %%%%%%%% %% /     %%  ___/'
    echo -E ' %%\   %% %% |  %% %% %% | %% %%   %% /%%   ____%% |     %% |'
    echo -E ' \%%%%%%  %%%%%%%  %% %% | %% %%%%%% / \%%%%%%%\\%%%%%%%\%%%%%%%\'
    echo -E '  \______/%%  ____/\__\__| \__\_____/   \_______|\_______\_______|'
    echo -E '          %% |'
    echo -E '          %% |'
    echo -E '          \__|'
    echo
    echo -e ' Modelling the Phase Field of Spinoidal Decomposition\n'
}

help_message () {
    echo "usage: spindec [-h]"
    echo "               [-c DEBUG]"
    echo "               [-t]"
    echo "               [-e]"
    echo
    echo "options:"
    echo "  -h, --help              show this help message and exit"
    echo "  -c, --compile DEBUG     compile the code with optional debug option"
    echo "                          (default=false)"
    echo "  -t, --test              run automated unit tests"
    echo "  -e, --example           run with example initialisation states"
}

### Argument Parser ###
# Requires util-linux
# getopt options list
options=$(getopt -o c::te::h -l compile::,test,example::,help -- "$@")

# Exit if error code
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# Help and exit if no options given
if [[ $# -lt 1 ]]; then
    ascii_art
    help_message
fi

eval set -- "$options"

# Commands associated with options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c | --compile)
            compile $2
            shift 2
            ;;
        -t | --test)
            auto_test
            shift
            ;;
        -e | --example)
            example
            shift 2
            ;;
        -h | --help)
            ascii_art
            help_message
            break
            ;;
        --) shift;
            break
            ;;
    esac
done
