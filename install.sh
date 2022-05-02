#!/bin/bash 
# install.sh
# Compilation/install script for SpinDec

# Exit if something fails
set -e

### Compilation ###
compile () {
    # Option for choosing the compiler
    if [[ $1 == "d" ]] || [[ $1 == "debug" ]]; then
        comp_line="gfortran -g -std=f2008 -Wall -fimplicit-none -fcheck=all -Wextra -pedantic -fbacktrace"
    elif [[ -z $1 ]]; then
        comp_line="gfortran -g"
    else
        echo -e "Not a valid option for -c/--compile\n"
        help_message
        exit 2
    fi

    # Add program files from prog_files.txt
    input_file="./prog_files.txt"
    prog_files=()
    while read -r line; do
        prog_files+=("$line")
    done < "$input_file"

    echo "prog file: $prog_files"

    # Binary name and location
    compd_file="./bin/spindec"

    # NetCDF flags
    flibs=`nf-config --flibs`
    fflags=`nf-config --fflags`

    echo -e "Compile line: $comp_line\n"

    # Compile
    $comp_line $fflags $prog_files $flibs -o $compd_file

    # Add binary to $PATH (with some checks)
    while true; do
        if [[ $debug == true ]]; then
            read -p -n 1 -r 'Add to $PATH? [Y/n] ' confirmation

            if [[ $confirmation =~ ^[Yy]$ ]] || [[ $confirmation == '' ]]; then
                if grep -Fq "spindec" $HOME/.bashrc; then
                    echo 'The binary is already in your $PATH'
                    break
                else
                    if [[ "$SHELL" == *"bash"* ]]; then
                        working_dir=`pwd`
                        echo "export PATH=`pwd`/bin/spindec:$PATH" >> $HOME/.bashrc
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
}

help_message () {
    echo "Haha good luck"
}

### Argument Parser ###
# Requires util-linux
options=$(getopt -o c::te::h -l compile::,test,example::,help -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

eval set -- "$options"

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
