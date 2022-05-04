#!/bin/bash 
# install.sh
# Compilation/install script for SpinDec2

# Exit if something fails
set -e

# Enable recursive globbing 
# shopt -s globstar

compile () {
    ### Compilation ###
    # Option for choosing the compiler
    if [[ "$1" == "d" ]] || [[ "$1" == "debug" ]]; then
        comp_line="gfortran -std=f2008 -Wall -fimplicit-none -fcheck=all -Wextra -pedantic -fbacktrace"
    elif [[ -z "$1" ]]; then
        comp_line="gfortran"
    else
        echo -e "$1 is not a valid option for -c/--compile\n"
        help_message
        exit 2
    fi

    # Add program files from src 
    prog_files=(src/*)

    # Binary name and location
    mod_files="bin/"
    compd_file="bin/spindec"

    # NetCDF flags
    flibs=`nf-config --flibs`
    fflags=`nf-config --fflags`

    # Compile
    echo -e "Compile line: $comp_line $fflags ${prog_files[@]} $flibs -J$mod_files -I$mod_files -o $compd_file\n"
    $comp_line $fflags ${prog_files[@]} $flibs -J$mod_files -I$mod_files -o $compd_file

    # Add binary to $PATH (with some checks)
    while true; do
        read -p 'Add binary to $PATH? [Y/n] ' confirmation

        if [[ "$confirmation" =~ ^[Yy]$ ]] || [[ "$confirmation" == '' ]]; then
            if grep -Fq "spindec" $HOME/.bashrc; then
                echo 'The binary is already in your $PATH'
                break
            else
                if [[ "$SHELL" == *"bash"* ]]; then
                    working_dir=`pwd`
                    echo "export PATH=`pwd`/bin/spindec:$PATH" >> $HOME/.bashrc
                    echo 'Binary added to $PATH and written to ~/.bashrc'
                    echo "You will need to restart your shell with 'source ~/.bashrc'"
                    break
                else
                    echo 'Unable to add to $PATH as bash is not your default shell'
                    break
                fi
            fi

        elif [[ "$confirmation" =~ ^[Nn]$ ]]; then
            echo 'Not added to $PATH'
            break
        else
            echo -e 'Not a valid option\n'
        fi
    done
}

clean () {
    ### Remove compiled binaries ###
    bins=("bin/*.mod" "bin/*spindec")

    # Check for binaries
    if [[ `ls */* | grep -E 'mod|spindec'` == '' ]]; then
        echo "No binaries found"
        exit 1
    fi

    # Remove globs from $bins if they aren't found
    for glob in ${bins[@]}; do 
        if [[ "$glob" == *'*'* ]]; then
            bins=("${bins[@]/$glob}")
        fi
    done     

    # Remove binaries
    while true; do
        echo -e "Removing the following files:\n"
        ls */* | grep -E 'mod|spindec'
        echo

        # Ask for confirmation before removing if option provided
        if [[ "$1" == "c" ]] || [[ "$1" == "confirm" ]]; then
            read -p 'Proceed? [Y/n] ' confirmation
        elif [[ -z "$1" ]]; then
            confirmation="y"
        else 
            echo -e "$1 is not a valid option for -C/--clean\n"
            help_message
            exit 2
        fi
      
        if [[ "$confirmation" =~ ^[Yy]$ ]] || [[ "$confirmation" == '' ]]; then
            for file in ${bins[@]}; do
                rm "$file"
            done
            echo "Cleaned successfully"
            break 
        elif [[ "$confirmation" =~ ^[Nn]$ ]]; then
            echo 'Files not removed'
            break
        else
            echo -e '\nNot a valid option'
        fi

    done
}

unit_test_compile () {
    ### Compile unit tests ###

    # Tests directory
    tests=(./test/*)
}

unit_test_run () {
    ### Run unit tests ###
    placeholder="placeholder"
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
    echo -e ' Modelling the Phase Field of Spinoidal Decomposition'
    echo -e ' _________________________________________________________________\n'
}

help_message () {
    echo "usage: spindec [-h]"
    echo "               [-c DEBUG]"
    echo "               [-C CONFIRM]"
    echo "               [-t OPTIONS]"
    echo "               [-e]"
    echo
    echo "options:"
    echo "  -h, --help              show this help message and exit"
    echo "  -c, --compile DEBUG     compile the code with optional debug option"
    echo "                          DEBUG options: [ none | d/debug ] (default=none)"
    echo
    echo "  -C, --clean CONFIRM     remove compiled binaries from repository"
    echo "                          CONFIRM options: [ none | c/confirm ] (default=none)"
    echo
    echo "  -t, --test RUN          run automated unit tests"
    echo "                          RUN options: [ c/compile | r/run | b/both ] (default=both)"
    echo
    echo "  -e, --example           run with example initialisation states"
}

### Argument Parser ###
# Requires util-linux (which should be installed with netcdf)
# getopt options list
options=$(getopt -o c::C::te::h -l compile::,clean::,test,example::,help -- "$@")

# Exit if error code
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# Help and exit if incorrect args given
if [[ $# -lt 1 ]]; then
    ascii_art
    help_message
    exit 0
elif [[ $# -gt 1 ]]; then
    echo -e "Only 1 argument may be specified\n"
    help_message
    exit 2
fi

eval set -- "$options"

# Commands associated with options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c | --compile)
            compile $2
            shift 2
            break
            ;;
        -C | --clean)
            clean $2
            shift 2
            break
            ;;
        -t | --test)
            if [[ "$2" == "b" ]] || [[ "$2" == "both" ]]; then
                unit_test_compile
                unit_test_run
            elif [[ "$2" == "c" ]] || [[ "$2" == "compile" ]]; then
                unit_test_compile
            elif [[ "$2" == "r" ]] || [[ "$2" == "run" ]]; then
                unit_test_run
            else
                echo -e "$2 is not a valid option for -t/--test\n"
                help_message
                exit 2
            fi
            break
            ;;
        -e | --example)
            example
            break
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
