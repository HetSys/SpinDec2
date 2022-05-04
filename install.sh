#!/bin/bash 
# install.sh
# Compilation/install script for SpinDec

# Exit if something fails
set -e

# Enable recursive globbing 
shopt -s globstar

compile () {
    ### Compilation ###
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

    # Add program files from src 
    prog_files=(./src/*)

    # Binary name and location
    mod_files="./bin"
    compd_file="./bin/spindec"

    # NetCDF flags
    flibs=`nf-config --flibs`
    fflags=`nf-config --fflags`

    # Compile
    echo -e "Compile line: $comp_line $fflags ${prog_files[@]} -J$flibs -o $compd_file\n"
    $comp_line $fflags ${prog_files[@]} -J$flibs -o $compd_file

    # Add binary to $PATH (with some checks)
    while true; do
        read -p 'Add binary to $PATH? [Y/n] ' confirmation

        if [[ "$confirmation" =~ ^[Yy]$ ]] || "[[ $confirmation" == '' ]]; then
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
            echo -e 'Not a valid option\n'
        fi
    done
}

clean () {
    ### Remove compiled binaries ###
    bins=("*/*.mod" "*/*.o" "*/spindec")

    # Check for binaries
    if [[ `echo ${bins[@]}` == "*/*.mod */*.o */spindec" ]]; then 
	echo "No binaries found"
	exit 0
    fi

    # Remove globs from $bins if they aren't found
    for glob in ${bins[@]}; do 
	if [[ "$glob" == *'*'* ]]; then
	    # Write thing to remove from array
	fi 
    done     

    exit 0

    # Remove binaries
    while true; do
	echo -e "Removing the following compiled binaries:\n"
    	ls */* | grep -E 'bin/*.mod|bin/*.o|bin/spindec'
	echo
	read -p 'Proceed? [Y/n] ' confirmation
	
        if [[ $confirmation =~ ^[Yy]$ ]] || [[ $confirmation == '' ]]; then
	    for i in ${bins[@]}; do 
		rm $i
	    done
	    break 
        elif [[ $confirmation =~ ^[Nn]$ ]]; then
            echo 'Binaries not removed'
            break
        else
            echo -e 'Not a valid option\n'
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
    echo "                          (default=none)"
    echo "  -C, --clean 	    remove compiled binaries from repository"
    echo "  -t, --test              run automated unit tests"
    echo "  -e, --example           run with example initialisation states"
}

### Argument Parser ###
# Requires util-linux
# getopt options list
options=$(getopt -o c::Cte::h -l compile::,clean,test,example::,help -- "$@")

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
	-C | --clean)
	    clean
	    ;;
        -t | --test)
            auto_test
            ;;
        -e | --example)
            example
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
