#!/bin/bash
# install_spectral.sh
# Compilation/install script for SpinDec2

# Exit if something fails
set -e

# Enable recursive globbing
# shopt -s globstar

compile () {
    ### Compilation ###
    # Option for choosing the compiler
    if [[ "$1" == "d" ]] || [[ "$1" == "debug" ]]; then
        comp_line="gfortran -fopenmp -O2 -std=f2008 -Wall -fimplicit-none -fcheck=all -Wextra -pedantic -fbacktrace"
    elif [[ -z "$1" ]]; then
        comp_line="gfortran -fopenmp -O2"
    else
        echo -e "$1 is not a valid option for -c/--compile\n"
        help_message
        exit 2
    fi

    #fftw_dir=$(find . -maxdepth 1 -name "fftw*" -type d)

    #if [[ "$fftw_dir" == *"fftw-3"* ]]; then
    #    echo "Found fftw-3 library"
    #    echo "Renaming to fftw3"
    #    mv "$fftw_dir" fftw3
    #elif [[ "$fftw_dir" == *"fftw3" ]]; then
    #    echo "fftw3 library already present"
    #elif [[ "$fftw_dir" == "" ]]; then
    #    echo "Required fftw-3 library not found"
    #    echo "Ensure this is installed before attempting compilation"
    #    echo "Exiting compilation"
    #    exit 1
    #fi

    # Add program files from src
    prog_files=(src/*.f90)

    # Move main to last item in array
    #main="src/spectral_main.f90"
    main="src/main.f90"
    prog_files=("${prog_files[@]/$main}")
    prog_files+=("$main")
    #prog_files=("${prog_files[@]/$mainb}")

    # Binary name and location
    bin_files="bin/"
    obj_files="bin/*.o"
    compd_file="bin/spindec"

    # NetCDF flags
    flibs="`nf-config --flibs` -lfftw3_omp -lfftw3 -lm -I/usr/include"
    fflags=`nf-config --fflags`

    # Compile
    echo "Compile line:"
    for file in ${prog_files[@]}; do
        file_nx=${file/.f90/}
        file_nx=${file_nx/src\//}  # Remove extension and upper file directory
        $comp_line -J$bin_files -c $file $fflags $flibs -o $bin_files$file_nx.o
        echo "$comp_line -J$bin_files -c $file $fflags $flibs -o $bin_files$file_nx.o"
    done

    $comp_line -o $compd_file $obj_files $fflags $flibs
    echo "$comp_line -o $compd_file $obj_files $fflags $flibs"

    # Old line:
    # $comp_line $fflags ${prog_files[@]} $flibs -J$mod_files -I$mod_files -o $compd_file

    # Don't prompt to add to $PATH if debug option specified
    if [[ "$1" == "d" ]] || [[ "$1" == "debug" ]]; then
        exit 0
    fi

    # Don't prompt to add to $PATH if already in path
    if grep -Fq "SpinDec2/bin" $HOME/.bashrc; then
        echo -e '\nBinary already in $PATH'
        exit 0
    fi

    # Add binary to $PATH if bash is default $SHELL
    while true; do
        echo
        read -p 'Add binary to $PATH? [Y/n] ' confirmation

        if [[ "$confirmation" =~ ^[Yy]$ ]] || [[ "$confirmation" == '' ]]; then
            if [[ "$SHELL" == *"bash"* ]]; then
                working_dir=`pwd`
                path_var='$PATH'
                echo '' >> $HOME/.bashrc
                echo "export PATH=$working_dir/bin/:$path_var" >> $HOME/.bashrc
                echo 'Binary added to $PATH and written to ~/.bashrc'
                echo "For changes to take effect, use the command 'source ~/.bashrc'"
                break
            else
                echo 'Unable to add to $PATH as bash is not your default shell'
                break
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
    bins=("bin/*.mod" "bin/*.o" "bin/*spindec")

    # Check for binaries
    if [[ `ls bin/* | grep -E 'mod|[.]o|spindec'` == '' ]]; then
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
        ls bin/* | grep -E 'mod|[.]o|spindec'
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
    # Compile line
    # Only the debug compile line from above to be used
    comp_line="gfortran -fopenmp -std=f2008 -Wall -fimplicit-none -fcheck=all -Wextra -pedantic -fbacktrace"

    # f90 file directories
    test_files=(test/*.f90)
    src_files=(src/*.f90)

    # Move test main to last item in array
    test_main="test/test.f90"
    test_files=("${test_files[@]/$test_main}")
    test_files+=("$test_main")

    # Remove actual main from src
    src_main="src/main.f90"
    src_files=("${src_files[@]/$src_main}")

    #src_mainb="src/spectral_main.f90"
    #src_files=("${src_files[@]/$src_mainb}")

    # Binary name and location
    bin_files="test/test_bin/"
    obj_files="test/test_bin/*.o"
    compd_file="test/test_bin/test_spindec"

    # NetCDF flags
    flibs="`nf-config --flibs` -lfftw3_omp -lfftw3 -lm"
    fflags=`nf-config --fflags`

    # C H O N K Y  compilation
    echo "Compile line:"

    # Files in src
    for file in ${src_files[@]}; do
        file_nx=${file/.f90/}
        file_nx=${file_nx/src\//}  # Remove extension and upper file directory
        $comp_line -J$bin_files -c $file $fflags $flibs -o $bin_files$file_nx.o
        echo "$comp_line -J$bin_files -c $file $fflags $flibs -o $bin_files$file_nx.o"
    done

    # Files in test
    for file in ${test_files[@]}; do
        file_nx=${file/.f90/}
        file_nx=${file_nx/test\//}  # Remove extension and upper file directory
        $comp_line -J$bin_files -c $file $fflags $flibs -o $bin_files$file_nx.o
        echo "$comp_line -J$bin_files -c $file $fflags $flibs -o $bin_files$file_nx.o"
    done

    $comp_line -o $compd_file $obj_files $fflags $flibs
    echo "$comp_line -o $compd_file $obj_files $fflags $flibs"
}

unit_test_clean () {
    ### Remove compiled binaries from unit testing ###
    bins=("test/test_bin/*.mod" "test/test_bin/*.o" "test/test_bin/*test_spindec")

    # Check for binaries
    if [[ `ls test/test_bin/* | grep -E 'mod|[.]o|test_spindec'` == '' ]]; then
        echo "No binaries found"
        exit 1
    fi

    # Remove globs from $bins if they aren't found
    for glob in ${bins[@]}; do
        if [[ "$glob" == *'*'* ]]; then
            bins=("${bins[@]/$glob}")
        fi
    done

    echo -e "Removing the following files:\n"
    ls test/test_bin/* | grep -E 'mod|[.]o|test_spindec'
    echo

    for file in ${bins[@]}; do
        rm "$file"
    done
    echo "Cleaned successfully"
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
    echo -e ' Modelling Spinoidal Decomposition Using a Phase Field Approach'
    echo -e ' _________________________________________________________________\n'
}

help_message () {
    echo "usage: spindec [-h]"
    echo "               [-c DEBUG]"
    echo "               [-C CONFIRM]"
    echo "               [-t RUN]"
    echo "               [-T]"
    echo
    echo "options:"
    echo "  -h, --help              show this help message and exit"
    echo "  -c, --compile DEBUG     compile the code with optional debug option"
    echo "                          optional DEBUG arguments: [ none | d/debug ] (default=none)"
    echo
    echo "  -C, --clean CONFIRM     remove compiled binaries from repository"
    echo "                          optional CONFIRM arguments: [ none | c/confirm ] (default=none)"
    echo
    echo "  -t, --test RUN          run automated unit tests"
    echo "                          required RUN arguments: [ c/compile | r/run | b/both ]"
    echo
    echo "  -T, --test-clean        clean test binaries"
}

### Argument Parser ###
# Requires util-linux (which should be installed automatically with netcdf)
# getopt options list
options=$(getopt -o c::C::t:Th -l compile::,clean::,test:,test-clean,help -- "$@")

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
                cd test/test_bin/
                echo
                ./test_spindec
                cd ../../
            elif [[ "$2" == "c" ]] || [[ "$2" == "compile" ]]; then
                unit_test_compile
            elif [[ "$2" == "r" ]] || [[ "$2" == "run" ]]; then
                cd test/test_bin/ && ./test_spindec
                cd ../../
            else
                echo -e "$2 is not a valid option for -t/--test\n"
                help_message
                exit 2
            fi
            break
            ;;
        -T | --test-clean)
            unit_test_clean
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