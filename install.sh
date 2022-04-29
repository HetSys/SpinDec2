#!/bin/bash 

# install.sh
# Compilation/install script for SpinDec
# Author: Dylan Morgan dylan.morgan@warwick.ac.uk 

### Compilation ###
# function compile{
#     # Compiler 
#     if [[ ]]
# }

echo "test"

### Argument Parser ###
# Requires util-linux
OPTIONS=$(getopt -o c::dteh -l compile::,debug,test,example,help -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi


eval set -- "$OPTIONS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c | --compile)
            echo "compile $2"
            shift 2
            ;;
        -d | --debug)
            echo "debug"
            shift
            ;;
        -t | --test)
            echo "test"
            shift
            ;;
        -e | --example)
            echo "example $2"
            shift 2
            ;;
        -h | --help) 
            echo "help"
            ;;
        --) shift;
            break
            ;;
    esac
done

