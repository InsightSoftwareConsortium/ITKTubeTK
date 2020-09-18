#!/bin/bash
if [[ ! -v TubeTK_BUILD_DIR ]]; then
    echo 'Please define TubeTK_BUILD_DIR appropriately before running this script' >&2
    return 1
fi

add() {
    [ $# -eq 2 ] || { echo 'add requires two arguments' >&2; exit 1; }
    if [[ -v $1 ]]; then
	temp=${!1}:$2
	eval "$1=\$temp"
    else
	eval "$1=\$2"
    fi
}

for p in ITK-build/Wrapping/Generators/Python \
	 ITK-build/lib \
	 TubeTK-build/ITKModules/TubeTKITK-build/Wrapping/Generators/Python \
	 TubeTK-build/lib; do
    add PYTHONPATH "$TubeTK_BUILD_DIR/$p"
done

unset -f add

export PYTHONPATH
