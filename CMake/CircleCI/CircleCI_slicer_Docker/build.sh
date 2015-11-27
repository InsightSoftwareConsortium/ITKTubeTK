#!/bin/sh

# Extract tag
script_dir="`cd $(dirname $0); pwd`"
image_tag="`echo $script_dir | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 2- | rev | cut -d'_' -f 2- | rev`"
# Run common build script
$script_dir/../CircleCI_Docker/build.sh $image_tag $script_dir
