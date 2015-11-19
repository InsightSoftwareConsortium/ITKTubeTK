#!/bin/sh

# Extract tag
script_dir="`cd $(dirname $0); pwd`"
image_tag="`echo $script_dir | rev | cut -d'/' -f 1 | rev | cut -d'-' -f 2-`"
# Run common run script
$script_dir/../Docker/run.sh $image_tag
