#!/bin/sh

# Extract tag
script_dir="`cd $(dirname $0); pwd`"
image_tag="`echo $script_dir | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 2- | rev | cut -d'_' -f 2- | rev`"
# Run common run script
$script_dir/../CircleCI_Docker/run.sh $image_tag
