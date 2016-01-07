#!/bin/sh

# Extract tag
script_dir="`cd $(dirname $0); pwd`"
image_tag="`echo $script_dir | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 2- | rev | cut -d'_' -f 2- | rev`"
# Run common push script
$script_dir/../CircleCI_Docker/push.sh $image_tag
