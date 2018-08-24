#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

die() {
  echo "Error: $@" 1>&2
  exit 1;
}

if [ ! $1 ];
then
  die "Empty Image Tag "
fi

lower_case_tag="`echo $1 | tr "[:upper:]" "[:lower:]" `"

docker run \
  --rm \
  -v $script_dir/../../..:/usr/src/ITKTubeTK \
    kitwaremedical/itktubetk:$lower_case_tag \
      /usr/src/ITKTubeTK/CMake/CircleCI/CircleCI_TubeTK_Experimental.sh $lower_case_tag
