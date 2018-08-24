#!/bin/sh

die() {
  echo "Error: $@" 1>&2
  exit 1;
}

if [ ! $1 ];
then
  die "Empty Image Tag "
fi

if [ ! $2 ];
then
  die "Empty Path to Dockerfile "
fi

lower_case_tag="`echo $1 | tr "[:upper:]" "[:lower:]" `"
docker build -t kitwaremedical/itktubetk:$lower_case_tag $2
