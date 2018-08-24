#!/bin/sh

die() {
  echo "Error: $@" 1>&2
  exit 1;
}

if [ ! $1 ];
then
  die "Empty Image Tag "
fi

lower_case_tag="`echo $1 | tr "[:upper:]" "[:lower:]" `"

docker push kitwaremedical/itktubetk:$lower_case_tag
