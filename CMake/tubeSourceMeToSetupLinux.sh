#! /bin/bash
SOURCE="${BASH_SOURCE[0]}"
DIR="$( dirname "$SOURCE" )"
while [ -h "$SOURCE" ]
do
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
  DIR="$( cd -P "$( dirname "$SOURCE"  )" && pwd )"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
echo $DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/lib/TubeTK/plugins
export PATH=$PATH:$PWD/bin:$DIR/lib/TubeTK/plugins
