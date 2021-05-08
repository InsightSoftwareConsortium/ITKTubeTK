#!/bin/bash

set -ex

# This is a script to build the project and run the test suite in the base
# Docker container.

# jsoncpp
cd /usr/src
git clone git://github.com/Slicer/jsoncpp.git
mkdir jsoncpp-build && cd $_
cmake -G Ninja \
  -DBUILD_TESTING:BOOL=OFF \
  -DJSONCPP_WITH_TESTS:BOOL=OFF \
  -DJSONCPP_WITH_POST_BUILD_UNITTEST:BOOL=OFF \
  -DJSONCPP_WITH_WARNING_AS_ERROR:BOOL=OFF \
  -DJSONCPP_WITH_PKGCONFIG_SUPPORT:BOOL=OFF \
  -DJSONCPP_WITH_CMAKE_PACKAGE:BOOL=ON \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DBUILD_STATIC_LIBS:BOOL=OFF \
  ../jsoncpp
ninja install

# ParameterSerializer
cd /usr/src
git clone git://github.com/Slicer/ParameterSerializer.git
mkdir ParameterSerializer-build && cd $_
cmake \
  -G Ninja \
  -DBUILD_TESTING:BOOL=OFF \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  /usr/src/ParameterSerializer
ninja

# SlicerExecutionModel
mkdir -p /usr/src/SlicerExecutionModel-build && cd $_
BUILDNAME=sem_use_serializer-on_$1
cmake \
  -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILDNAME:STRING=$BUILDNAME \
  -DParameterSerializer_DIR:PATH=/usr/src/ParameterSerializer-build \
  -DSlicerExecutionModel_USE_SERIALIZER:BOOL=ON \
  -DSlicerExecutionModel_USE_JSONCPP:BOOL=ON \
    /usr/src/SlicerExecutionModel
ctest -VV -D Experimental
