FROM debian:8
MAINTAINER Insight Software Consortium <community@itk.org>

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  git \
  ninja-build \
  wget \
  vim \
  subversion \
  qt4-dev-tools \
  libxt-dev \
  libboost-all-dev \
  valgrind \
  python2.7-dev


WORKDIR /usr/src

ENV Slicer_GIT_TAG 30636dce872c3f51c8f70d0af2b311a5eaf601e9
RUN git clone https://github.com/Slicer/Slicer.git && \
  cd Slicer && \
  git checkout ${Slicer_GIT_TAG} && \
  cd ../ && \
  mkdir Slicer-Release && \
  cd Slicer-Release && \
  cmake \
    -G Ninja \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    ../Slicer && \
  ninja && \
  find . -name '*.o' -delete && \
  cd ../
