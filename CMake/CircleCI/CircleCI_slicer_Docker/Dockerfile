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

ENV Slicer_GIT_TAG 0fa5202aadf70664dbdbf6aa0f1f5265a573333b
RUN git clone https://github.com/Slicer/Slicer.git && \
  cd Slicer && \
  git checkout ${Slicer_GIT_TAG} && \
  cd ../ && \
  mkdir Slicer-Release && \
  cd Slicer-Release && \
  cmake \
    -G "Unix Makefiles" \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    ../Slicer && \
  make -j 8 && \
  find . -name '*.o' -delete && \
  cd ../
