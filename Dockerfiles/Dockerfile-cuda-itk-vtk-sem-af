FROM nvidia/cuda:8.0-cudnn5-devel-ubuntu16.04
MAINTAINER Deepak Roy Chittajallu <deepk.chittajallu@kitware.com>

# Install system pre-requisites
RUN apt-get update && \
    apt-get install -y \
    build-essential wget git \
    make cmake cmake-curses-gui ninja-build \
    libxt-dev libgl1-mesa-dev libcupti-dev \
    libboost-all-dev libfftw3-dev liblapack-dev liblapacke-dev libopenblas-dev \
    libfontconfig1-dev libfreeimage-dev xorg-dev \
    ocl-icd-opencl-dev opencl-headers && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    echo "50.58.123.189 data.kitware.com" >> /etc/hosts && \
    echo "50.58.123.181 midas3.kitware.com" >> /etc/hosts

# Setting up symlinks for libcuda and OpenCL ICD
RUN ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/lib/libcuda.so.1 && \
    ln -s /usr/lib/libcuda.so.1 /usr/lib/libcuda.so && \
    mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd && \
    echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf && \
    ldconfig
ENV PATH=/usr/local/nvidia/bin:/usr/local/cuda/bin:${PATH}

# Libraries build path
ENV BUILD_PATH /build

# Install miniconda
RUN mkdir -p $BUILD_PATH && \
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh \
    -O $BUILD_PATH/install_miniconda.sh && \
    bash $BUILD_PATH/install_miniconda.sh -b -p $BUILD_PATH/miniconda && \
    rm $BUILD_PATH/install_miniconda.sh && \
    chmod -R +r $BUILD_PATH && \
    chmod +x $BUILD_PATH/miniconda/bin/python
ENV PATH $BUILD_PATH/miniconda/bin:${PATH}

# Install CMake
ENV CMAKE_ARCHIVE_SHA256 10ca0e25b7159a03da0c1ec627e686562dc2a40aad5985fd2088eb684b08e491
ENV CMAKE_VERSION_MAJOR 3
ENV CMAKE_VERSION_MINOR 8
ENV CMAKE_VERSION_PATCH 1
ENV CMAKE_VERSION ${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}.${CMAKE_VERSION_PATCH}
RUN cd $BUILD_PATH && \
  wget https://cmake.org/files/v${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
  hash=$(sha256sum ./cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz | awk '{ print $1 }') && \
  [ $hash = "${CMAKE_ARCHIVE_SHA256}" ] && \
  tar -xzvf cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
  rm cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz
ENV PATH $BUILD_PATH/cmake-${CMAKE_VERSION}-Linux-x86_64/bin:${PATH}

# Disable "You are in 'detached HEAD' state." warning
RUN git config --global advice.detachedHead false

# Download/configure/build/install ITK
ENV ITK_GIT_TAG v4.11.1
ENV ITK_BUILD_DIR $BUILD_PATH/ITK-build
RUN cd $BUILD_PATH && \
    git clone --depth 1 -b ${ITK_GIT_TAG} https://github.com/InsightSoftwareConsortium/ITK.git && \
    mkdir ITK-build && \
    cd ITK-build && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DBUILD_EXAMPLES:BOOL=OFF \
        -DBUILD_TESTING:BOOL=OFF \
        -DITKV3_COMPATIBILITY:BOOL=ON \
        -DITK_BUILD_DEFAULT_MODULES:BOOL=ON \
        -DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON \
        -DITK_LEGACY_REMOVE:BOOL=OFF \
        -DITK_LEGACY_SILENT:BOOL=ON \
        -DITK_WRAP_PYTHON:BOOL=ON \
        -DModule_MinimalPathExtraction:BOOL=ON \
        -DKWSYS_USE_MD5:BOOL=ON \
        -DModule_ITKReview:BOOL=ON \
        ../ITK && \
    ninja && \
    cp Wrapping/Generators/Python/WrapITK.pth $BUILD_PATH/miniconda/lib/python2.7/site-packages && \
    python -c "import itk" && \
    find . -name '*.o' -delete && \
    find ../ITK* -depth -name .git -exec rm -rf '{}' \;

# Download/configure/build/install VTK
ENV VTK_GIT_TAG v7.1.1
ENV VTK_BUILD_DIR $BUILD_PATH/VTK-build
RUN cd $BUILD_PATH && \
    git clone --depth 1 -b ${VTK_GIT_TAG} https://github.com/Kitware/VTK.git && \
    mkdir VTK-build && \
    cd VTK-build && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DBUILD_EXAMPLES:BOOL=OFF \
        -DBUILD_TESTING:BOOL=OFF \
        -DVTK_LEGACY_REMOVE:BOOL=ON \
        -DVTK_WRAP_PYTHON:BOOL=ON \
        ../VTK && \
    ninja && \
    echo "${VTK_BUILD_DIR}/lib" > $BUILD_PATH/miniconda/lib/python2.7/site-packages/WrapVTK.pth && \
    echo "${VTK_BUILD_DIR}/Wrapping/Python" >> $BUILD_PATH/miniconda/lib/python2.7/site-packages/WrapVTK.pth && \
    python -c "import vtk" && \
    find . -name '*.o' -delete && \
    find ../VTK* -depth -name .git -exec rm -rf '{}' \;

# Download/configure/build/install SlicerExecutionModel
ENV SEM_GIT_TAG master
ENV SEM_BUILD_DIR $BUILD_PATH/SEM-build
RUN cd $BUILD_PATH && \
    git clone --depth 1 -b ${SEM_GIT_TAG} https://github.com/Slicer/SlicerExecutionModel.git SEM && \
    mkdir SEM-build && cd SEM-build && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DBUILD_TESTING:BOOL=OFF \
        -DITK_DIR:PATH=$BUILD_PATH/ITK-build \
        ../SEM && \
    ninja && \
    find . -name '*.o' -delete && \
    find ../SEM* -depth -name .git -exec rm -rf '{}' \;

# Download and install slicer_cli_web
ENV SLICER_CLI_WEB_GIT_TAG master
RUN cd $BUILD_PATH && \
    pip install --upgrade 'git+https://github.com/cdeepakroy/ctk-cli' && \
    git clone --depth 1 -b ${SLICER_CLI_WEB_GIT_TAG} https://github.com/girder/slicer_cli_web.git && \
    cd slicer_cli_web && \
    find . -depth -name .git -exec rm -rf '{}' \;

# Download/configure/build/install GLFW
ENV GLFW_GIT_TAG 3.2.1
RUN cd $BUILD_PATH && \
    git clone --depth 1 -b ${GLFW_GIT_TAG} https://github.com/glfw/glfw.git && \
    mkdir glfw-build && cd glfw-build && \
    cmake \
          -G Ninja \
          -DCMAKE_INSTALL_PREFIX=/usr \
          ../glfw && \
    ninja install && \
    cd .. && rm -rf glfw*

# Download/configure/build/install arrayfire
ENV AF_GIT_TAG devel
ENV AF_PATH=/usr/local/arrayfire AF_DISABLE_GRAPHICS=1
RUN cd $BUILD_PATH && \
    git clone --recursive --depth 1 -b ${AF_GIT_TAG} https://github.com/arrayfire/arrayfire.git && \
    mkdir arrayfire-build && cd arrayfire-build && \
    cmake \
          -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_CPU=ON \
          -DBUILD_CUDA=ON \
          -DBUILD_OPENCL=OFF \
          -DBUILD_UNIFIED=ON \
          -DBUILD_GRAPHICS=OFF \
          -DBUILD_NONFREE=OFF \
          -DBUILD_EXAMPLES=OFF \
          -DBUILD_TEST=OFF \
          -DBUILD_DOCS=OFF \
          -DUSE_FREEIMAGE_STATIC=OFF \
          ../arrayfire && \
    make && make install && \
    echo "${AF_PATH}/lib" >> /etc/ld.so.conf.d/arrayfire.conf && \
    echo "/usr/local/cuda/nvvm/lib64" >> /etc/ld.so.conf.d/arrayfire.conf && \
    ldconfig && \
    cd $BUILD_PATH && rm -rf arrayfire*
