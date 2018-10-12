FROM kitwaremedical/cuda-itk-vtk-sem-af
MAINTAINER Deepak Roy Chittajallu <deepk.chittajallu@kitware.com>

# Download/configure/build/install TubeTK
ENV TubeTK_SRC_DIR=${BUILD_PATH}/TubeTK
ENV TubeTK_BUILD_DIR=${BUILD_PATH}/TubeTK-build
RUN mkdir -p ${TubeTK_SRC_DIR} && mkdir -p ${TubeTK_BUILD_DIR}
COPY . ${TubeTK_SRC_DIR}
WORKDIR ${TubeTK_BUILD_DIR}

RUN cmake \
        -G Ninja \
        -DTubeTK_BUILD_APPLICATIONS:BOOL=ON \
        -DTubeTK_USE_PYTHON:BOOL=ON \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DBUILD_TESTING:BOOL=ON \
        -DTubeTK_BUILD_USING_SLICER:BOOL=OFF \
        -DTubeTK_USE_ARRAYFIRE:BOOL=ON \
        -DTubeTK_USE_EXAMPLES_AS_TESTS:BOOL=OFF \
        -DTubeTK_USE_BOOST:BOOL=OFF \
        -DTubeTK_USE_PYQTGRAPH:BOOL=OFF \
        -DTubeTK_USE_NUMPY_STACK:BOOL=OFF \
        -DITK_DIR:PATH=$ITK_BUILD_DIR \
        -DTubeTK_USE_VTK:BOOL=ON \
        -DVTK_DIR:PATH=$VTK_BUILD_DIR \
        -DUSE_SYSTEM_SLICER_EXECUTION_MODEL:BOOL=ON \
        -DSlicerExecutionModel_DIR:PATH=$SEM_BUILD_DIR \
        ../TubeTK && \
    ninja

RUN cd ${TubeTK_BUILD_DIR} && \
    cp TubeTK-build/ITKModules/TubeTKITK-build/Wrapping/Generators/Python/WrapITK.pth $BUILD_PATH/miniconda/lib/python2.7/site-packages/WrapITKTubeTK.pth && \
    cp TubeTK-build/PythonModules/tubetk.pth $BUILD_PATH/miniconda/lib/python2.7/site-packages && \
    find . -name "*.o" -delete && \
    find ../TubeTK* -depth -name .git -exec rm -rf '{}' \; && \
    cd $TubeTK_SRC_DIR && \
    pip install -r requirements-ml.txt && \
    cd $TubeTK_SRC_DIR/Applications && \
    python generate_slicer_cli_list_json.py
ENV PATH="${TubeTK_BUILD_DIR}/TubeTK-build/bin:${PATH}"

# Set workdir to TubeTK Applications
WORKDIR $TubeTK_SRC_DIR/Applications

# Test slicer_cli_web entrypoint
RUN python /build/slicer_cli_web/server/cli_list_entrypoint.py --list_cli

# Set entrypoint
ENTRYPOINT ["/build/miniconda/bin/python", "/build/slicer_cli_web/server/cli_list_entrypoint.py"]
