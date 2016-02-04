find_package( ConsistentPython REQUIRED )

if( NOT NUMPY_INCLUDE_DIR )
  find_package( PythonInterp )
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE _numpy_include
    OUTPUT_STRIP_TRAILING_WHITESPACE )
  find_path( NUMPY_INCLUDE_DIR numpy/arrayobject.h
    HINTS ${_numpy_include} )
endif( NOT NUMPY_INCLUDE_DIR )

include_directories( ${NUMPY_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS} )

python_add_module( _tubetk_numpy _tubetk_numpy.cxx )
set_target_properties( _tubetk_numpy
  PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

target_link_libraries( _tubetk_numpy
  ${ITK_LIBRARIES} TubeTKFiltering )
  
file( COPY ./ DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
  FILES_MATCHING PATTERN "*.py" )
