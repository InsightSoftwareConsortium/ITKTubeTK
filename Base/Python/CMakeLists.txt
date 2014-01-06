# Configure tubetk.pth.
# Add a symlink to this file in the directory that results from::
#
#   python -m site --user-site
#
# to make TubeTK Python packages and modules importable.
set( TubeTK_C_PYTHON_EXTENSION_DIR ${CMAKE_CURRENT_BINARY_DIR} )
configure_file( tubetk.pth.in tubetk.pth )
if( BUILD_TESTING AND NOT TubeTK_USE_SLICER )
  # Install into the testing virtualenv.
  execute_process( COMMAND ${PYTHON_TESTING_EXECUTABLE}
      -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    RESULT_VARIABLE _result_variable
    ERROR_VARIABLE _error_variable
    OUTPUT_VARIABLE site_packages )
  if( ${_result_variable} )
    message( SEND_ERROR "Could not get virtual packages directory:
      ${_result_variable}
      ${_error_variable}" )
  else( ${_result_variable} )
    message( STATUS "Detected testing virtualenv packages directory: ${site_packages}" )
    string( STRIP "${site_packages}" site_packages )
    configure_file( tubetk.pth.in "${site_packages}/tubetk.pth" )
  endif( ${_result_variable} )
endif( BUILD_TESTING AND NOT TubeTK_USE_SLICER )

include_directories( ${TubeTK_SOURCE_DIR}/Base/Filtering )

add_subdirectory( tubetk )

if( BUILD_TESTING AND NOT TubeTK_USE_SLICER )
  add_subdirectory( Testing )
endif( BUILD_TESTING AND NOT TubeTK_USE_SLICER )
