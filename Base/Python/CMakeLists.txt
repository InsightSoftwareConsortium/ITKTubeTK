# Configure tubetk.pth.
# Add a symlink to this file in the directory that results from::
#
#   python -m site --user-site
#
# to make TubeTK Python packages and modules importable.
set( TubeTK_C_PYTHON_EXTENSION_DIR ${CMAKE_CURRENT_BINARY_DIR} )
configure_file( tubetk.pth.in tubetk.pth )
if( BUILD_TESTING )
  # Install into the testing virtualenv.
  execute_process( COMMAND ${PYTHON_TESTING_EXECUTABLE}
      -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE site_packages )
  string( STRIP ${site_packages} site_packages )
  configure_file( tubetk.pth.in "${site_packages}/tubetk.pth" )
endif( BUILD_TESTING )

include_directories( ${TubeTK_SOURCE_DIR}/Base/Filtering )

add_subdirectory( tubetk )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
