# Configure tubetk.pth.
# Add a symlink to this file in the directory that results from::
#
#   python -m site --user-site
#
# to make TubeTK Python packages and modules importable.
set( TubeTK_C_PYTHON_EXTENSION_DIR ${CMAKE_CURRENT_BINARY_DIR} )
configure_file( tubetk.pth.in tubetk.pth )

include_directories( ${TubeTK_SOURCE_DIR}/Base/Filtering )

if( TubeTK_USE_NumPy )
  add_subdirectory( tubetk )
endif( TubeTK_USE_NumPy )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
