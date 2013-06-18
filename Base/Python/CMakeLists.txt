# Configure tubetk.pth.
# Add a symlink to this file in the directory that results from::
#
#   python -m site --user-site
#
# to make TubeTK Python packages and modules importable.
set( TubeTK_PURE_PYTHON_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
configure_file( tubetk.pth.in tubetk.pth )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
