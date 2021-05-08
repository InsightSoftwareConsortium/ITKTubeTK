# Generate the TCLAPConfig.cmake file in the build tree and configure one
# the installation tree.

# Settings specific to build trees
#
#    tclap include files referenced as "tclap/foo.h"
set(TCLAP_INCLUDE_DIRS_CONFIG ${TCLAP_INCLUDE_DIRS})
set(TCLAP_LIBRARY_DIRS_CONFIG
  ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  )
set(TCLAP_USE_FILE_CONFIG ${TCLAP_BINARY_DIR}/UseTCLAP.cmake)


# Configure TCLAPConfig.cmake for the install tree.
configure_file(${TCLAP_SOURCE_DIR}/TCLAPConfig.cmake.in
               ${TCLAP_BINARY_DIR}/TCLAPConfig.cmake @ONLY)



# Settings specific for installation trees
#
#    tclap include files referenced as "tclap/foo.h"


# TODO- Configure TCLAPInstallConfig.cmake for the install tree.
#configure_file(${TCLAP_SOURCE_DIR}/TCLAPInstallConfig.cmake.in
#               ${TCLAP_BINARY_DIR}/install/TCLAPConfig.cmake @ONLY)
