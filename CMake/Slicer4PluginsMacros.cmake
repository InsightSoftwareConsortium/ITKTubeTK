#-----------------------------------------------------------------------------
# Set the default output paths for one or more plugins/CLP
#
macro(slicer4_set_plugins_output_path)
  set_target_properties(${ARGN}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${Slicer4_INSTALL_PLUGINS_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${Slicer4_INSTALL_PLUGINS_LIB_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${Slicer4_INSTALL_PLUGINS_LIB_DIR}"
    )
endmacro(slicer4_set_plugins_output_path)

#-----------------------------------------------------------------------------
# Install one or more plugins to the default plugin location
#
macro(slicer4_install_plugins)
  install(TARGETS ${ARGN}
    RUNTIME DESTINATION ${Slicer4_INSTALL_PLUGINS_BIN_DIR} COMPONENT RuntimeLibraries
    LIBRARY DESTINATION ${Slicer4_INSTALL_PLUGINS_LIB_DIR} COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION ${Slicer4_INSTALL_PLUGINS_LIB_DIR} COMPONENT RuntimeLibraries
    )
endmacro(slicer4_install_plugins)
