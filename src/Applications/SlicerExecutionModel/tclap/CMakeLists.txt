project(TCLAP LANGUAGES CXX)

# --------------------------------------------------------------------------
# Install headers
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  set(${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT ON)
endif()
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  file(GLOB headers "${TCLAP_SOURCE_DIR}/include/tclap/*.h")
  # Note "COPYING" is a filename containing the license file
  install(
    FILES ${headers} COPYING
    DESTINATION include/tclap COMPONENT Development
    )
endif()

# --------------------------------------------------------------------------
# Install support files
# --------------------------------------------------------------------------
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES
    ${TCLAP_BINARY_DIR}/UseTCLAP.cmake
    ${TCLAP_BINARY_DIR}/install/TCLAPConfig.cmake
    DESTINATION lib/tclap COMPONENT Development
    )
endif()

# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
  CACHE INTERNAL "${PROJECT_NAME} include dirs" FORCE)

# --------------------------------------------------------------------------
# Configure support file
# --------------------------------------------------------------------------
include(GenerateTCLAPConfig.cmake)

configure_file(
  ${TCLAP_SOURCE_DIR}/UseTCLAP.cmake.in
  ${TCLAP_BINARY_DIR}/UseTCLAP.cmake COPYONLY
  )
