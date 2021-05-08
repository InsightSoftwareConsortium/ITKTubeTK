#-----------------------------------------------------------------------------
project(GenerateCLP LANGUAGES C CXX)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# --------------------------------------------------------------------------
# Prerequisites
# --------------------------------------------------------------------------
find_package(ModuleDescriptionParser REQUIRED)
include(${ModuleDescriptionParser_USE_FILE})

find_package(TCLAP REQUIRED)
include(${TCLAP_USE_FILE})

#
# ITK
#
set(${PROJECT_NAME}_ITK_COMPONENTS
  ${ModuleDescriptionParser_ITK_COMPONENTS}
  )
find_package(ITK 4.3 COMPONENTS ${${PROJECT_NAME}_ITK_COMPONENTS} REQUIRED)
set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # See Slicer/Libs/ITKFactoryRegistration/CMakeLists.txt
include(${ITK_USE_FILE})

# --------------------------------------------------------------------------
# Option(s)
# --------------------------------------------------------------------------
if(NOT DEFINED BUILD_SHARED_LIBS)
  option(BUILD_SHARED_LIBS "Build with shared libraries." ON)
endif()
if(NOT DEFINED GENERATECLP_USE_MD5)
  option(GENERATECLP_USE_MD5 "Enabling that option will teach GenerateCLP to overwrite output file only if needed. Make sure to build ITK with KWSYS_USE_MD5 enabled." OFF)
endif()
if(NOT DEFINED GenerateCLP_USE_JSONCPP)
  option(GenerateCLP_USE_JSONCPP "Automatic flags are added to serialize and deserialized the parameters." OFF)
  mark_as_advanced(GenerateCLP_USE_JSONCPP)
endif()
if(GenerateCLP_USE_JSONCPP)
  include_directories(${JsonCpp_INCLUDE_DIRS})
endif()

if(NOT DEFINED GenerateCLP_USE_SERIALIZER)
  option(GenerateCLP_USE_SERIALIZER "Automatic flags are added to dump the json-ld schema for the parameters." OFF)
  mark_as_advanced(GenerateCLP_USE_SERIALIZER)
endif()
if(GenerateCLP_USE_SERIALIZER)
  if(NOT ModuleDescriptionParser_USE_SERIALIZER)
    message(SEND_ERROR "ModuleDescriptionParser_USE_SERIALIZER must be ON if GenerateCLP_USE_SERIALIZER is ON.")
  endif()
  include_directories(${ParameterSerializer_INCLUDE_DIRS})
endif()

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
include(CTest)

#-----------------------------------------------------------------------------
# Output directories.
#-----------------------------------------------------------------------------
foreach(type LIBRARY RUNTIME ARCHIVE)
  # Make sure the directory exists
  if(DEFINED GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY
     AND NOT EXISTS ${GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY})
    message(FATAL_ERROR "GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY is set to a non-existing directory [${GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY}]")
  endif()

  if(NOT DEFINED GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY)
    set(GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY ${GenerateCLP_BINARY_DIR}/bin)
  endif()
  set(CMAKE_${type}_OUTPUT_DIRECTORY ${GenerateCLP_CMAKE_${type}_OUTPUT_DIRECTORY} CACHE INTERNAL "Single output directory for building all libraries.")
endforeach()

# --------------------------------------------------------------------------
# Configure header(s)
# --------------------------------------------------------------------------
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GenerateCLPConfig.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/GenerateCLPConfig.h
  )

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  )

# --------------------------------------------------------------------------
# List of targets to export
# --------------------------------------------------------------------------
set(targets_to_export)

# --------------------------------------------------------------------------
# Build executable
# --------------------------------------------------------------------------
set(executable_name GenerateCLP)
set(GENERATECLP_SOURCE
  GenerateCLP.cxx
  )
sem_add_executable(${executable_name} ${GENERATECLP_SOURCE})
list(APPEND targets_to_export ${executable_name})

target_link_libraries(${executable_name}
  # Appearently windows does not like static libs mixed with shared libs ModuleDescriptionParser-static
  # A different solution will have to be investigated for makeing GenerateCLP work without
  # shared libs.
  ${JsonCpp_LIBRARIES}
  ${ParameterSerializer_LIBRARIES}
  ModuleDescriptionParser
  ${ITK_LIBRARIES}
  )

# --------------------------------------------------------------------------
# Install executable
# --------------------------------------------------------------------------
if(NOT DEFINED GenerateCLP_INSTALL_BIN_DIR)
  set(GenerateCLP_INSTALL_BIN_DIR bin)
endif()
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  set(${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT ON)
endif()
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(TARGETS GenerateCLP RUNTIME
    DESTINATION ${GenerateCLP_INSTALL_BIN_DIR} COMPONENT Development
    )
endif()

# --------------------------------------------------------------------------
# Install support files
# --------------------------------------------------------------------------
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES ${GenerateCLP_BINARY_DIR}/GenerateCLPConfig.cmake_install
    DESTINATION lib/GenerateCLP
    COMPONENT Development
    RENAME GenerateCLPConfig.cmake
    )
  install(FILES ${GenerateCLP_BINARY_DIR}/UseGenerateCLP.cmake_install
    DESTINATION lib/GenerateCLP
    COMPONENT Development
    RENAME UseGenerateCLP.cmake
    )
endif()

# --------------------------------------------------------------------------
# Launcher
# --------------------------------------------------------------------------
# Enable shared link forwarding support if needed.
# This is required so that GenerateCLP can be run from either a build or an installed tree,
# where the binary is actually not in the same directory as the shared
# libraries it depends on.
#

#Add additional forward paths
set(additional_forward_paths_build "")
foreach(lib_path ${SlicerExecutionModel_ADDITIONAL_FORWARD_PATHS_BUILD})
  string(APPEND additional_forward_paths_build
      "\"${lib_path}\", "
    )
endforeach()

set(additional_forward_paths_install "")
foreach(lib_path ${SlicerExecutionModel_ADDITIONAL_FORWARD_PATHS_INSTALL})
  string(APPEND additional_forward_paths_install
      "\"${lib_path}\", "
    )
endforeach()

if(JsonCpp_FOUND)
  get_filename_component(json_library_dir ${JsonCpp_LIBRARY} DIRECTORY)
  set(json_forward_path_build
     "\"${json_library_dir}\" CONFIG_DIR_POST, \\
     \"${json_library_dir}\","
    )
  set(json_forward_path_install "\"${json_library_dir}\",")
endif()

set(GenerateCLP_FORWARD_DIR_BUILD "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(GenerateCLP_FORWARD_PATH_BUILD " \\
   \".\" CONFIG_DIR_POST, \\
   \"${ITK_DIR}/bin\" CONFIG_DIR_POST, \\
   ${json_forward_path_build} \\
   ${additional_forward_paths_build} \\
   \"${ModuleDescriptionParser_DIR}/bin\" CONFIG_DIR_POST \\
")
set(GenerateCLP_FORWARD_PATH_INSTALL " \\
   \".\", \\
   \"${ITK_DIR}/bin\", \\
   ${json_forward_path_install} \\
   ${additional_forward_paths_install} \\
   \"${ModuleDescriptionParser_DIR}/bin\" \\
")
set(GenerateCLP_FORWARD_EXE GenerateCLP)
set(GenerateCLP_FORWARD_DIR_INSTALL "..")


# --------------------------------------------------------------------------
# Build launcher
# --------------------------------------------------------------------------
set(GenerateCLPLauncher_SOURCE
  ${CMAKE_CURRENT_BINARY_DIR}/GenerateCLPLauncher.c
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GenerateCLPLauncher.c.in
  ${GenerateCLPLauncher_SOURCE}
  @ONLY)
add_custom_command(TARGET GenerateCLP POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E touch ${GenerateCLPLauncher_SOURCE}
  COMMENT "Force GenerateCLPLauncher to rebuild after GenerateCLP was modified"
  )
sem_add_executable(GenerateCLPLauncher
  ${GenerateCLPLauncher_SOURCE}
  )
list(APPEND targets_to_export GenerateCLPLauncher)
add_dependencies(GenerateCLPLauncher GenerateCLP)

# --------------------------------------------------------------------------
# Install launcher
# --------------------------------------------------------------------------
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(TARGETS GenerateCLPLauncher RUNTIME
    DESTINATION bin COMPONENT Development
    )
endif()

# --------------------------------------------------------------------------
# Export target
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_TARGETS ${targets_to_export} CACHE INTERNAL "${PROJECT_NAME} targets" FORCE)
export(TARGETS ${targets_to_export} FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)

# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include dirs" FORCE)

# --------------------------------------------------------------------------
# Configure support files
# --------------------------------------------------------------------------
configure_file(
  GenerateCLP.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/GenerateCLP.cmake
  @ONLY
  )
include(GenerateGenerateCLPConfig.cmake)

# --------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

