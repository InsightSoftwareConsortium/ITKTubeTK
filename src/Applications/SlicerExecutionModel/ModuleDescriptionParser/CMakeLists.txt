project(ModuleDescriptionParser LANGUAGES CXX)

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
# Option(s)
# --------------------------------------------------------------------------
if(NOT DEFINED BUILD_SHARED_LIBS)
  option(BUILD_SHARED_LIBS "Build with shared libraries." ON)
endif()

if(NOT DEFINED ModuleDescriptionParser_USE_SERIALIZER)
  option(ModuleDescriptionParser_USE_SERIALIZER "Serialization and deserialization support is added for ModuleDescriptionParser objects." OFF)
  mark_as_advanced(ModuleDescriptionParser_USE_SERIALIZER)
endif()
if(ModuleDescriptionParser_USE_SERIALIZER)
  find_package(ParameterSerializer REQUIRED)
  include_directories(${ParameterSerializer_INCLUDE_DIRS})
endif()

# --------------------------------------------------------------------------
# Prerequisites
# --------------------------------------------------------------------------

#
# ITK
#
set(${PROJECT_NAME}_ITK_COMPONENTS
  ITKCommon# For itksys
  ITKIOXML # For ITKEXPAT
  ITKExpat # For Expat library
  )

find_package(ITK 4.3 COMPONENTS ${${PROJECT_NAME}_ITK_COMPONENTS} REQUIRED)
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
include(CTest)

#-----------------------------------------------------------------------------
# Output directories.
#-----------------------------------------------------------------------------
foreach(type LIBRARY RUNTIME ARCHIVE)
  # Make sure the directory exists
  if(DEFINED ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY
     AND NOT EXISTS ${ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY})
    message(FATAL_ERROR "ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY is set to a non-existing directory [${ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY}]")
  endif()

  if(NOT DEFINED ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY)
    set(ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY ${ModuleDescriptionParser_BINARY_DIR}/bin)
  endif()
  set(CMAKE_${type}_OUTPUT_DIRECTORY ${ModuleDescriptionParser_CMAKE_${type}_OUTPUT_DIRECTORY} CACHE INTERNAL "Single output directory for building all libraries.")
endforeach()

# --------------------------------------------------------------------------
# Include dirs
# --------------------------------------------------------------------------
set(include_dirs
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

include_directories(${include_dirs})

# --------------------------------------------------------------------------
# Configure headers
# --------------------------------------------------------------------------
set(configure_header_file ModuleDescriptionParserConfigure.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
  )

# --------------------------------------------------------------------------
# Install headers
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  set(${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT ON)
endif()
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES
    BatchMakeUtilities.h
    itkSEMModuleDescriptionSerializer.h
    itkSEMModuleParameterGroupSerializer.h
    itkSEMModuleParameterSerializer.h
    JsonSerializationUtilities.h
    ModuleDescription.h
    ModuleDescriptionParser.h
    ModuleDescriptionUtilities.h
    ModuleLogo.h
    ModuleParameterGroup.h
    ModuleParameter.h
    ModuleProcessInformation.h
    ${CMAKE_CURRENT_BINARY_DIR}/ModuleDescriptionParserExport.h
    ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
    DESTINATION include/${PROJECT_NAME}
    COMPONENT Development
    )
endif()

# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------
set(ModuleDescriptionParser_SRCS
  ModuleParameter.cxx
  ModuleParameterGroup.cxx
  ModuleDescriptionUtilities.cxx
  ModuleDescription.cxx
  ModuleDescriptionParser.cxx
  ModuleProcessInformation.cxx
  ModuleLogo.cxx
  BatchMakeUtilities.cxx
  )
if(ModuleDescriptionParser_USE_SERIALIZER)
  list(APPEND ModuleDescriptionParser_SRCS
    itkSEMModuleDescriptionSerializer.cxx
    itkSEMModuleParameterSerializer.cxx
    itkSEMModuleParameterGroupSerializer.cxx
    )
endif()

# --------------------------------------------------------------------------
# Author warnings for unavailable features
# --------------------------------------------------------------------------
if(DEFINED USE_BFD)
  message(AUTHOR_WARNING
    "Option USE_BFD has no effect. BinaryFileDescriptor support has been removed.")
endif()
if(DEFINED ModuleDescriptionParser_USE_PYTHON)
  message(AUTHOR_WARNING
    "Option ModuleDescriptionParser_USE_PYTHON has no effect. Python support has been removed.")
endif()

# --------------------------------------------------------------------------
# Build library
# --------------------------------------------------------------------------
set(lib_name ModuleDescriptionParser)
include(GenerateExportHeader)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/ModuleDescriptionParserPragmas.h.in ModuleDescriptionParser_PRAGMAS)
add_library(${lib_name} ${ModuleDescriptionParser_SRCS})
generate_export_header(${lib_name}
  BASE_NAME ${lib_name}
  EXPORT_MACRO_NAME ${lib_name}_EXPORT
  EXPORT_FILE_NAME ${lib_name}Export.h
  STATIC_DEFINE ${lib_name}_STATIC
  CUSTOM_CONTENT_FROM_VARIABLE ModuleDescriptionParser_PRAGMAS
  )

## Always build an explicitly static library for linking against GenerateCLP so that
## GenerateCLP can be run from CMake without having to set DYLD_LIBRARY_PATH or LD_LIBRARY_PATH
## to the future location of of libModuleDescriptionParser.so
## add_library(${lib_name}-static STATIC
  ## ${ModuleDescriptionParser_SRCS}
  ## )

set(libs
  ${ITK_LIBRARIES}
  )

#
# Append extra platform dependent libraries required for linking
#

if(ModuleDescriptionParser_USE_SERIALIZER)
  list(APPEND libs ${ParameterSerializer_LIBRARIES})
endif()

target_link_libraries(${lib_name} ${libs})
## target_link_libraries(${lib_name}-static ${link_libs})
set(USE_COMPILER_HIDDEN_VISIBILITY ON)
if (USE_COMPILER_HIDDEN_VISIBILITY) 
  # Prefer to use target properties supported by newer cmake 
  set_target_properties(${lib_name} PROPERTIES CXX_VISIBILITY_PRESET hidden) 
  set_target_properties(${lib_name} PROPERTIES C_VISIBILITY_PRESET hidden) 
  set_target_properties(${lib_name} PROPERTIES VISIBILITY_INLINES_HIDDEN 1) 
endif()

option(USE_COMPILER_HIDDEN_VISIBILITY "Use HIDDEN visibility support if available." ON)
mark_as_advanced(USE_COMPILER_HIDDEN_VISIBILITY)
if(USE_COMPILER_HIDDEN_VISIBILITY)
  set_target_properties(${lib_name} PROPERTIES CXX_VISIBILITY_PRESET hidden)
  set_target_properties(${lib_name} PROPERTIES C_VISIBILITY_PRESET hidden)
  set_target_properties(${lib_name} PROPERTIES VISIBILITY_INLINES_HIDDEN 1)
endif()

#
# Apply user-defined properties to the library target.
#
if(ModuleDescriptionParser_LIBRARY_PROPERTIES)
  set_target_properties(${lib_name} PROPERTIES ${ModuleDescriptionParser_LIBRARY_PROPERTIES} )
endif()

# --------------------------------------------------------------------------
# Export target
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_TARGETS ${lib_name} CACHE INTERNAL "${PROJECT_NAME} targets" FORCE)
export(TARGETS ${lib_name} FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)

# --------------------------------------------------------------------------
# Install library
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_BIN_DIR)
  set(${PROJECT_NAME}_INSTALL_BIN_DIR bin)
endif()

if(NOT DEFINED ${PROJECT_NAME}_INSTALL_LIB_DIR)
  set(${PROJECT_NAME}_INSTALL_LIB_DIR lib/${PROJECT_NAME})
endif()

install(TARGETS ${lib_name}
  RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_BIN_DIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${${PROJECT_NAME}_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${${PROJECT_NAME}_INSTALL_LIB_DIR} COMPONENT Development
  )

# --------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

# --------------------------------------------------------------------------
# Install support files
# --------------------------------------------------------------------------
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/UseModuleDescriptionParser.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/install/ModuleDescriptionParserConfig.cmake
    DESTINATION lib/${lib_name} COMPONENT Development)
endif()

# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include dirs" FORCE)

# --------------------------------------------------------------------------
# Configure support files
# --------------------------------------------------------------------------
include(GenerateModuleDescriptionParserConfig.cmake)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/UseModuleDescriptionParser.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/UseModuleDescriptionParser.cmake COPYONLY)

