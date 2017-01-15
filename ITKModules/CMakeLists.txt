##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

##############################################################################
#
# Hint: To add an ITK module update CMake/Superbuild/ExternalProjectConfig.cmake
#
##############################################################################

# Including module 'ExternalProjectsConfig' set the
# following variables:
#
#   TubeTK_ITK_MODULES
#   <module>_URL
#   <module>_HASH_OR_TAG
#
include( ${TubeTK_SOURCE_DIR}/CMake/Superbuild/ExternalProjectsConfig.cmake )

# Include ITK. This ensures ITK_CMAKE_DIR and ITK_VERSION are set.
find_package( ITK REQUIRED )

# Disable "Factory Registration Manager".
#
# * Ensure TestDriver associated with ITKModules do NOT attempt to include
#     'itk(Image|Transform)IOFactoryRegisterManager.h'.
#
# * The driver does not rely on these files to register its IO factories.
#
# * Keeping the manager enabled is NOT supported by ITK build system and
#     causes link errors.
#
# * Simply setting ITK_NO_IO_FACTORY_REGISTER_MANAGER to 1 is not sufficient
#     because the compilation was possibly added in a parent directory.
#
get_property( _definitions DIRECTORY PROPERTY COMPILE_DEFINITIONS)
list(REMOVE_ITEM _definitions ITK_IO_FACTORY_REGISTER_MANAGER)
set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${_definitions}")

# Check if ITK support building external module.
set( ITK_WITH_EXTERNAL_MODULE_SUPPORT 1 )
set( _reason )
if( ITK_VERSION VERSION_LESS "4.10.0" )
  # Version of ITK with at least commit InsightSoftwareConsortium/ITK@a766d34
  # is required to support building of TubeTK.
  set( _reason
    "ITK version >= 4.10 is needed and current version is ${ITK_VERSION}" )
  set( ITK_WITH_EXTERNAL_MODULE_SUPPORT 0 )
endif()
if( NOT EXISTS ${ITK_CMAKE_DIR}/ITKModuleMacros.cmake )
  set( _reason
    "'ITKModuleMacros' module not found in [ITK_CMAKE_DIR:${ITK_CMAKE_DIR}]" )
  set( ITK_WITH_EXTERNAL_MODULE_SUPPORT 0 )
endif()

if(NOT ITK_WITH_EXTERNAL_MODULE_SUPPORT)
  message( STATUS
    "ITKModules: 'ITK External Module' not supported: ${_reason}" )
endif()

# XXX Ensure CMake modules provided by ITK are included first
#     when building using Slicer. This was motivated by the fact
#     ITK modules uses a modified version of GenerateExportHeader
#     not available in VTK.
set(CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} ${CMAKE_MODULE_PATH})

foreach( itk-module ${TubeTK_ITK_MODULES} )

  set( ${itk-module}_SOURCE_DIR ${${itk-module}_URL} )

  if( NOT EXISTS ${${itk-module}_URL} )

    message( STATUS "ITKModules: Downloading ${itk-module}" )

    set( ${itk-module}_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${itk-module} )

    # Create build directory
    set( _download_project_dir
      ${CMAKE_CURRENT_BINARY_DIR}/${itk-module}-download-project )
    file( MAKE_DIRECTORY ${_download_project_dir} )

    # Configure
    execute_process(
      COMMAND ${CMAKE_COMMAND}
        -G ${CMAKE_GENERATOR}
        -DPROJECT_NAME=${itk-module}
        -DGIT_REPOSITORY=${${itk-module}_URL}
        -DGIT_TAG=${${itk-module}_HASH_OR_TAG}
        -DMY_SOURCE_DIR=${${itk-module}_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/DownloadSource
      WORKING_DIRECTORY ${_download_project_dir}
      OUTPUT_QUIET
      RESULT_VARIABLE _result )
    if( NOT _result EQUAL 0 )
      message( FATAL_ERROR
        "Failed to configure '${itk-module}DownloadSource' project" )
    endif()

    # Build
    execute_process(
      COMMAND ${CMAKE_COMMAND}
        --build ${_download_project_dir}
        --config Release
      WORKING_DIRECTORY ${_download_project_dir}
      OUTPUT_QUIET
      RESULT_VARIABLE _result )
    if( NOT _result EQUAL 0 )
      message( FATAL_ERROR
        "Failed to build '${itk-module}DownloadSource' project" )
    endif()

  endif()

  set( ${itk-module}_BINARY_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/${itk-module}-build )

  if( ITK_WITH_EXTERNAL_MODULE_SUPPORT AND NOT TARGET ${itk-module} )
    message( STATUS "ITKModules: Adding '${itk-module}' directory" )

    # Note that adding the ITK module source directory will create
    # a target only if there is a "src" directory.
    add_subdirectory(
      ${${itk-module}_SOURCE_DIR}
      ${${itk-module}_BINARY_DIR} )

  endif()

  # Here we add a target only if needed. This could happen
  # in two scenarios:
  #
  # (1) The ITK module directory has been added above but it
  #     does NOT have a src folder.
  #     XXX Revisit this when ITK will use 'Modern CMake'. In that case an
  #       interface library will already be available.
  #
  # (2) ITK_DIR used to build TubeTK doesn't support External Module.
  #

  set( _target_type PUBLIC )
  if( NOT TARGET ${itk-module} )

    if( EXISTS ${${itk-module}_SOURCE_DIR}/src )

      find_package( ITK REQUIRED )
      include( ${ITK_USE_FILE} )

      # Note: Using GLOB_RECURSE is needed since we don't know
      #       the source files to add.

      file( GLOB_RECURSE _module_H_Files
        ${${itk-module}_SOURCE_DIR}/include/*.h )

      file( GLOB_RECURSE _module_HXX_Files
        ${${itk-module}_SOURCE_DIR}/include/*.hxx )

      file( GLOB_RECURSE _module_CXX_Files
        ${${itk-module}_SOURCE_DIR}/src/*.cxx )

      message( STATUS "ITKModules: Adding ${itk-module} STATIC library" )

      add_library( ${itk-module} STATIC
        ${_module_H_Files}
        ${_module_HXX_Files}
        ${_module_CXX_Files} )
      target_link_libraries( ${itk-module} ${ITK_LIBRARIES} )
      # Ensure TubeTKITK link against the required base libraries
      if( itk-module STREQUAL "TubeTKITK" )
        # Include depends file to set TubeTKITK_LIBRARIES in current scope
        include( ${itk-module}/tubetk-depends.cmake )
        target_link_libraries( ${itk-module} ${TubeTKITK_LIBRARIES} )
      endif()

      # Generate the export macro header for symbol visibility/Windows DLL
      #   declspec
      include( GenerateExportHeader )
      set( _export_header_file
        "${${itk-module}_BINARY_DIR}/include/${itk-module}Export.h" )
      # Note: The generate_export_header macro needs a relative path to
      #   generate the export header file. This path is created using the
      #   difference between the desired full path and the
      #   CMAKE_CURRENT_BINARY_DIR.
      file( RELATIVE_PATH _export_header_relative_path
        ${CMAKE_CURRENT_BINARY_DIR} ${_export_header_file} )
      generate_export_header( ${itk-module}
        EXPORT_FILE_NAME ${_export_header_relative_path}
        EXPORT_MACRO_NAME ${itk-module}_EXPORT
        # XXX Revisit this when `TEMPLATE_EXPORT_MACRO_NAME` feature will
        #   either:
        #     * be provided by 'GenerateExportHeader' module available in the
        #       minimum required CMake version.
        #     * or be available in an install tree of ITK and older version of
        #       ITK are not supported.
        #
        #TEMPLATE_EXPORT_MACRO_NAME ${itk-module}_TEMPLATE_EXPORT
        NO_EXPORT_MACRO_NAME ${itk-module}_HIDDEN
        STATIC_DEFINE ITK_STATIC )

    else()

      message( STATUS "ITKModules: Adding ${itk-module} INTERFACE library" )

      set( _target_type INTERFACE )
      add_library( ${itk-module} ${_target_type} )

    endif()

    TubeTKMacroInstallPlugins( ${itk-module} )

  endif()

  get_target_property( _imported ${itk-module} IMPORTED )
  if( NOT _imported )
    set( ${itk-module}_INC_DIRS ${${itk-module}_SOURCE_DIR}/include )
    if( ${itk-module}_ENABLE_SHARED )
      list( APPEND
        ${itk-module}_INC_DIRS ${${itk-module}_BINARY_DIR}/include )
    endif()
    target_include_directories( ${itk-module} ${_target_type}
      ${${itk-module}_INC_DIRS})
  endif()

endforeach()
