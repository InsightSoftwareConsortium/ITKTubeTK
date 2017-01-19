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
# Record ITK, VTK, and GenerateCLP libraries used by TubeTK for
# packaging.
file( GLOB_RECURSE ITK_BINARIES "${ITK_DIR}/lib/lib*" )
if( TubeTK_USE_VTK )
  file( GLOB_RECURSE VTK_BINARIES "${VTK_DIR}/bin/lib*" )
endif( TubeTK_USE_VTK )
file( GLOB_RECURSE SEM_BINARIES "${GenerateCLP_DIR}/bin/lib*" )
file( GLOB_RECURSE JsonCpp_BINARIES "${JsonCpp_DIR}/lib/lib*" )
list( APPEND TubeTK_SUPPORT_BINARIES
  ${SEM_BINARIES} ${VTK_BINARIES} ${ITK_BINARIES} )
if( TubeTK_USE_JsonCpp )
  list( APPEND TubeTK_SUPPORT_BINARIES ${JsonCpp_BINARIES} )
endif( TubeTK_USE_JsonCpp )
install( FILES ${TubeTK_SUPPORT_BINARIES}
  DESTINATION ${TubeTK_INSTALL_LIB_DIR}
  COMPONENT Runtime)

# Record the environment setup files for packaging.
set( TubeTK_LINUX_SETUP
  ${TubeTK_SOURCE_DIR}/CMake/TubeTKSourceMeToSetupLinux.sh )
set( TubeTK_WINDOWS_SETUP
  ${TubeTK_SOURCE_DIR}/CMake/TubeTKRunMeToSetupWindows.bat )

include( InstallRequiredSystemLibraries )

set( CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "TubeTK for the analysis of geometric structures in images." )
set( CPACK_PACKAGE_VENDOR "Kitware, Inc." )
set( CPACK_PACKAGE_DESCRIPTION_FILE
  "${TubeTK_SOURCE_DIR}/CMake/CPack/DESCRIPTION.txt" )
set( CPACK_RESOURCE_FILE_LICENSE
  "${TubeTK_SOURCE_DIR}/CMake/CPack/LICENSE.txt" )

set( CPACK_PACKAGE_VERSION_MAJOR ${TubeTK_VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${TubeTK_VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH ${TubeTK_VERSION_PATCH} )
set( CPACK_PACKAGE_VERSION ${TubeTK_VERSION_BUILD} )

set( CPACK_PACKAGE_CONTACT "stephen.aylward@kitware.com" )

set( CPACK_PACKAGE_INSTALL_DIRECTORY
  "TubeTK_${TubeTK_VERSION_MAJOR}.${TubeTK_VERSION_MINOR}" )

if( WIN32 AND NOT UNIX )
  # There is a bug in NSI that does not handle full UNIX paths properly. Make
  # sure there is at least one set of four backslashes.
  set( CPACK_PACKAGE_ICON "${TubeTK_SOURCE_DIR}\\\\Icon.png" )
  set( CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\TubeTK.exe" )
  set( CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} TubeTK" )
  set( CPACK_NSIS_HELP_LINK "http:\\\\\\\\tubetk.org" )
  set( CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\tubetk.org" )
  set( CPACK_NSIS_CONTACT "stephen.aylward@kitware.com" )
  set( CPACK_NSIS_MODIFY_PATH ON )
else( WIN32 AND NOT UNIX )
  set( CPACK_STRIP_FILES true )
  set( CPACK_SOURCE_STRIP_FILES true )
endif( WIN32 AND NOT UNIX )

set( CMAKE_INSTALL_RPATH "$ORIGIN/../..;." )

# Standard CPack setup.
option( CMAKE_INSTALL_DEBUG_LIBRARIES
  "Install the Microsoft runtime debug libraries with CMake." FALSE )
mark_as_advanced( CMAKE_INSTALL_DEBUG_LIBRARIES )

# By default, do not warn when built on machines using only Microsoft Visual
# Studio Express.
if( NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS )
  set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON )
endif( NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS )

# Make this explicit here, rather than accepting the CPack default value, so
# we can refer to it.
set( CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}" )

# Installers for 32- vs. 64-bit CMake:
#   - Root install directory (displayed to end user at installer runtime).
#   - NSIS package/display name (text used in the installer GUI).
#   - Registry key used to store information about the installation.
if( CMAKE_CL_64 )
  set( CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64" )
  set( CPACK_NSIS_PACKAGE_NAME
    "${CPACK_PACKAGE_INSTALL_DIRECTORY} ( Win64 )" )
  set( CPACK_PACKAGE_INSTALL_REGISTRY_KEY
    "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION} ( Win64 )" )
else( CMAKE_CL_64 )
  set( CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES" )
  set( CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}" )
  set( CPACK_PACKAGE_INSTALL_REGISTRY_KEY
    "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}" )
endif( CMAKE_CL_64 )

if( NOT DEFINED CPACK_SYSTEM_NAME )
  set( CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} )
endif( NOT DEFINED CPACK_SYSTEM_NAME )

if( ${CPACK_SYSTEM_NAME} MATCHES Windows )
  if( CMAKE_CL_64 )
    set( CPACK_SYSTEM_NAME win64-x64 )
  else( CMAKE_CL_64 )
    set( CPACK_SYSTEM_NAME win32-x86 )
  endif( CMAKE_CL_64 )
endif( ${CPACK_SYSTEM_NAME} MATCHES Windows )

set( CPACK_SOURCE_PACKAGE_FILE_NAME "TubeTK-${TubeTK_VERSION_BUILD}" )
set( CPACK_PACKAGE_FILE_NAME
  "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}" )

# Include CPack model once all variables are set.
include( CPack )
