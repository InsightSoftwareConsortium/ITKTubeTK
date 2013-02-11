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
#
#  We use Slicer4 variable names to simplify porting to Slicer4
set( Slicer4_INSTALL_BIN_DIR "bin" )
set( Slicer4_INSTALL_LIB_DIR "lib/TubeTK" )
set( Slicer4_INSTALL_INCLUDE_DIR "include/TubeTK" )
set( Slicer4_INSTALL_SHARE_DIR "share/TubeTK" )
set( Slicer4_INSTALL_ITKFACTORIES_DIR "${Slicer4_INSTALL_LIB_DIR}/ITKFactories" )

set( Slicer4_INSTALL_PLUGINS_BIN_DIR "${Slicer4_INSTALL_BIN_DIR}" )
set( Slicer4_INSTALL_PLUGINS_LIB_DIR "${Slicer4_INSTALL_LIB_DIR}" )
set( Slicer4_INSTALL_PLUGINS_INCLUDE_DIR "${Slicer4_INSTALL_INCLUDE_DIR}" )
set( Slicer4_INSTALL_PLUGINS_SHARE_DIR "${Slicer4_INSTALL_SHARE_DIR}" )
set( Slicer4_INSTALL_PLUGINS_CACHE_DIR "${Slicer4_INSTALL_LIB_DIR}/CACHE" )

set( Slicer4_HOME "${Slicer4_BINARY_DIR}" )

include( ${TubeTK_SOURCE_DIR}/CMake/Slicer4PluginsMacros.cmake )
include( ${SlicerExecutionModel_DIR}/CMake/SEMMacroBuildCLI.cmake )

set( Slicer_SOURCE_DIR "${TubeTK_BINARY_DIR}/Slicer4" )
set( Slicer_BINARY_DIR "${TubeTK_BINARY_DIR}/Slicer4-Build" )
make_directory( ${Slicer_SOURCE_DIR} )
make_directory( ${Slicer_BINARY_DIR} )

# Needed by SlicerBaseCLI
set( Slicer_INSTALL_BIN_DIR ${Slicer4_INSTALL_BIN_DIR} )
set( Slicer_INSTALL_LIB_DIR ${Slicer4_INSTALL_LIB_DIR} )

# Define MACROs for downloading Slicer4 applications and libs used by TubeTK
if( NOT SVNCOMMAND )
  find_package( Subversion REQUIRED )
  set( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
endif( NOT SVNCOMMAND )

macro( svnGetSlicerCLI FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer4/trunk/Modules/CLI/" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${Slicer_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult
    )
  if( svnError )
    message( "Error while fetching CLI module from Slicer4 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${Slicer_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( "Error = ${svnError}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerCLI FILENAME )

macro( svnGetSlicerBase FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer4/trunk/" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${Slicer_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult
    )
  if( svnError )
    message( "Error while fetching module from Slicer4 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${Slicer_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerBase FILENAME )
