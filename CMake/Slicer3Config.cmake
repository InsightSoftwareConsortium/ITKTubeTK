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
#  We use Slicer3 variable names to simplify porting to Slicer3
set( Slicer3_INSTALL_BIN_DIR "bin" )
set( Slicer3_INSTALL_LIB_DIR "lib/TubeTK" )
set( Slicer3_INSTALL_INCLUDE_DIR "include/TubeTK" )
set( Slicer3_INSTALL_SHARE_DIR "share/TubeTK" )
set( Slicer3_INSTALL_ITKFACTORIES_DIR "${Slicer3_INSTALL_LIB_DIR}/ITKFactories" )

set( Slicer3_INSTALL_PLUGINS_BIN_DIR "${Slicer3_INSTALL_BIN_DIR}" )
set( Slicer3_INSTALL_PLUGINS_LIB_DIR "${Slicer3_INSTALL_LIB_DIR}" )
set( Slicer3_INSTALL_PLUGINS_INCLUDE_DIR "${Slicer3_INSTALL_INCLUDE_DIR}" )
set( Slicer3_INSTALL_PLUGINS_SHARE_DIR "${Slicer3_INSTALL_SHARE_DIR}" )
set( Slicer3_INSTALL_PLUGINS_CACHE_DIR "${Slicer3_INSTALL_LIB_DIR}/CACHE" )

set( Slicer3_HOME "${Slicer3_BINARY_DIR}" )

include( ${TubeTK_SOURCE_DIR}/CMake/Slicer3PluginsMacros.cmake )

# CLI_SOURCE_DIR is a Slicer3 variable.  Needs to point to the top of
#   the Slicer3 CLI modules.
set( CLI_SOURCE_DIR "${TubeTK_BINARY_DIR}/Slicer3" )
set( CLI_BINARY_DIR "${TubeTK_BINARY_DIR}/Slicer3-Build" )
make_directory( ${CLI_SOURCE_DIR} )
make_directory( ${CLI_BINARY_DIR} )
set( Slicer3_SOURCE_DIR ${CLI_SOURCE_DIR} )

# Define MACROs for downloading Slicer3 applications and libs used by TubeTK
if( NOT SVNCOMMAND )
  find_package( Subversion REQUIRED )
  set( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
endif( NOT SVNCOMMAND )

macro( svnGetSlicerCLI FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer3/trunk/Applications/CLI" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${CLI_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult
    )
  if( svnError )
    message( "Error while fetching CLI module from Slicer3 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${CLI_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerCLI FILENAME )

macro( svnGetSlicerBase FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer3/trunk/" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${CLI_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult
    )
  if( svnError )
    message( "Error while fetching module from Slicer3 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${CLI_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerBase FILENAME )
