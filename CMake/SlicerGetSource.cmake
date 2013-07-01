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

set( TubeTK_Slicer_SOURCE_DIR "${TubeTK_BINARY_DIR}/Slicer" )
set( TubeTK_Slicer_BINARY_DIR "${TubeTK_BINARY_DIR}/Slicer-build" )
make_directory( ${TubeTK_Slicer_SOURCE_DIR} )
make_directory( ${TubeTK_Slicer_BINARY_DIR} )

# Define macros for downloading Slicer4 applications and libs used by TubeTK
if( NOT SVNCOMMAND )
  find_package( Subversion REQUIRED )
  set( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
endif( NOT SVNCOMMAND )

macro( svnGetSlicerCLI FILENAME )
  svnGetSlicerBase( "Modules/CLI/${FILENAME}" )
endmacro( svnGetSlicerCLI FILENAME )

macro( svnGetSlicerBase FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer4/trunk/" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${TubeTK_Slicer_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult )
  if( svnError )
    message( "Error while fetching module from Slicer4 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${TubeTK_Slicer_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( svnGetSlicerBase FILENAME )
