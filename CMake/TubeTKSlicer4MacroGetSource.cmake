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

# Define macros for downloading Slicer4 applications and libs used by TubeTK
if( NOT SVNCOMMAND )
  find_package( Subversion REQUIRED )
  set( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
endif( NOT SVNCOMMAND )

macro( TubeTKSlicer4MacroGetSource FILENAME )
  set( svnCmd "co" )
  set( svnSite "http://svn.slicer.org/Slicer4/trunk/" )
  execute_process(
    COMMAND "${SVNCOMMAND}" ${svnCmd} ${svnSite}/${FILENAME} ${TubeTK_SLICER_SOURCE_DIR}/${FILENAME}
    OUTPUT_VARIABLE svnOutput
    ERROR_VARIABLE svnError
    RESULT_VARIABLE svnResult )
  if( svnError )
    message( "Error while fetching module from Slicer4 SVN..." )
    message( "Svn = ${SVNCOMMAND}" )
    message( "Command = ${svnCmd}" )
    message( "Site = ${svnSite}" )
    message( "File = ${FILENAME}" )
    message( "Dir = ${TubeTK_SLICER_SOURCE_DIR}" )
    message( "Output = ${svnOutput}" )
    message( FATAL "Result = ${svnResult}" )
  endif( svnError )
endmacro( TubeTKSlicer4MacroGetSource FILENAME )

macro( TubeTKSlicer4MacroGetCLISource FILENAME )
  TubeTKSlicer4MacroGetSource( "Modules/CLI/${FILENAME}" )
endmacro( TubeTKSlicer4MacroGetCLISource FILENAME )
