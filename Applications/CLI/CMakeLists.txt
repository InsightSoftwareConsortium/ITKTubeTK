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

add_subdirectory( SampleCLIApplication )
add_subdirectory( HessianTubeness2D )
add_subdirectory( ExtractCurves2D )
add_subdirectory( SimulateTomosynthesis )
add_subdirectory( ComputeImageSimilarity )
add_subdirectory( Crop )
add_subdirectory( Merge )
add_subdirectory( MatchImageWithPrior )
add_subdirectory( CompareImageWithPrior )
add_subdirectory( RandomlyModifyBinaryMask )
add_subdirectory( DeblendImagesUsingPrior )
add_subdirectory( OtsuThreshold )

# The rest of the CMakeLists.txt file is dedicated to downloading, configuring,
#  building, and installing CLI modules that are distributed with Slicer3 and
#  that we would like to access from TubeTK.
IF( NOT SVNCOMMAND )
  FIND_PACKAGE( Subversion REQUIRED )
  SET( SVNCOMMAND ${Subversion_SVN_EXECUTABLE} )
ENDIF( NOT SVNCOMMAND )

# CLI_SOURCE_DIR is a Slicer3 variable.  Needs to point to the top of
#   the Slicer3 CLI modules.
set( CLI_SOURCE_DIR "${TubeTK_BINARY_DIR}/Applications/CLI/Slicer3-Src" )
set( CLI_BINARY_DIR "${TubeTK_BINARY_DIR}/Applications/CLI/Slicer3-Bin" )
make_directory( ${CLI_SOURCE_DIR} )
make_directory( ${CLI_BINARY_DIR} )
set( Slicer3_SOURCE_DIR ${CLI_SOURCE_DIR} )

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

svnGetSlicerBase( Resources ${CLI_SOURCE_DIR} )
svnGetSlicerCLI( Templates ${CLI_SOURCE_DIR} )
svnGetSlicerCLI( RegisterImagesModule ${CLI_SOURCE_DIR} )

set( CXX_TEST_PATH  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR} )
set( LAUNCH_EXE  "" )

set( SUPPORT_2D_IMAGES 1 )
add_subdirectory( ${CLI_SOURCE_DIR}/RegisterImagesModule
                  ${CLI_BINARY_DIR}/RegisterImagesModule )

