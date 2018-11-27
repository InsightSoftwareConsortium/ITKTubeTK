##############################################################################
#
# Library:   TubeTKLib
#
# Copyright Kitware Inc.
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

project( TubeTKLibMetaIO )

set( TubeTKLib_Base_MetaIO_H_Files
  itktubeMetaLDA.h
  itktubeMetaNJetLDA.h
  itktubeMetaClassPDF.h
  itktubeMetaRidgeSeed.h
  itktubeMetaTubeExtractor.h )

set( TubeTKLib_Base_MetaIO_HXX_Files )

set( TubeTKLib_Base_MetaIO_SRCS
  itktubeMetaClassPDF.cxx
  itktubeMetaLDA.cxx
  itktubeMetaNJetLDA.cxx
  itktubeMetaRidgeSeed.cxx
  itktubeMetaTubeExtractor.cxx )

add_library( ${PROJECT_NAME} STATIC
  ${TubeTKLib_Base_MetaIO_H_Files}
  ${TubeTKLib_Base_MetaIO_HXX_Files}
  ${TubeTKLib_Base_MetaIO_SRCS} )

target_link_libraries( ${PROJECT_NAME} PUBLIC
  TubeTKLibCommon )

target_include_directories( ${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

if( TubeTKLib_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTKLib_BUILD_TESTING )

if( TubeTKLib_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTKLib_Base_Numerics_H_Files}
    ${TubeTKLib_Base_Numerics_HXX_Files}
    DESTINATION ${TubeTKLib_INSTALL_INCLUDE_DIR}
    COMPONENT Development )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS ${PROJECT_NAME} )
