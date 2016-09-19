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

project( TubeTKMetaIO )

set( TubeTK_Base_MetaIO_H_Files
  itktubeMetaLDA.h
  itktubeMetaNJetLDA.h
  itktubeMetaClassPDF.h
  itktubeMetaRidgeSeed.h
  itktubeMetaTubeExtractor.h )

set( TubeTK_Base_MetaIO_HXX_Files )

set( TubeTK_Base_MetaIO_SRCS
  itktubeMetaClassPDF.cxx
  itktubeMetaLDA.cxx
  itktubeMetaNJetLDA.cxx
  itktubeMetaRidgeSeed.cxx
  itktubeMetaTubeExtractor.cxx )

add_library( ${PROJECT_NAME} STATIC
  ${TubeTK_Base_MetaIO_H_Files}
  ${TubeTK_Base_MetaIO_HXX_Files}
  ${TubeTK_Base_MetaIO_SRCS} )

target_link_libraries( ${PROJECT_NAME} PUBLIC
  TubeTKCommon )

target_include_directories( ${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

if( TubeTK_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTK_BUILD_TESTING )

if( TubeTK_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTK_Base_Numerics_H_Files}
    ${TubeTK_Base_Numerics_HXX_Files}
    DESTINATION ${TubeTK_INSTALL_INCLUDE_DIR}
    COMPONENT Development )

  TubeTKMacroInstallPlugins( ${PROJECT_NAME} )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS ${PROJECT_NAME} )
