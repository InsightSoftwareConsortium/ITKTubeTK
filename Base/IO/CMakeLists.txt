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

project( TubeTKIO )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_IO_H_Files
  itktubeMetaLDA.h
  itktubeMetaNJetLDA.h
  itktubeMetaClassPDF.h
  itktubeMetaRidgeSeed.h
  itktubeMetaTubeParams.h
  itktubeRidgeSeedFilterIO.h )

set( TubeTK_Base_IO_HXX_Files
  itktubeRidgeSeedFilterIO.hxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Common )

set( TubeTK_Base_IO_SRCS
  itktubeMetaLDA.cxx
  itktubeMetaNJetLDA.cxx
  itktubeMetaRidgeSeed.cxx
  itktubeMetaClassPDF.cxx
  itktubeMetaTubeParams.cxx )

add_library( TubeTKIO STATIC
  ${TubeTK_Base_IO_H_Files}
  ${TubeTK_Base_IO_HXX_Files}
  ${TubeTK_Base_IO_SRCS} )

target_link_libraries( TubeTKIO ${ITK_LIBRARIES} )

set( TARGETS TubeTKIO )
TubeTKMacroInstallPlugins( ${TARGETS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Numerics_H_Files}
  ${TubeTK_Base_Numerics_HXX_Files}
  DESTINATION include )

# Export target
set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKIO )
