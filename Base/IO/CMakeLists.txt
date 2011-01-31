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

project( TubeIO )

FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

set( TubeTK_Base_IO_H_Files
  itkTubeMetaLDA.h
  itkTubeMetaNJetLDA.h
  )

set( TubeTK_Base_IO_TXX_Files
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  )

SET( TubeTK_Base_IO_SRCS
  itkTubeMetaLDA.cxx
  itkTubeMetaNJetLDA.cxx
  )

ADD_LIBRARY( TubeIO STATIC ${TubeTK_Base_IO_SRCS} )
TARGET_LINK_LIBRARIES( TubeIO ${ITK_LIBRARIES} )

slicer3_set_plugins_output_path( TubeIO )
set( TARGETS TubeIO )
slicer3_install_plugins( ${TARGETS} )

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )

#install(FILES 
  #${TubeTK_Base_Numerics_H_Files}
  #${TubeTK_Base_Numerics_TXX_Files}
  #DESTINATION include)

