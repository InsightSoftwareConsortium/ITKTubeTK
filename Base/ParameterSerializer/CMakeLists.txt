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

project( ParameterSerializer )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

find_package( JsonCpp REQUIRED )
include_directories( ${JsonCpp_INCLUDE_DIRS} )

set( TubeTK_Base_ParameterSerializer_H_Files
  USTK/itkAcousticImpulseResponseImageFilterSerializer.h
  )

set( TubeTK_Base_ParameterSerializer_TXX_Files
  USTK/itkAcousticImpulseResponseImageFilterSerializer.txx
  )

set( TubeTK_Base_ParameterSerializer_SRCS
  )

#add_library( TubeParameterSerializer STATIC
  #${TubeTK_Base_ParameterSerializer_SRCS} )
#target_link_libraries( TubeParameterSerializer
  #${ITK_LIBRARIES}
  #${JsonCpp_LIBRARIES}
  #)

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_ParameterSerializer_H_Files}
  ${TubeTK_Base_ParameterSerializer_TXX_Files}
  DESTINATION include
  )
