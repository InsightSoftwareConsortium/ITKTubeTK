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

project( TubeFiltering )

FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

set( TubeTK_Base_Numerics_H_Files
  itkNJetImageFunction.h
  itkJointHistogramImageFunction.h
  tubeMatrixMath.h
  tubeOptParabolicFit1D.h
  tubeSpline1D.h
  tubeOptBrent1D.h
  tubeSplineApproximation1D.h
  tubeOptGoldenMean1D.h
  tubeSplineND.h
  tubeOptimizer1D.h
  tubeTubeMath.h
  tubeOptimizerND.h
  tubeUserFunc.h
  )

set( TubeTK_Base_Numerics_TXX_Files
  itkNJetImageFunction.txx
  itkJointHistogramImageFunction.txx
  tubeTubeMath.txx
  tubeMatrixMath.txx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  )

SET( TubeTK_Base_Numerics_SRCS
  itkTubeNJetLDAGenerator2F.cxx
  tubeOptBrent1D.cxx
  tubeOptGoldenMean1D.cxx
  tubeOptParabolicFit1D.cxx
  tubeOptimizer1D.cxx
  tubeOptimizerND.cxx
  tubeSpline1D.cxx
  tubeSplineApproximation1D.cxx
  tubeSplineND.cxx
  )

ADD_LIBRARY( TubeNumerics STATIC ${TubeTK_Base_Numerics_SRCS} )
TARGET_LINK_LIBRARIES( TubeNumerics ${ITK_LIBRARIES} )

if( NOT TubeTK_USE_ITKV4 )
slicer3_set_plugins_output_path( TubeNumerics )
set( TARGETS TubeNumerics )
slicer3_install_plugins( ${TARGETS} )
endif( NOT TubeTK_USE_ITKV4 )

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )

#install(FILES 
  #${TubeTK_Base_Numerics_H_Files}
  #${TubeTK_Base_Numerics_TXX_Files}
  #DESTINATION include)

