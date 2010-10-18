##############################################################################
# 
# Library:   TubeTK/VTree
# 
# Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt
# 
# Original implementation:
# Copyright University of North Carolina, Chapel Hill, NC, USA.
# 
# Revised implementation:
# Copyright Kitware Inc., Carrboro, NC, USA.
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

PROJECT( VTree )

FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

SET( VTree_SRCS
  itkOptBrent1D.cxx
  itkOptGoldenMean1D.cxx
  itkOptParabolicFit1D.cxx
  itkOptimizer1D.cxx
  itkOptimizerND.cxx
  itkSpline1D.cxx
  itkSplineApproximation1D.cxx
  itkSplineND.cxx
  )

ADD_LIBRARY( VTree STATIC ${VTree_SRCS} )
TARGET_LINK_LIBRARIES( VTree ${ITK_LIBRARIES} )

slicer3_set_plugins_output_path( VTree )
set( TARGETS VTree )
slicer3_install_plugins( ${TARGETS} )

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )


