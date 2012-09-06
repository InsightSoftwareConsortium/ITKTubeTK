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

project(Registration)

set( TubeTK_Base_Registration_H_Files
  itkAnisotropicDiffusiveRegistrationFilter.h
  itkAnisotropicDiffusiveRegistrationFunction.h
  itkAnisotropicDiffusiveSparseRegistrationFilter.h
  itkDiffusiveRegistrationFilter.h
  itkDiffusiveRegistrationFilterUtils.h
  itkImageToTubeRigidMetric2.h
  itkImageToTubeRigidMetric.h
  itkImageToTubeRigidRegistration2.h
  itkImageToTubeRigidRegistration.h
  itkMeanSquareRegistrationFunction.h
  itkTubeToTubeTransformFilter.h
  )

set( TubeTK_Base_Registration_TXX_Files
  itkAnisotropicDiffusiveRegistrationFilter.txx
  itkAnisotropicDiffusiveRegistrationFunction.txx
  itkAnisotropicDiffusiveSparseRegistrationFilter.txx
  itkDiffusiveRegistrationFilter.txx
  itkDiffusiveRegistrationFilterUtils.txx
  itkImageToTubeRigidMetric2.txx
  itkImageToTubeRigidMetric.txx
  itkImageToTubeRigidRegistration2.txx
  itkImageToTubeRigidRegistration.txx
  itkMeanSquareRegistrationFunction.txx
  itkTubeToTubeTransformFilter.txx
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install(FILES
  ${TubeTK_Base_Segmentation_H_Files}
  ${TubeTK_Base_Segmentation_TXX_Files}
  DESTINATION include
  )
