##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################
set( DOCUMENTATION 
"TubeTK provides methods for extracting, registering, and quantifying the shape of tubular structures (e.g., vessels, nerves, bones, roads) in 2D and 3D images (e.g., from MRI, CT, ultrasound, or satelite imaging).")

itk_module( TubeTK
  ENABLE_SHARED
  DEPENDS
    ITKCommon
    ITKFFT
    ITKHDF5
    ITKIOImageBase
    ITKIOTransformBase
    ITKMetaIO
    ITKOptimizers
    ITKRegionGrowing
    ITKSmoothing
    ITKStatistics
    ITKSpatialObjects
    ITKTransform
    ITKThresholding
    MinimalPathExtraction
  COMPILE_DEPENDS
    ITKAnisotropicSmoothing
    ITKBinaryMathematicalMorphology
    ITKDistanceMap
    ITKImageFilterBase
    ITKImageFunction
    ITKImageFeature
    ITKImageGradient
    ITKImageIntensity
    ITKImageSources
    ITKIOMeta
    ITKIOSpatialObjects
    ITKIOCSV
    ITKLabelVoting
    ITKPDEDeformableRegistration
    ITKRegistrationCommon
    ITKTestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
)
