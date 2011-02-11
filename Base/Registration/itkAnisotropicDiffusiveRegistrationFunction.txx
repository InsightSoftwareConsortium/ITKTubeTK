/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itkAnisotropicDiffusiveRegistrationFunction_txx
#define __itkAnisotropicDiffusiveRegistrationFunction_txx

#include "itkAnisotropicDiffusiveRegistrationFunction.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
AnisotropicDiffusiveRegistrationFunction
 < TFixedImage, TMovingImage, TDeformationField >
::AnisotropicDiffusiveRegistrationFunction()
{
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}

/**
 * Called at the beginning of each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeIteration()
{
  Superclass::InitializeIteration();
}

/**
  * Computes the update term
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::PixelType
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeUpdate(
    const NeighborhoodType &neighborhood,
    const NormalVectorNeighborhoodType
        &normalVectorNeighborhood,
    const DiffusionTensorNeighborhoodType
        &tangentialTensorNeighborhood,
    const TensorDerivativeImageRegionType
        &tangentialTensorDerivativeRegion,
    const DeformationVectorComponentNeighborhoodArrayType
        &tangentialDeformationComponentNeighborhoods,
    const DiffusionTensorNeighborhoodType
        &normalTensorNeighborhood,
    const TensorDerivativeImageRegionType
        &normalTensorDerivativeRegion,
    const DeformationVectorComponentNeighborhoodArrayType
        &normalDeformationComponentNeighborhoods,
    const SpacingType &spacing,
    void *globalData,
    const FloatOffsetType &offset )
{
  return Superclass::ComputeUpdate( neighborhood,
                                    normalVectorNeighborhood,
                                    tangentialTensorNeighborhood,
                                    tangentialTensorDerivativeRegion,
                                    tangentialDeformationComponentNeighborhoods,
                                    normalTensorNeighborhood,
                                    normalTensorDerivativeRegion,
                                    normalDeformationComponentNeighborhoods,
                                    spacing,
                                    globalData,
                                    offset );
}

} // end namespace itk

#endif
