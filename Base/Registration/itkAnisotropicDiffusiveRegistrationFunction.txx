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
  m_UseAnisotropicRegularization = true;
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

  os << indent << "Use anisotropic regularization: "
     << ( m_UseAnisotropicRegularization ? "on" : "off" ) << std::endl;
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
    const DiffusionTensorNeighborhoodType
        &tangentialTensorNeighborhood,
    const TensorDerivativeImageRegionType
        &tangentialTensorDerivativeRegion,
    const DeformationVectorComponentNeighborhoodArrayType
        &tangentialDeformationComponentNeighborhoods,
    const NormalVectorNeighborhoodType
        &normalVectorNeighborhood,
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
  // Get the global data structure
  GlobalDataStruct * gd = ( GlobalDataStruct * ) globalData;

  // The superclass computes the intensity term and regularization in the
  // tangential plane ( i.e. div(P^P \grad(u_l))(e_l) )
  PixelType intensityTermPlusTangentialRegularizationTerm
      = Superclass::ComputeUpdate(neighborhood,
                                  tangentialTensorNeighborhood,
                                  tangentialTensorDerivativeRegion,
                                  tangentialDeformationComponentNeighborhoods,
                                  spacing,
                                  globalData,
                                  offset );

  // Compute the normal component of the regularization update term
  PixelType normalRegularizationTerm;
  normalRegularizationTerm.Fill(0);
  if( this->GetComputeRegularizationTerm()
    && this->GetUseAnisotropicRegularization() )
    {
    NormalVectorType                  normalVector;
    DeformationVectorComponentType    intermediateNormalRegularizationComponent;
    PixelType                         intermediateNormalRegularizationTerm;
    NormalVectorType                  nln; // n(l)n

    // Get the normal at this pixel once
    const typename FixedImageType::IndexType index = neighborhood.GetIndex();
    normalVector
        = normalVectorNeighborhood.GetImagePointer()->GetPixel( index );

    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      // Compute the regularization in the normal direction
      // Compute div(w^2nn^T grad(u_l^\perp))
      intermediateNormalRegularizationComponent
          = this->GetRegularizationFunctionPointer()->ComputeUpdate(
              normalDeformationComponentNeighborhoods[i],
              normalTensorNeighborhood,
              normalTensorDerivativeRegion,
              spacing,
              gd->m_RegularizationGlobalDataStruct,
              offset );

      // The actual update term for the normal component is
      // div(w^2nn^T grad(u_l^\perp))n_ln
      nln = normalVector[i] * normalVector;
      intermediateNormalRegularizationTerm
          = intermediateNormalRegularizationComponent * nln;
      normalRegularizationTerm
          = normalRegularizationTerm + intermediateNormalRegularizationTerm;
      }
    }

  // Will hold the final update term, including the regularization and intensity
  // distance terms
  PixelType updateTerm;
  updateTerm.Fill(0);
  updateTerm = intensityTermPlusTangentialRegularizationTerm
               + normalRegularizationTerm;
  return updateTerm;

}

} // end namespace itk

#endif
