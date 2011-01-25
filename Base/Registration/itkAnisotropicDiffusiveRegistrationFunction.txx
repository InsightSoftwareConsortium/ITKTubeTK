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
  typename Superclass::RadiusType r;
  r.Fill(1);
  this->SetRadius(r);

  m_ComputeRegularizationTerm = true;
  m_ComputeIntensityDistanceTerm = true;
  m_UseAnisotropicRegularization = true;

  m_RegularizationFunction = RegularizationFunctionType::New();
  m_IntensityDistanceFunction = IntensityDistanceFunctionType::New();
  this->SetTimeStep( 1.0 );

  this->SetMovingImage(0);
  this->SetFixedImage(0);
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

  os << indent << "TimeStep: " << m_TimeStep;
  os << indent << "ComputeRegularizationTerm: "
      << this->GetComputeRegularizationTerm() << std::endl;
  os << indent << "ComputeIntensityDistanceTerm: "
      << this->GetComputeIntensityDistanceTerm() << std::endl;
  os << indent << "UseAnisotropicRegularization: "
      << this->GetUseAnisotropicRegularization() << std::endl;
  if ( m_RegularizationFunction )
    {
    os << indent << "RegularizationFunction: " << std::endl;
    m_RegularizationFunction->Print( os, indent );
    }
  if ( m_IntensityDistanceFunction )
    {
    os << indent << "IntensityDistanceFunction: " << std::endl;
    m_IntensityDistanceFunction->Print( os, indent );
    }
}

/**
 * Creates a pointer to the data structure used to manage global values
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void *
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::GetGlobalDataPointer() const
{
  GlobalDataStruct * ans = new GlobalDataStruct();

  // Create the component global data pointers
  if( this->GetComputeRegularizationTerm() )
    {
    ans->m_RegularizationGlobalDataStruct
        = m_RegularizationFunction->GetGlobalDataPointer();
    }
  if( this->GetComputeIntensityDistanceTerm() )
    {
    ans->m_IntensityDistanceGlobalDataStruct
        = m_IntensityDistanceFunction->GetGlobalDataPointer();
    }

  return ans;
}

/**
 * Deletes the global data structure
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ReleaseGlobalDataPointer(void *GlobalData) const
{
  GlobalDataStruct * gd = ( GlobalDataStruct * ) GlobalData;

  // Release the component data structures
  if( this->GetComputeRegularizationTerm() )
    {
    m_RegularizationFunction->ReleaseGlobalDataPointer(
        gd->m_RegularizationGlobalDataStruct );
    }
  if( this->GetComputeIntensityDistanceTerm() )
    {
    m_IntensityDistanceFunction->ReleaseGlobalDataPointer(
        gd->m_IntensityDistanceGlobalDataStruct );
    }

  delete gd;
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
  if( !this->GetMovingImage() || !this->GetFixedImage()
    || !this->GetDeformationField() )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or deformation field "
                       << "not set" );
    }

  // Setup and initialize the component functions
  if( this->GetComputeIntensityDistanceTerm() )
    {
    m_IntensityDistanceFunction->SetMovingImage( this->GetMovingImage() );
    m_IntensityDistanceFunction->SetFixedImage( this->GetFixedImage() );
    m_IntensityDistanceFunction->SetDeformationField(
        this->GetDeformationField() );
    m_IntensityDistanceFunction->InitializeIteration();
    }
  if( this->GetComputeRegularizationTerm() )
    {
    m_RegularizationFunction->InitializeIteration();
    }
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
::ComputeUpdate(const NeighborhoodType &, void *, const FloatOffsetType & )
{
  // This function should never be called!
  itkExceptionMacro( << "ComputeUpdate(neighborhood, gd, offset) should never"
                     << "be called.  Use another ComputeUpdate() defined in"
                     << "itkAnisotropicDiffusiveRegistrationFunction instead" );
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
  // Get the global data structure
  GlobalDataStruct * gd = ( GlobalDataStruct * ) globalData;

  // Iterate over the deformation field components to compute the regularization
  // and intensity distance terms - note that PixelType corresponds to a
  // deformation vector
  PixelType                         tangentialRegularizationTerm;
  tangentialRegularizationTerm.Fill(0);
  PixelType                         normalRegularizationTerm;
  normalRegularizationTerm.Fill(0);
  PixelType                         intensityDistanceTerm;
  intensityDistanceTerm.Fill(0);
  PixelType                         updateTerm;
  updateTerm.Fill(0);

  // Compute the intensity distance update update term
  if (this->GetComputeIntensityDistanceTerm() )
    {
    intensityDistanceTerm = m_IntensityDistanceFunction->ComputeUpdate(
        neighborhood, gd->m_IntensityDistanceGlobalDataStruct, offset );
    }

  // Compute the motion field regularization update term
  if ( this->GetComputeRegularizationTerm() )
    {
    NormalVectorType                  normalVector;
    DeformationVectorComponentType    intermediateNormalRegularizationComponent;
    PixelType                         intermediateNormalRegularizationTerm;
    NormalVectorType                  nln; // n(l)n

    // Get the normal at this pixel once
    if( this->GetUseAnisotropicRegularization() )
      {
      const typename FixedImageType::IndexType index = neighborhood.GetIndex();
      normalVector
          = normalVectorNeighborhood.GetImagePointer()->GetPixel( index );
      }

    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      // Compute the regularization in the tangential plane (this will be in the
      // entire 3D space if we are using the Gaussian regularization)
      tangentialRegularizationTerm[i]
          = m_RegularizationFunction->ComputeUpdate(
              tangentialDeformationComponentNeighborhoods[i],
              tangentialTensorNeighborhood,
              tangentialTensorDerivativeRegion,
              spacing,
              gd->m_RegularizationGlobalDataStruct,
              offset );

      if( this->GetUseAnisotropicRegularization() )
        {
        // Compute the regularization in the normal direction
        intermediateNormalRegularizationComponent
            = m_RegularizationFunction->ComputeUpdate(
                normalDeformationComponentNeighborhoods[i],
                normalTensorNeighborhood,
                normalTensorDerivativeRegion,
                spacing,
                gd->m_RegularizationGlobalDataStruct,
                offset );

        nln = normalVector[i] * normalVector;

        intermediateNormalRegularizationTerm
            = intermediateNormalRegularizationComponent * nln;
        normalRegularizationTerm
            = normalRegularizationTerm + intermediateNormalRegularizationTerm;
        }
      }
    }

  updateTerm = intensityDistanceTerm
               + tangentialRegularizationTerm + normalRegularizationTerm;
  return updateTerm;
}

} // end namespace itk

#endif
