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
#ifndef _itkImageToImageDiffusiveDeformableRegistrationFunction_txx
#define _itkImageToImageDiffusiveDeformableRegistrationFunction_txx

#include "itkImageToImageDiffusiveDeformableRegistrationFunction.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ImageToImageDiffusiveDeformableRegistrationFunction()
{
  // TODO will be computed
  m_TimeStep = 0.05;

  RadiusType r;
  r.Fill(1);
  this->SetRadius(r);

  m_RegularizationFunction = RegularizationFunctionType::New();
  m_RegularizationFunction->SetTimeStep( m_TimeStep );

  m_IntensityDistanceFunction = IntensityDistanceFunctionType::New();
  // TODO computes global timestep - is this a problem?
  //m_IntensityDistanceFunction->SetTimeStep( m_TimeStep );

  this->SetMovingImage(0);
  this->SetFixedImage(0);

}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "TimeStep: " << m_TimeStep;
  std::cout << std::endl;
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
 * Sets the fixed image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::SetFixedImage( const FixedImageType * ptr )
{
  Superclass::SetFixedImage( ptr );

  if ( m_IntensityDistanceFunction )
    {
    m_IntensityDistanceFunction->SetFixedImage( ptr );
    }
}

/**
 * Sets the moving image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::SetMovingImage( const MovingImageType * ptr )
{
  Superclass::SetMovingImage( ptr );

  if ( m_IntensityDistanceFunction )
    {
    m_IntensityDistanceFunction->SetMovingImage( ptr );
    }
}

/**
 * Sets the deformation field
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::SetDeformationField( const DeformationFieldTypePointer ptr )
{
  Superclass::SetDeformationField( ptr );

  if ( m_IntensityDistanceFunction )
    {
    m_IntensityDistanceFunction->SetDeformationField( ptr );
    }

}

/**
 * Creates a pointer to the data structure used to manage global values
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void *
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::GetGlobalDataPointer() const
{
  GlobalDataStruct * ans = new GlobalDataStruct();

  // Create the component global data pointers
  ans->m_RegularizationGlobalDataStruct =
                            m_RegularizationFunction->GetGlobalDataPointer();
  ans->m_IntensityDistanceGlobalDataStruct =
                            m_IntensityDistanceFunction->GetGlobalDataPointer();

  return ans;
}

/**
  * Deletes the global data structure
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ReleaseGlobalDataPointer(void *GlobalData) const
{
  GlobalDataStruct * gd = ( GlobalDataStruct * ) GlobalData;

  // Release the component data structures
  m_RegularizationFunction->ReleaseGlobalDataPointer(
                                    gd->m_RegularizationGlobalDataStruct );
  m_IntensityDistanceFunction->ReleaseGlobalDataPointer(
                                    gd->m_IntensityDistanceGlobalDataStruct );

  delete gd;
}


/**
  * Called at the beginning of each iteration
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::InitializeIteration()
{
  std::cout << "InitializeIteration for FUNCTION" << std::endl;

  if( !this->GetMovingImage()
      || !this->GetFixedImage()
      /*|| !m_MovingImageInterpolator*/ ) // m_TangentalDiffusionTensorImage interpolator
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
    }

  // Initialize the component functions
  m_RegularizationFunction->InitializeIteration();
  // TODO put back
  //m_IntensityDistanceFunction->InitializeIteration();
}

/**
  * Computes the update term
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename ImageToImageDiffusiveDeformableRegistrationFunction
                                < TFixedImage, TMovingImage, TDeformationField >
::PixelType
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ComputeUpdate(const NeighborhoodType &neighborhood,
                void *gd,
                const FloatOffsetType& offset)
{
  // TODO this likely won't work, but shouldn't be called anyways
  NormalVectorImageNeighborhoodType   normalVectorImageNeighborhood;
  DiffusionTensorNeighborhoodType     tangentialDiffusionTensor;
  DeformationFieldComponentNeighborhoodArrayType
                                      tangentialDeformationFieldComponentArray;
  DiffusionTensorNeighborhoodType     normalDiffusionTensor;
  DeformationFieldComponentNeighborhoodArrayType
                                      normalDeformationFieldComponentArray;
  return this->ComputeUpdate( neighborhood,
                              normalVectorImageNeighborhood,
                              tangentialDiffusionTensor,
                              tangentialDeformationFieldComponentArray,
                              normalDiffusionTensor,
                              normalDeformationFieldComponentArray,
                              gd,
                              offset);
}

/**
  * Computes the update term
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename ImageToImageDiffusiveDeformableRegistrationFunction
                                < TFixedImage, TMovingImage, TDeformationField >
::PixelType
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ComputeUpdate(const NeighborhoodType &neighborhood,
                const NormalVectorImageNeighborhoodType
                              &normalVectorImageNeighborhood,
                const DiffusionTensorNeighborhoodType
                              &tangentialNeighborhoodTensor,
                const DeformationFieldComponentNeighborhoodArrayType
                              &tangentialNeighborhoodDeformationFieldComponents,
                const DiffusionTensorNeighborhoodType
                              &normalNeighborhoodTensor,
                const DeformationFieldComponentNeighborhoodArrayType
                              &normalNeighborhoodDeformationFieldComponents,
                void *gd,
                const FloatOffsetType& offset)
{  
  // Get the global data structure
  GlobalDataStruct * globalData = ( GlobalDataStruct * ) gd;

  // Get the normal at this pixel
  const IndexType index = normalVectorImageNeighborhood.GetIndex();
  NormalVectorType normalVector = normalVectorImageNeighborhood.GetImagePointer()->GetPixel( index );

  // Iterate over the deformation field components to compute the regularization
  // and intensity distance terms
  PixelType                         tangentialRegularizationTerm;
  PixelType                         normalRegularizationTerm;
  DeformationFieldScalarType        intermediateNormalRegularizationComponent;
  PixelType                         intermediateNormalRegularizationTerm;
  PixelType                         intensityDistanceTerm;
  PixelType                         updateTerm;
  NormalVectorType                  nln; // n(l)n
  normalRegularizationTerm.Fill(0);
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    // TODO compute the intensity distance function
    intensityDistanceTerm[i] = 0;

    // compute the regularization in the tangential plane
    tangentialRegularizationTerm[i]
                            = m_RegularizationFunction->ComputeUpdate(
                            tangentialNeighborhoodDeformationFieldComponents[i],
                            tangentialNeighborhoodTensor,
                            globalData->m_RegularizationGlobalDataStruct,
                            offset );

    // compute the regularization in the normal direction
    intermediateNormalRegularizationComponent
                            = m_RegularizationFunction->ComputeUpdate(
                            normalNeighborhoodDeformationFieldComponents[i],
                            normalNeighborhoodTensor,
                            globalData->m_RegularizationGlobalDataStruct,
                            offset );

    nln = normalVector[i] * normalVector;
    intermediateNormalRegularizationTerm
              = intermediateNormalRegularizationComponent * nln;
    normalRegularizationTerm
              = normalRegularizationTerm + intermediateNormalRegularizationTerm;
    }

  // TODO in future, add result of intensity difference
  updateTerm = intensityDistanceTerm
                      + tangentialRegularizationTerm + normalRegularizationTerm;
  return updateTerm;

}

} // end namespace itk

#endif
