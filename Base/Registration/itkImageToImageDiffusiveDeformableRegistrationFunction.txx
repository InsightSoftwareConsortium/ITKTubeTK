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

  this->SetMovingImage(0);
  this->SetFixedImage(0);

  m_RegularizationFunction = RegularizationFunctionType::New();
  m_RegularizationFunction->SetTimeStep(m_TimeStep);

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
  ans->m_RegularizationGlobalDataStruct = ( RegularizationGlobalDataStruct *)
                              m_RegularizationFunction->GetGlobalDataPointer();

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
  m_RegularizationFunction->ReleaseGlobalDataPointer( gd->m_RegularizationGlobalDataStruct );

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
  DiffusionTensorNeighborhoodType diffusionTensor;
  DeformationFieldComponentNeighborhoodArrayType deformationFieldComponentArray;
  return this->ComputeUpdate( neighborhood,
                              diffusionTensor,
                              deformationFieldComponentArray,
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
                const DiffusionTensorNeighborhoodType &neighborhoodTensor,
                const DeformationFieldComponentNeighborhoodArrayType
                                        &neighborhoodDeformationFieldComponents,
                void *gd,
                const FloatOffsetType& offset)
{  

  // Assertion to make sure the deformation field components are setup properly
  // TODO don't do this processing if not in debug mode
  DeformationFieldConstPointer deformationField
                                              = neighborhood.GetImagePointer();
  itk::FixedArray< DeformationFieldComponentImageConstPointer, ImageDimension >
                               deformationFieldComponentImageArray;
  typedef Index< ImageDimension > IndexType;
  const IndexType index = neighborhood.GetIndex();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    deformationFieldComponentImageArray[i]
                  = neighborhoodDeformationFieldComponents[i].GetImagePointer();
    assert( deformationField->GetPixel( index )[i]
                == deformationFieldComponentImageArray[i]->GetPixel( index ) );
    }

  // Get the global data structure
  GlobalDataStruct * globalData = ( GlobalDataStruct * ) gd;

  // Get the global data structure for the regularization
  RegularizationGlobalDataStruct * regularizationGlobalData =
                                  globalData->m_RegularizationGlobalDataStruct;

  // Iterate over the deformation field components to compute the regularization
  // term
  PixelType regularizationTerm;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    regularizationTerm[i] = m_RegularizationFunction->ComputeUpdate(
                                    neighborhoodDeformationFieldComponents[i],
                                    neighborhoodTensor,
                                    regularizationGlobalData,
                                    offset );
    }

  // TODO in future, add result of intensity difference
  return regularizationTerm;

}

} // end namespace itk

#endif
