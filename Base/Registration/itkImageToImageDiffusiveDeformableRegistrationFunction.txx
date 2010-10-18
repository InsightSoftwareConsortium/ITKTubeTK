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

  m_ComputeRegularizationTerm = true;
  m_ComputeIntensityDistanceTerm = true;

  m_RegularizationFunction = RegularizationFunctionType::New();
  m_RegularizationFunction->SetTimeStep( m_TimeStep );

  m_IntensityDistanceFunction = IntensityDistanceFunctionType::New();
  // TODO more parameters for the intensity distance function?
  // TODO make these options for this function
  // TODO NOW
  // ... inherited for finite
  //  m_IntensityDistanceFunction->SetScaleCoefficients( vals );
  // ... inherited from PDE function
  //  m_IntensityDistanceFunction->SetGradientStep( 0.0 );
  //  m_IntensityDistanceFunction->SetNormalizeGradient( false );




  // TODO default for now is LinearInterpolateImageFilter - can change later
  //  m_IntensityDistanceFunction->SetMovingImageInterpolator( interp );


  // TODO computes global timestep - will need to change once time step changes
  //   MeanSquareRegistration uses constant timestep of 1
  // but maybe ok because we're just calling compute update, which doesn't
  // use m_TimeStep
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
  os << indent << "ComputeRegularizationTerm: "
                                << m_ComputeRegularizationTerm << std::endl;
  os << indent << "ComputeIntensityDistanceTerm: "
                                << m_ComputeIntensityDistanceTerm << std::endl;
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
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::SetComputeRegularizationTerm( bool compute )
{
  m_ComputeRegularizationTerm = compute;
}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
bool
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::GetComputeRegularizationTerm() const
{
  return m_ComputeRegularizationTerm;
}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::SetComputeIntensityDistanceTerm( bool compute )
{
  m_ComputeIntensityDistanceTerm = compute;
}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
bool
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::GetComputeIntensityDistanceTerm() const
{
  return m_ComputeIntensityDistanceTerm;
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
  std::cout << "\tInitializeIteration for FUNCTION" << std::endl;

  if( !this->GetMovingImage()
      || !this->GetFixedImage()
      /*|| !this->GetDeformationField()*/ ) // TODO put back
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or deformation field not set" );
    }

  // Setup the component functions
  // TODO moving image and fixed image can be set only once
  // TODO more set methods to be applied to the intensity distance function?
  // TODO NOW
  m_IntensityDistanceFunction->SetMovingImage( this->GetMovingImage() ) ;
  m_IntensityDistanceFunction->SetFixedImage( this->GetFixedImage() );
  m_IntensityDistanceFunction->SetDeformationField( this->GetDeformationField() );

  // Initialize the component functions
  m_RegularizationFunction->InitializeIteration();
  m_IntensityDistanceFunction->InitializeIteration();
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
  // TODO this likely won't work, but this function shouldn't be called anyways
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
  NormalVectorType normalVector
          = normalVectorImageNeighborhood.GetImagePointer()->GetPixel( index );

  // Iterate over the deformation field components to compute the regularization
  // and intensity distance terms
  PixelType                         tangentialRegularizationTerm;
  PixelType                         normalRegularizationTerm;
  DeformationFieldScalarType        intermediateNormalRegularizationComponent;
  PixelType                         intermediateNormalRegularizationTerm;
  PixelType                         intensityDistanceTerm;
  PixelType                         updateTerm;
  NormalVectorType                  nln; // n(l)n

  intensityDistanceTerm.Fill(0);
  tangentialRegularizationTerm.Fill(0);
  normalRegularizationTerm.Fill(0); // essential because incremented in loop

  // Compute the intensity distance update
  // TODO make sure that there is no normalization in the intensity distance
  // function - check that smooth gradient is off, and that the normalize
  // metric doesn't do anything we don't want
  // TODO weighting between the intensity distance and regularization terms?

  if ( m_ComputeIntensityDistanceTerm )
    {
    intensityDistanceTerm = m_IntensityDistanceFunction->ComputeUpdate(
                              neighborhood,
                              globalData->m_IntensityDistanceGlobalDataStruct,
                              offset );
    }

  // Compute the motion field regularization
  if (m_ComputeRegularizationTerm )
    {
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      // Compute the regularization in the tangential plane
      tangentialRegularizationTerm[i]
                              = m_RegularizationFunction->ComputeUpdate(
                              tangentialNeighborhoodDeformationFieldComponents[i],
                              tangentialNeighborhoodTensor,
                              globalData->m_RegularizationGlobalDataStruct,
                              offset );

      // Compute the regularization in the normal direction
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

    }



  // TODO weighting here?  Don't worry about weighting if one term is not
  // computed because of boolean settings
  updateTerm = intensityDistanceTerm
                      + tangentialRegularizationTerm + normalRegularizationTerm;

  // TODO take out
//  if ( updateTerm[0] != 0 || updateTerm[1] != 0 || updateTerm[2] != 0 )
//    {
//    std::cout << "update term " << updateTerm[0] << " "
//                                << updateTerm[1] << " "
//                                << updateTerm[2] << std::endl;
//    }

  return updateTerm;

}

} // end namespace itk

#endif
