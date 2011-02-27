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

  m_RegularizationFunction = RegularizationFunctionType::New();
  m_IntensityDistanceFunction = IntensityDistanceFunctionType::New();
  this->SetTimeStep( 1.0 );

  this->SetMovingImage(0);
  this->SetFixedImage(0);

  m_SumOfSquaredChange = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_RMSChange = NumericTraits< double >::max();
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

  os << indent << "Time step: " << m_TimeStep << std::endl;
  os << indent << "Compute regularization term: "
     << ( m_ComputeRegularizationTerm ? "on" : "off" ) << std::endl;
  os << indent << "Compute intensity distance term: "
     << ( m_ComputeIntensityDistanceTerm ? "on" : "off" ) << std::endl;
  if ( m_RegularizationFunction )
    {
    os << indent << "Regularization function: " << std::endl;
    m_RegularizationFunction->Print( os, indent );
    }
  if ( m_IntensityDistanceFunction )
    {
    os << indent << "Intensity distance function: " << std::endl;
    m_IntensityDistanceFunction->Print( os, indent );
    }
  os << indent << "Sum of squared change: " << m_SumOfSquaredChange
      << std::endl;
  os << indent << "Number of pixels processed: " << m_NumberOfPixelsProcessed
      << std::endl;
  os << indent << "RMS change: " << m_RMSChange << std::endl;
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

  ans->m_SumOfSquaredChange = 0;
  ans->m_NumberOfPixelsProcessed = 0L;

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

  // Update the variables used to calculate RMS change
  m_SumOfSquaredChange += gd->m_SumOfSquaredChange;
  m_NumberOfPixelsProcessed += gd->m_NumberOfPixelsProcessed;
  if( m_NumberOfPixelsProcessed )
    {
    m_RMSChange
        = vcl_sqrt( m_SumOfSquaredChange
                    / static_cast< double >( m_NumberOfPixelsProcessed ) );
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

  // initialize RMS calculation variables
  m_SumOfSquaredChange = 0.0;
  m_NumberOfPixelsProcessed = 0L;
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
                     << "be called.  Use the other ComputeUpdate() defined in"
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
    const DiffusionTensorNeighborhoodVectorType & tensorNeighborhoods,
    const ScalarDerivativeImageRegionArrayVectorType
        & deformationComponentFirstOrderDerivativeRegions,
    const TensorDerivativeImageRegionArrayVectorType
        & deformationComponentSecondOrderDerivativeRegions,
    const TensorDerivativeImageRegionVectorType & tensorDerivativeRegions,
    const DeformationVectorImageRegionArrayVectorType
        & multiplicationVectorRegionArrays,
    void * globalData,
    const FloatOffsetType & offset )
{
  // Get the global data structure
  GlobalDataStruct * gd = ( GlobalDataStruct * ) globalData;
  assert( gd );

  // Iterate over the deformation field components to compute the regularization
  // and intensity distance terms - note that PixelType corresponds to a
  // deformation vector

  // Compute the intensity distance update update term
  PixelType intensityDistanceTerm;
  intensityDistanceTerm.Fill(0);
  if ( this->GetComputeIntensityDistanceTerm() )
    {
    intensityDistanceTerm = m_IntensityDistanceFunction->ComputeUpdate(
        neighborhood, gd->m_IntensityDistanceGlobalDataStruct, offset );
    }

  // Compute the motion field regularization update term
  PixelType regularizationTerm;
  regularizationTerm.Fill(0);
  if ( this->GetComputeRegularizationTerm() )
    {
    int numTerms = tensorNeighborhoods.size();
    assert( (int) tensorDerivativeRegions.size() == numTerms );

    DeformationVectorComponentType intermediateComponent = 0;
    PixelType intermediateVector;
    intermediateVector.Fill(0);

    // Iterate over each div(T \grad(u))v term
    for ( int term = 0; term < numTerms; term++ )
      {
      assert( tensorNeighborhoods[term].GetImagePointer() );
      assert( tensorDerivativeRegions[term].GetImage() );
      // we don't necessarily have vectors to multiply, so no assert required

      // Iterate over each dimension
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        assert( deformationComponentFirstOrderDerivativeRegions[term][i].
                GetImage() );
        assert( deformationComponentSecondOrderDerivativeRegions[term][i].
                GetImage() );

        // Compute div(T \grad(u))
        intermediateComponent = m_RegularizationFunction->ComputeUpdate(
            tensorNeighborhoods[term],
            deformationComponentFirstOrderDerivativeRegions[term][i],
            deformationComponentSecondOrderDerivativeRegions[term][i],
            tensorDerivativeRegions[term],
            gd->m_RegularizationGlobalDataStruct,
            offset );

        // Multiply by the vector, if given
        intermediateVector.Fill(0);
        if( multiplicationVectorRegionArrays[term][i].GetImage() )
          {
          intermediateVector = intermediateComponent *
                               multiplicationVectorRegionArrays[term][i].Get();
          }
        else
          {
          intermediateVector[i] = intermediateComponent;
          }
        regularizationTerm += intermediateVector;
        }
      }
    }

  // Compute the final update term
  PixelType updateTerm;
  updateTerm.Fill(0);
  updateTerm = intensityDistanceTerm + regularizationTerm;

  // Update the variables used to calculate RMS change
  gd->m_NumberOfPixelsProcessed += 1;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    gd->m_SumOfSquaredChange += vnl_math_sqr( updateTerm[i] );
    }

  return updateTerm;
}

} // end namespace itk

#endif
