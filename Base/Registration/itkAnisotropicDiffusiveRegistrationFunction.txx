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

  m_IntensityDistanceWeighting = 1.0;
  m_RegularizationWeighting = 1.0;

  m_NumberOfPixelsProcessed = 0L;
  m_SumOfSquaredTotalChange = 0.0;
  m_SumOfSquaredIntensityDistanceChange = 0.0;
  m_SumOfSquaredRegularizationChange = 0.0;
  m_RMSTotalChange = 0.0;
  m_RMSIntensityDistanceChange = 0.0;
  m_RMSRegularizationChange = 0.0;
  m_RegularizationEnergy = 0.0;
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
  os << indent << "Intensity distance weighting: "
     << m_IntensityDistanceWeighting << std::endl;
  os << indent << "Regularization weighting: " << m_RegularizationWeighting
     << std::endl;
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
  os << indent << "Number of pixels processed: " << m_NumberOfPixelsProcessed
     << std::endl;
  os << indent << "Sum of squared total change: " << m_SumOfSquaredTotalChange
     << std::endl;
  os << indent << "RMS total change: " << m_RMSTotalChange << std::endl;
  os << indent << "Sum of squared intensity distance change: "
     << m_SumOfSquaredIntensityDistanceChange << std::endl;
  os << indent << "RMS intensity distance change: "
     << m_RMSIntensityDistanceChange << std::endl;
  os << indent << "Sum of squared regularization change: "
     << m_SumOfSquaredRegularizationChange << std::endl;
  os << indent << "RMS regularization change: " << m_RMSRegularizationChange
     << std::endl;
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

  ans->m_SumOfSquaredTotalChange = 0;
  ans->m_SumOfSquaredIntensityDistanceChange = 0;
  ans->m_SumOfSquaredRegularizationChange = 0;
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
  m_MetricCalculationLock.Lock();
  m_SumOfSquaredTotalChange += gd->m_SumOfSquaredTotalChange;
  m_SumOfSquaredIntensityDistanceChange
      += gd->m_SumOfSquaredIntensityDistanceChange;
  m_SumOfSquaredRegularizationChange += gd->m_SumOfSquaredRegularizationChange;
  m_NumberOfPixelsProcessed += gd->m_NumberOfPixelsProcessed;
  if( m_NumberOfPixelsProcessed )
    {
    m_RMSTotalChange
        = vcl_sqrt( m_SumOfSquaredTotalChange
                    / static_cast< double >( m_NumberOfPixelsProcessed ) );
    m_RMSIntensityDistanceChange
        = vcl_sqrt( m_SumOfSquaredIntensityDistanceChange
                    / static_cast< double >( m_NumberOfPixelsProcessed ) );
    m_RMSRegularizationChange
        = vcl_sqrt( m_SumOfSquaredRegularizationChange
                    / static_cast< double >( m_NumberOfPixelsProcessed ) );
    }
  m_MetricCalculationLock.Unlock();

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

  // Initialize RMS calculation variables
  m_SumOfSquaredTotalChange = 0.0;
  m_SumOfSquaredIntensityDistanceChange = 0.0;
  m_SumOfSquaredRegularizationChange = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_RegularizationEnergy = 0.0;
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
    bool includeInMetricComputations,
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
        neighborhood,
        gd->m_IntensityDistanceGlobalDataStruct,
        includeInMetricComputations,
        offset );
    }

  // Setup for computing the regularization energy at this voxel
  // Since we are iterating over terms before iterating over x,y,z
  // we need to store the sum for each dimension
  std::vector< DeformationVectorType > termRegularizationEnergies;
  if ( includeInMetricComputations )
    {
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      termRegularizationEnergies.push_back( DeformationVectorType(0.0) );
      }
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

        // Compute this part of the regularization energy
        if( includeInMetricComputations )
          {
          DiffusionTensorType diffusionTensor
              = tensorNeighborhoods[term].GetCenterPixel();
          ScalarDerivativeType deformationComponentFirstOrderDerivative
              = deformationComponentFirstOrderDerivativeRegions[term][i].Get();
          itk::Vector< double, ImageDimension > multVector(0.0);
          for( int row = 0; row < ImageDimension; row++ )
            {
            for( int col = 0; col < ImageDimension; col++ )
              {
              multVector[row]
                  += diffusionTensor(row,col) * deformationComponentFirstOrderDerivative[col];
              }
            }
          termRegularizationEnergies[i] += multVector;
          }

        }
      }
    }

  // Weight the intensity and regularization terms
  intensityDistanceTerm *= m_IntensityDistanceWeighting;
  regularizationTerm *= m_RegularizationWeighting;

  // Compute the final update term
  PixelType updateTerm;
  updateTerm.Fill(0);
  updateTerm = intensityDistanceTerm + regularizationTerm;

  // Finish computing the regularization energies
  if( includeInMetricComputations )
    {
    m_EnergyCalculationLock.Lock();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      m_RegularizationEnergy +=
          0.5 * ( termRegularizationEnergies[i] * termRegularizationEnergies[i] );
      }
    m_EnergyCalculationLock.Unlock();

    // Update the variables used to calculate RMS change
    gd->m_NumberOfPixelsProcessed += 1;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      gd->m_SumOfSquaredTotalChange += vnl_math_sqr( updateTerm[i] );
      gd->m_SumOfSquaredIntensityDistanceChange
          += vnl_math_sqr( intensityDistanceTerm[i] );
      gd->m_SumOfSquaredRegularizationChange
          += vnl_math_sqr( regularizationTerm[i] );
      }
    }

  return updateTerm;
}

} // end namespace itk

#endif
