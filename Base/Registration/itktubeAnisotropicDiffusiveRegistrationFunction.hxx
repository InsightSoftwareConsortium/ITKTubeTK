/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeAnisotropicDiffusiveRegistrationFunction_hxx
#define __itktubeAnisotropicDiffusiveRegistrationFunction_hxx

#include "itktubeAnisotropicDiffusiveRegistrationFunction.h"

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
AnisotropicDiffusiveRegistrationFunction
 < TFixedImage, TMovingImage, TDeformationField >
::AnisotropicDiffusiveRegistrationFunction( void )
{
  typename Superclass::RadiusType r;
  r.Fill( 1 );
  this->SetRadius( r );

  m_ComputeRegularizationTerm = true;
  m_ComputeIntensityDistanceTerm = true;

  m_RegularizationFunction = RegularizationFunctionType::New();

  m_IntensityDistanceFunction = IntensityDistanceFunctionType::New();
  m_IntensityDistanceFunction->SetBackgroundIntensity( 0.0 );
  m_IntensityDistanceFunction->SetNormalizeGradient( false );

  this->SetTimeStep( 1.0 );
  this->SetMovingImage( 0 );
  this->SetFixedImage( 0 );

  m_RegularizationWeighting = 1.0;
  m_RegularizationEnergy = 0.0;
}

/**
 * PrintSelf
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Time step: " << m_TimeStep << std::endl;
  os << indent << "Compute regularization term: "
     << ( m_ComputeRegularizationTerm ? "on" : "off" ) << std::endl;
  os << indent << "Compute intensity distance term: "
     << ( m_ComputeIntensityDistanceTerm ? "on" : "off" ) << std::endl;
  os << indent << "Regularization weighting: " << m_RegularizationWeighting
     << std::endl;
  if( m_RegularizationFunction )
    {
    os << indent << "Regularization function: " << std::endl;
    m_RegularizationFunction->Print( os, indent );
    }
  if( m_IntensityDistanceFunction )
    {
    os << indent << "Intensity distance function: " << std::endl;
    m_IntensityDistanceFunction->Print( os, indent );
    }
}

/**
 * Creates a pointer to the data structure used to manage global values
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void *
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::GetGlobalDataPointer( void ) const
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
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ReleaseGlobalDataPointer( void *GlobalData ) const
{
  GlobalDataStruct * gd = static_cast<GlobalDataStruct *>( GlobalData );

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
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeIteration( void )
{
  if( !this->GetMovingImage() || !this->GetFixedImage()
    || !this->GetDisplacementField() )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or deformation field "
                       << "not set" );
    }

  // Setup and initialize the component functions
  if( this->GetComputeIntensityDistanceTerm() )
    {
    m_IntensityDistanceFunction->SetMovingImage( this->GetMovingImage() );
    m_IntensityDistanceFunction->SetFixedImage( this->GetFixedImage() );
    m_IntensityDistanceFunction->SetDisplacementField(
        this->GetDisplacementField() );
    m_IntensityDistanceFunction->InitializeIteration();
    }
  if( this->GetComputeRegularizationTerm() )
    {
    m_RegularizationFunction->InitializeIteration();
    }
  m_RegularizationEnergy = 0.0;
}

/**
 * Computes the update term
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::PixelType
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeUpdate( const NeighborhoodType &, void *, const FloatOffsetType & )
{
  // This function should never be called!
  itkExceptionMacro(
    << "ComputeUpdate( neighborhood, gd, offset ) should never"
    << "be called.  Use the other ComputeUpdate() defined in"
    << "itktubeAnisotropicDiffusiveRegistrationFunction instead" );
}

/**
  * Computes the update term
  */
template< class TFixedImage, class TMovingImage, class TDeformationField >
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
    PixelType & intensityDistanceTerm,
    PixelType & regularizationTerm,
    const FloatOffsetType & offset )
{
  // Get the global data structure
  GlobalDataStruct * gd = static_cast<GlobalDataStruct *>( globalData );
  assert( gd );

  // Iterate over the deformation field components to compute the regularization
  // and intensity distance terms - note that PixelType corresponds to a
  // deformation vector

  // Compute the intensity distance update update term
  intensityDistanceTerm.Fill( 0 );
  if( this->GetComputeIntensityDistanceTerm() )
    {
    intensityDistanceTerm = m_IntensityDistanceFunction->ComputeUpdate(
        neighborhood,
        gd->m_IntensityDistanceGlobalDataStruct,
        offset );
    }

  // Compute the ( weighted ) motion field regularization update term
  regularizationTerm.Fill( 0 );
  if( this->GetComputeRegularizationTerm() )
    {
    regularizationTerm = this->ComputeRegularizationUpdate(
          tensorNeighborhoods,
          deformationComponentFirstOrderDerivativeRegions,
          deformationComponentSecondOrderDerivativeRegions,
          tensorDerivativeRegions,
          multiplicationVectorRegionArrays,
          gd->m_RegularizationGlobalDataStruct,
          offset );
    }

  // Compute the final update term - incorporates weighting between intensity
  // distance term and regularization term, but not global timestep scaling
  PixelType updateTerm;
  updateTerm.Fill( 0 );
  updateTerm = intensityDistanceTerm + regularizationTerm;

  return updateTerm;
}

/**
  * Computes the update term for the regularization
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::PixelType
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeRegularizationUpdate(
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
  PixelType regularizationTerm;
  regularizationTerm.Fill( 0 );

  DeformationVectorComponentType intermediateComponent = 0;
  PixelType intermediateVector;
  intermediateVector.Fill( 0 );

  int numTerms = tensorNeighborhoods.size();
  assert( ( int ) tensorDerivativeRegions.size() == numTerms );

  // Iterate over each div( T \grad( u ) )v term
  for( int term = 0; term < numTerms; term++ )
    {
    assert( tensorNeighborhoods[term].GetImagePointer() );
    assert( tensorDerivativeRegions[term].GetImage() );
    // we don't necessarily have vectors to multiply, so no assert required

    // Iterate over each dimension
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      assert( deformationComponentFirstOrderDerivativeRegions[term][i].
              GetImage() );
      assert( deformationComponentSecondOrderDerivativeRegions[term][i].
              GetImage() );

      // Compute div( T \grad( u ) )
      intermediateComponent = m_RegularizationFunction->ComputeUpdate(
          tensorNeighborhoods[term],
          deformationComponentFirstOrderDerivativeRegions[term][i],
          deformationComponentSecondOrderDerivativeRegions[term][i],
          tensorDerivativeRegions[term],
          globalData,
          offset );

      // Multiply by the vector, if given
      intermediateVector.Fill( 0 );
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

  // Weight the regularization term, but don't multiply by the timestep
  // because the PDE framework will do that for you
  regularizationTerm *= m_RegularizationWeighting;

  return regularizationTerm;
}

/**
  * Computes the intensity distance energy
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
double
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeIntensityDistanceEnergy(
  const typename NeighborhoodType::IndexType index,
  const DeformationVectorType & update )
{
  return vnl_math_sqr( m_IntensityDistanceFunction->ComputeIntensityDifference(
                         index, update ) );
}

/**
  * Computes the regularization energy
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
double
AnisotropicDiffusiveRegistrationFunction
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeRegularizationEnergy(
    const DiffusionTensorNeighborhoodVectorType & tensorNeighborhoods,
    const ScalarDerivativeImageRegionArrayVectorType
        & deformationComponentFirstOrderDerivativeRegions )
{
  // Since we are iterating over terms before iterating over x,y,z
  // we need to store the sum for each dimension
  std::vector< DeformationVectorType > termRegularizationEnergies;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    termRegularizationEnergies.push_back( DeformationVectorType( 0.0 ) );
    }

  int numTerms = tensorNeighborhoods.size();
  for( int term = 0; term < numTerms; term++ )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      DiffusionTensorType diffusionTensor
          = tensorNeighborhoods[term].GetCenterPixel();
      ScalarDerivativeType deformationComponentFirstOrderDerivative
          = deformationComponentFirstOrderDerivativeRegions[term][i].Get();
      itk::Vector< double, ImageDimension > multVector( 0.0 );
      for( unsigned int row = 0; row < ImageDimension; row++ )
        {
        for( unsigned int col = 0; col < ImageDimension; col++ )
          {
          multVector[row] += diffusionTensor( row, col )
            * deformationComponentFirstOrderDerivative[col];
          }
        }
      termRegularizationEnergies[i] += multVector;
      }
    }

  // Calculate and weight the regularization energy
  double regularizationEnergy = 0.0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    regularizationEnergy += 0.5 *
      ( termRegularizationEnergies[i] * termRegularizationEnergies[i] );
    }
  regularizationEnergy *= m_RegularizationWeighting;

  return regularizationEnergy;
}

} // End namespace tube

} // End namespace itk

// End !defined( __itktubeAnisotropicDiffusiveRegistrationFunction_hxx )
#endif
