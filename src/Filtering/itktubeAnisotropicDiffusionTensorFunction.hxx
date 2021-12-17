/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#ifndef __itktubeAnisotropicDiffusionTensorFunction_hxx
#define __itktubeAnisotropicDiffusionTensorFunction_hxx


#include <itkTextOutput.h>

namespace itk
{

namespace tube
{

template< class TImageType >
AnisotropicDiffusionTensorFunction< TImageType>
::AnisotropicDiffusionTensorFunction( void )
{
  typename Superclass::RadiusType r;
  r.Fill( 1 );
  this->SetRadius( r );

  itk::OutputWindow::SetInstance( itk::TextOutput::New() );

  // Dummy neighborhood.
  NeighborhoodType it;
  it.SetRadius( r );

  // Find the center index of the neighborhood.
  m_Center = static_cast< unsigned int >( it.Size() / 2 );

  // Get the stride length for each axis.
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_xStride[i] = static_cast< unsigned int >( it.GetStride( i ) );
    }

  // Calculate the required indexes surrounding the center position once.
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_positionA[i] = m_Center + m_xStride[i];
    m_positionB[i] = m_Center - m_xStride[i];
    for( unsigned int j = i+1; j < ImageDimension; j++ )
      {
      m_positionAa[i][j] = m_Center - m_xStride[i] - m_xStride[j];
      m_positionBa[i][j] = m_Center - m_xStride[i] + m_xStride[j];
      m_positionCa[i][j] = m_Center + m_xStride[i] - m_xStride[j];
      m_positionDa[i][j] = m_Center + m_xStride[i] + m_xStride[j];
      }
    }

  // Whether or not to use the image spacing when computing the derivatives
  this->m_UseImageSpacing = true;
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction<TImageType>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate( const NeighborhoodType &neighborhood,
                void *globalData,
                const FloatOffsetType& offset )
{
  DiffusionTensorNeighborhoodType tensorNeighborhood;
  SpacingType                     spacing;
  return this->ComputeUpdate( neighborhood,
                              tensorNeighborhood,
                              spacing,
                              globalData,
                              offset );
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate( const NeighborhoodType &neighborhood,
                const DiffusionTensorNeighborhoodType &tensorNeighborhood,
                const SpacingType &spacing,
                void *globalData,
                const FloatOffsetType& itkNotUsed( offset ) )
{
  // Global data structure
  GlobalDataStruct *gd = static_cast<GlobalDataStruct *>( globalData );
  assert( gd );
  gd->m_GradMagSqr = 1.0e-6;

  // Compute intensity first and second order partial derivatives
  this->ComputeIntensityFirstAndSecondOrderPartialDerivatives(
      neighborhood, gd->m_dx, gd->m_dxy, spacing );

  // Compute diffusion tensor matrix first order partial derivatives
  this->ComputeDiffusionTensorFirstOrderPartialDerivatives(
      tensorNeighborhood, gd->m_DT_dxy, spacing );

  // Compute the update term
  return this->ComputeFinalUpdateTerm( tensorNeighborhood, gd );
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate(
    const DiffusionTensorNeighborhoodType &tensorNeighborhood,
    const ScalarDerivativeImageRegionType &intensityFirstDerivatives,
    const TensorDerivativeImageRegionType &intensitySecondDerivatives,
    const TensorDerivativeImageRegionType &tensorFirstDerivatives,
    void *globalData,
    const FloatOffsetType& itkNotUsed( offset ) )
{
  // Global data structure
  GlobalDataStruct *gd = static_cast<GlobalDataStruct *>( globalData );
  assert( gd );
  gd->m_GradMagSqr = 1.0e-6;

  // Copy the intensity first and second order partial derivatives into the
  // global data struct
  gd->m_dx = intensityFirstDerivatives.Get();
  gd->m_dxy = intensitySecondDerivatives.Get();

  // Copy the diffusion tensor matrix first order partial derivatives into the
  // global data struct
  gd->m_DT_dxy = tensorFirstDerivatives.Get();

  // Compute the update term
  return this->ComputeFinalUpdateTerm( tensorNeighborhood, gd );
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeDiffusionTensorFirstOrderPartialDerivatives(
    const DiffusionTensorNeighborhoodType &tensorNeighborhood,
    TensorDerivativeType &firstOrderResult,
    const SpacingType &spacing ) const
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    DiffusionTensorType positionA_Tensor_value
        = tensorNeighborhood.GetPixel( m_positionA[i] );
    DiffusionTensorType positionB_Tensor_value
        = tensorNeighborhood.GetPixel( m_positionB[i] );

    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      firstOrderResult( i, j ) = 0.5 * ( positionA_Tensor_value( i, j )
        - positionB_Tensor_value( i, j ) );

      // Handle image spacing
      if( m_UseImageSpacing )
        {
        firstOrderResult( i, j ) /= spacing[i];
        }
      }
    }
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeDiffusionTensorFirstOrderPartialDerivatives(
    const DiffusionTensorNeighborhoodType &tensorNeighborhood,
    TensorDerivativeImageRegionType &firstOrderResult,
    const SpacingType &spacing ) const
{
  this->ComputeDiffusionTensorFirstOrderPartialDerivatives(
      tensorNeighborhood, firstOrderResult.Value(), spacing );
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeIntensityFirstAndSecondOrderPartialDerivatives(
    const NeighborhoodType &neighborhood,
    ScalarDerivativeType &firstOrderResult,
    TensorDerivativeType &secondOrderResult,
    const SpacingType &spacing ) const
{
  const ScalarValueType center_value = neighborhood.GetCenterPixel();

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    const ScalarValueType it_positionA
        = neighborhood.GetPixel( m_positionA[i] );
    const ScalarValueType it_positionB
        = neighborhood.GetPixel( m_positionB[i] );

    // First order partial derivatives
    firstOrderResult[i] = 0.5 * ( it_positionA - it_positionB );

    // Second order partial derivatives
    secondOrderResult[i][i] = ( it_positionA + it_positionB -
      2.0 * center_value );

    for( unsigned int j = i+1; j < ImageDimension; j++ )
      {
      secondOrderResult[i][j]
          = secondOrderResult[j][i] // Guaranteed symmetric
            = 0.25 * ( neighborhood.GetPixel( m_positionAa[i][j] )
                - neighborhood.GetPixel( m_positionBa[i][j] )
                - neighborhood.GetPixel( m_positionCa[i][j] )
                + neighborhood.GetPixel( m_positionDa[i][j] ) );
      }
    }

  // Handle image spacing
  if( m_UseImageSpacing )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      firstOrderResult[i] /= spacing[i];
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        secondOrderResult[i][j] /= ( spacing[i] * spacing[j] );
        }
      }
    }
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeIntensityFirstAndSecondOrderPartialDerivatives(
    const NeighborhoodType &neighborhood,
    ScalarDerivativeImageRegionType &firstOrderResult,
    TensorDerivativeImageRegionType &secondOrderResult,
    const SpacingType &spacing ) const
{
 this->ComputeIntensityFirstAndSecondOrderPartialDerivatives(
     neighborhood,
     firstOrderResult.Value(),
     secondOrderResult.Value(),
     spacing );
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeFinalUpdateTerm(
    const DiffusionTensorNeighborhoodType &tensorNeighborhood,
    const GlobalDataStruct *gd ) const
{
  DiffusionTensorType center_Tensor_value = tensorNeighborhood.GetCenterPixel();

  ScalarValueType pdWrtDiffusion1 = gd->m_DT_dxy[0][0] * gd->m_dx[0]
                                    + gd->m_DT_dxy[0][1] * gd->m_dx[1]
                                    + gd->m_DT_dxy[0][2] * gd->m_dx[2];

  ScalarValueType pdWrtDiffusion2 = gd->m_DT_dxy[1][0] * gd->m_dx[0]
                                    + gd->m_DT_dxy[1][1] * gd->m_dx[1]
                                    + gd->m_DT_dxy[1][2] * gd->m_dx[2];

  ScalarValueType pdWrtDiffusion3 = gd->m_DT_dxy[2][0] * gd->m_dx[0]
                                    + gd->m_DT_dxy[2][1] * gd->m_dx[1]
                                    + gd->m_DT_dxy[2][2] * gd->m_dx[2];

  ScalarValueType pdWrtImageIntensity1
      = center_Tensor_value( 0, 0 ) * gd->m_dxy[0][0]
        + center_Tensor_value( 0, 1 ) * gd->m_dxy[0][1]
        + center_Tensor_value( 0, 2 ) * gd->m_dxy[0][2];

  ScalarValueType pdWrtImageIntensity2
      = center_Tensor_value( 1, 0 ) * gd->m_dxy[1][0]
        + center_Tensor_value( 1, 1 ) * gd->m_dxy[1][1]
        + center_Tensor_value( 1, 2 ) * gd->m_dxy[1][2];

  ScalarValueType pdWrtImageIntensity3
      = center_Tensor_value( 2, 0 ) * gd->m_dxy[2][0]
        + center_Tensor_value( 2, 1 ) * gd->m_dxy[2][1]
        + center_Tensor_value( 2, 2 ) * gd->m_dxy[2][2];

  ScalarValueType total
  = pdWrtDiffusion1 + pdWrtDiffusion2 + pdWrtDiffusion3
        + pdWrtImageIntensity1 + pdWrtImageIntensity2 + pdWrtImageIntensity3;

  return ( PixelType ) ( total );
}

template< class TImageType >
template< class TPixel, unsigned int VImageDimension >
void
AnisotropicDiffusionTensorFunction<TImageType>
::CheckTimeStepStability(
    const itk::Image< TPixel, VImageDimension > * input,
    bool useImageSpacing )
{
  double minSpacing;
  if( useImageSpacing )
    {
    minSpacing = input->GetSpacing()[0];
    for( unsigned int i = 1; i < ImageDimension; i++ )
      {
      if( input->GetSpacing()[i] < minSpacing )
        {
        minSpacing = input->GetSpacing()[i];
        }
      }
    }
  else
    {
    minSpacing = 1.0;
    }

  // TODO plus 1?
  double ratio =
    minSpacing / std::pow( 2.0, static_cast<double>( ImageDimension ) + 1 );

  if( m_TimeStep > ratio )
    {
    itkWarningMacro( << std::endl
      << "Anisotropic diffusion unstable time step:"
      << m_TimeStep << std::endl << "Minimum stable time step"
      << "for this image is " << ratio );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeAnisotropicDiffusionTensorFunction_hxx )
