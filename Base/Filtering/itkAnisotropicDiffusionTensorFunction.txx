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
#ifndef __itkAnisotropicDiffusionTensorFunction_txx
#define __itkAnisotropicDiffusionTensorFunction_txx

#include "itkAnisotropicDiffusionTensorFunction.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

namespace itk {

template< class TImageType >
AnisotropicDiffusionTensorFunction< TImageType>
::AnisotropicDiffusionTensorFunction()
{
  typename Superclass::RadiusType r;

  for( unsigned int i=0; i < ImageDimension; i++ )
    {
    r[i] = 1;
    }

  this->SetRadius(r);

  // Dummy neighborhood.
  NeighborhoodType it;
  it.SetRadius( r );

  // Find the center index of the neighborhood.
  m_Center = static_cast<unsigned int>( it.Size() / 2 );

  // Get the stride length for each axis.
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {  m_xStride[i] = static_cast<unsigned int>( it.GetStride(i) ); }

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
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >
::TimeStepType
AnisotropicDiffusionTensorFunction<TImageType>
::ComputeGlobalTimeStep(void *) const
{
  /* returns the time step supplied by the user. We don't need
     to use the global data supplied since we are returning a fixed value
  */

  return this->GetTimeStep();
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate(const NeighborhoodType &it, void *globalData,
                const FloatOffsetType& offset)
{
  DiffusionTensorNeighborhoodType diffusionTensor;
  return this->ComputeUpdate( it, diffusionTensor, globalData, offset );
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate(const NeighborhoodType &it,
                const DiffusionTensorNeighborhoodType &gt,
                void *globalData,
                const FloatOffsetType& itkNotUsed(offset) )
{
  // Global data structure
  // m_dx -> Intensity first derivative
  // m_dxy -> Intensity second derivative
  // m_DT_dxy -> Diffusion tensor first derivative
  GlobalDataStruct *gd = (GlobalDataStruct *)globalData;

  // Compute the first and 2nd derivative for the intensity images
  this->ComputeIntensityFirstAndSecondDerivatives( it, gd );

  // Compute the diffusion tensor matrix first derivatives if not provided
  this->ComputeDiffusionFirstDerivative( gt, gd );

  // Compute the update term
  return this->ComputeFinalUpdateTerm( gt, gd );
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate(const NeighborhoodType &it,
                const DiffusionTensorNeighborhoodType &gt,
                const DerivativeMatrixImageRegionType &tensorDerivative,
                void *globalData,
                const FloatOffsetType& itkNotUsed(offset) )
{
  // Global data structure
  // m_dx -> Intensity first derivative
  // m_dxy -> Intensity second derivative
  // m_DT_dxy -> Diffusion tensor first derivative
  GlobalDataStruct *gd = (GlobalDataStruct *)globalData;

  // Compute the first and 2nd derivative for the intensity images
  this->ComputeIntensityFirstAndSecondDerivatives( it, gd );

  // We are provided the diffusion tensor matrix first derivatives, so
  // copy them into the global data struct
  this->CopyDerivativeMatrixToGlobalData( tensorDerivative, gd );

  // Compute the update term
  return this->ComputeFinalUpdateTerm( gt, gd );
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::CopyDerivativeMatrixToGlobalData(
    const DerivativeMatrixImageRegionType &derivative,
     GlobalDataStruct *gd) const
{
  assert( gd );
  DerivativeMatrixType derivativePixel = derivative.Get();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      gd->m_DT_dxy[i][j] = derivativePixel(i, j);
      }
    }
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::DerivativeMatrixType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeDiffusionFirstDerivative(const DiffusionTensorNeighborhoodType &gt,
                                  GlobalDataStruct *gd) const
{
  DerivativeMatrixType derivativePixel;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    DiffusionTensorType positionA_Tensor_value = gt.GetPixel( m_positionA[i] );
    DiffusionTensorType positionB_Tensor_value = gt.GetPixel( m_positionB[i] );

    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      derivativePixel(i,j) = 0.5 *
                ( positionA_Tensor_value(i,j) - positionB_Tensor_value(i,j) );
      if( gd )
        {
        gd->m_DT_dxy[i][j] = derivativePixel(i,j);
        }
      }
    }
  return derivativePixel;
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeDiffusionFirstDerivative(const DiffusionTensorNeighborhoodType &gt,
                                  DerivativeMatrixImageRegionType &output ) const
{
  output.Set( this->ComputeDiffusionFirstDerivative( gt, NULL ) );
}

template< class TImageType >
void
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeIntensityFirstAndSecondDerivatives(const NeighborhoodType &it,
                                            GlobalDataStruct *gd) const
{
  assert( gd );
  gd->m_GradMagSqr = 1.0e-6;
  const ScalarValueType center_value  = it.GetCenterPixel();
  for( unsigned int i = 0; i < ImageDimension; i++)
    {
    const ScalarValueType it_positionA = it.GetPixel( m_positionA[i] );
    const ScalarValueType it_positionB = it.GetPixel( m_positionB[i] );

    gd->m_dx[i] = 0.5 * ( it_positionA - it_positionB );

    gd->m_dxy[i][i] = it_positionA + it_positionB - 2.0 * center_value;

    for( unsigned int j = i+1; j < ImageDimension; j++ )
      {
      gd->m_dxy[i][j] = gd->m_dxy[j][i] = 0.25 *
                                          ( it.GetPixel( m_positionAa[i][j] )
                                          - it.GetPixel( m_positionBa[i][j] )
                                          - it.GetPixel( m_positionCa[i][j] )
                                          + it.GetPixel( m_positionDa[i][j] )
        );
      }
    }
}

template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeFinalUpdateTerm(const DiffusionTensorNeighborhoodType &gt,
                         const GlobalDataStruct *gd) const
{
  DiffusionTensorType center_Tensor_value = gt.GetCenterPixel();

  ScalarValueType   pdWrtDiffusion1;

  pdWrtDiffusion1 = gd->m_DT_dxy[0][0] * gd->m_dx[0]
                    + gd->m_DT_dxy[0][1] * gd->m_dx[1]
                    + gd->m_DT_dxy[0][2] * gd->m_dx[2];

  ScalarValueType   pdWrtDiffusion2;

  pdWrtDiffusion2 = gd->m_DT_dxy[1][0] * gd->m_dx[0]
                    + gd->m_DT_dxy[1][1] * gd->m_dx[1]
                    + gd->m_DT_dxy[1][2] * gd->m_dx[2];

  ScalarValueType  pdWrtDiffusion3;

  pdWrtDiffusion3 = gd->m_DT_dxy[2][0] * gd->m_dx[0]
                    + gd->m_DT_dxy[2][1] * gd->m_dx[1]
                    + gd->m_DT_dxy[2][2] * gd->m_dx[2];

  ScalarValueType   pdWrtImageIntensity1;

  pdWrtImageIntensity1 = center_Tensor_value(0,0) * gd->m_dxy[0][0]  +
                    center_Tensor_value(0,1) * gd->m_dxy[0][1] +
                    center_Tensor_value(0,2) * gd->m_dxy[0][2];

  ScalarValueType   pdWrtImageIntensity2;

  pdWrtImageIntensity2 = center_Tensor_value(1,0) * gd->m_dxy[1][0]  +
                    center_Tensor_value(1,1) * gd->m_dxy[1][1] +
                    center_Tensor_value(1,2) * gd->m_dxy[1][2];

  ScalarValueType   pdWrtImageIntensity3;

  pdWrtImageIntensity3 = center_Tensor_value(2,0) * gd->m_dxy[2][0]  +
                    center_Tensor_value(2,1) * gd->m_dxy[2][1] +
                    center_Tensor_value(2,2) * gd->m_dxy[2][2];

  ScalarValueType   total;

  total = pdWrtDiffusion1 + pdWrtDiffusion2 + pdWrtDiffusion3 +
         pdWrtImageIntensity1 + pdWrtImageIntensity2 + pdWrtImageIntensity3;

  return ( PixelType ) ( total );
}

template <class TImageType>
void
AnisotropicDiffusionTensorFunction<TImageType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TImageType>
template < class TPixel, unsigned int VImageDimension >
void
AnisotropicDiffusionTensorFunction<TImageType>
::CheckTimeStepStability( const itk::Image< TPixel, VImageDimension > *
  input, bool useImageSpacing )
{
  double minSpacing;
  if( useImageSpacing )
    {
    minSpacing = input->GetSpacing()[0];
    for( unsigned int i = 1; i < ImageDimension; i++ )
      {
      if( input->GetSpacing()[i] < minSpacing)
        {
        minSpacing = input->GetSpacing()[i];
        }
      }
    }
  else
    {
    minSpacing = 1.0;
    }

  double ratio = minSpacing / vcl_pow(2.0,
    static_cast<double>(ImageDimension) + 1); // plus 1?

  if( m_TimeStep > ratio )
    {
    itkWarningMacro(<< std::endl
      << "Anisotropic diffusion unstable time step:"
      << m_TimeStep << std::endl << "Minimum stable time step"
      << "for this image is " << ratio );
    }
}

} // end namespace itk

#endif
