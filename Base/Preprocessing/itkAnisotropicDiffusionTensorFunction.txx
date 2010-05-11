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
  RadiusType r;

  for( unsigned int i=0; i < ImageDimension; i++ )
    {
    r[i] = 1;
    } 

  this->SetRadius(r);
  
  // Dummy neighborhood.
  NeighborhoodType it;
  it.SetRadius( r );
  
  // Find the center index of the neighborhood.
  m_Center =  it.Size() / 2;

  // Get the stride length for each axis.
  for(unsigned int i = 0; i < ImageDimension; i++)
    {  m_xStride[i] = it.GetStride(i); }
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
  return this->ComputeUpdate( it, diffusionTensor,globalData, offset ); 
}
template< class TImageType >
typename AnisotropicDiffusionTensorFunction< TImageType >::PixelType
AnisotropicDiffusionTensorFunction< TImageType >
::ComputeUpdate(const NeighborhoodType &it, 
                const DiffusionTensorNeighborhoodType &gt,
                void *globalData,
                const FloatOffsetType& offset)
{
  const ScalarValueType center_value  = it.GetCenterPixel();

  FloatOffsetType offsetCum = offset;

  // Global data structure
  GlobalDataStruct *gd = (GlobalDataStruct *)globalData;

  // m_dx -> Intensity first derivative 
  // m_dxy -> Intensity second derivative
  // m_DT_dxy -> Diffusion tensor first derivative

  // Compute the first and 2nd derivative 
  gd->m_GradMagSqr = 1.0e-6;
  for( unsigned int i = 0; i < ImageDimension; i++)
    {
    const unsigned int positionA = 
      static_cast<unsigned int>( m_Center + m_xStride[i]);
    const unsigned int positionB = 
      static_cast<unsigned int>( m_Center - m_xStride[i]);

    gd->m_dx[i] = 0.5 * (it.GetPixel( positionA ) - 
                     it.GetPixel( positionB )    );

    gd->m_dxy[i][i] = it.GetPixel( positionA )
      + it.GetPixel( positionB ) - 2.0 * center_value;
    
    for( unsigned int j = i+1; j < ImageDimension; j++ )
      {
      const unsigned int positionAa = static_cast<unsigned int>( 
        m_Center - m_xStride[i] - m_xStride[j] );
      const unsigned int positionBa = static_cast<unsigned int>( 
        m_Center - m_xStride[i] + m_xStride[j] );
      const unsigned int positionCa = static_cast<unsigned int>( 
        m_Center + m_xStride[i] - m_xStride[j] );
      const unsigned int positionDa = static_cast<unsigned int>( 
        m_Center + m_xStride[i] + m_xStride[j] );

      gd->m_dxy[i][j] = gd->m_dxy[j][i] = 0.25 *( it.GetPixel( positionAa )
                                          - it.GetPixel( positionBa )
                                          - it.GetPixel( positionCa )
                                          + it.GetPixel( positionDa )
        );
      }
    }

  // Compute the diffusion tensor matrix first derivatives 
  TensorPixelType center_Tensor_value  = gt.GetCenterPixel();

  for( unsigned i = 0; i < ImageDimension; i++)
    {
    const unsigned int positionA = 
      static_cast<unsigned int>( m_Center + m_xStride[i]);
    const unsigned int positionB = 
      static_cast<unsigned int>( m_Center - m_xStride[i]);
    
    TensorPixelType positionA_Tensor_value = gt.GetPixel( positionA );
    TensorPixelType positionB_Tensor_value = gt.GetPixel( positionB );

    for( unsigned int j = 0; j < ImageDimension; j++)
      { 
      gd->m_DT_dxy[i][j] = 0.5 *  ( positionA_Tensor_value(i,j) - 
                                positionB_Tensor_value(i,j) ); 
      }
    }

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

  pdWrtImageIntensity1 = center_Tensor_value(0,0) *  gd->m_dxy[0][0]  + 
                    center_Tensor_value(0,1) *  gd->m_dxy[0][1] +
                    center_Tensor_value(0,2) *  gd->m_dxy[0][2];
  
  ScalarValueType   pdWrtImageIntensity2;

  pdWrtImageIntensity2 = center_Tensor_value(1,0) *  gd->m_dxy[1][0]  + 
                    center_Tensor_value(1,1) *  gd->m_dxy[1][1] +
                    center_Tensor_value(1,2) *  gd->m_dxy[1][2];
 
  ScalarValueType   pdWrtImageIntensity3;

  pdWrtImageIntensity3 = center_Tensor_value(2,0) *  gd->m_dxy[2][0]  + 
                    center_Tensor_value(2,1) *  gd->m_dxy[2][1] +
                    center_Tensor_value(2,2) *  gd->m_dxy[2][2];
 
  ScalarValueType   total;

  total = pdWrtDiffusion1 + pdWrtDiffusion2 + pdWrtDiffusion3 +
         pdWrtImageIntensity1 + pdWrtImageIntensity2 + pdWrtImageIntensity3;

  return ( PixelType ) ( total );
} 

template <class TImageType>
void
AnisotropicDiffusionTensorFunction<TImageType>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk

#endif
