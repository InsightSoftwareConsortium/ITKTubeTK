/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeConvertSpatialGraphToImage_hxx
#define __tubeConvertSpatialGraphToImage_hxx

#include "tubeConvertSpatialGraphToImage.h"


namespace tube
{
template< typename TInputImage, typename TOutputImage >
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::ConvertSpatialGraphToImage( void )
{
  m_ConvertSpatialGraphToImageFilter =
    ConvertSpatialGraphToImageFilterType::New();
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::GetAdjacencyMatrixImage()
{
  return m_ConvertSpatialGraphToImageFilter->GetAdjacencyMatrixImage();
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::GetBranchnessImage()
{
  return m_ConvertSpatialGraphToImageFilter->GetBranchnessImage();
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::GetRadiusImage()
{
  return m_ConvertSpatialGraphToImageFilter->GetRadiusImage();
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::GetCentralityImage()
{
  return m_ConvertSpatialGraphToImageFilter->GetCentralityImage();
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::SetAdjacencyMatrix( vnl_matrix< double > a)
{
  m_ConvertSpatialGraphToImageFilter->SetAdjacencyMatrix( a );
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::SetBranchnessVector( vnl_vector< double > b )
{
  m_ConvertSpatialGraphToImageFilter->SetBranchnessVector( b );
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::SetRadiusVector( vnl_vector< double > r )
{
  m_ConvertSpatialGraphToImageFilter->SetRadiusVector( r );
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::SetCentralityVector( vnl_vector< double > c )
{
  m_ConvertSpatialGraphToImageFilter->SetCentralityVector( c );
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::SetInput( const TInputImage *inputImage )
{
  m_ConvertSpatialGraphToImageFilter->SetInput( inputImage );
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::Update()
{
  m_ConvertSpatialGraphToImageFilter->Update();
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::GetOutput()
{
  return m_ConvertSpatialGraphToImageFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
ConvertSpatialGraphToImage< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  m_ConvertSpatialGraphToImageFilter->PrintSelf( os, indent);
}

} // end namespace tube

#endif
