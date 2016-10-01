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
#ifndef __tubeCropImage_hxx
#define __tubeCropImage_hxx

#include "tubeCropImage.h"


namespace tube
{
template< typename TInputImage, typename TOutputImage >
CropImage< TInputImage, TOutputImage >
::CropImage( void )
{
  m_CropFilter = CropFilterType::New();
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetMin( typename ImageType::IndexType roiMin )
{
  m_CropFilter->SetMin( roiMin );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetMax( typename ImageType::IndexType roiMax )
{
  m_CropFilter->SetMax( roiMax );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetSize( typename ImageType::SizeType roiSize )
{
  m_CropFilter->SetSize( roiSize );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetCenter( typename ImageType::IndexType roiCenter )
{
  m_CropFilter->SetCenter( roiCenter );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetBoundary( typename ImageType::IndexType roiBoundary )
{
  m_CropFilter->SetBoundary( roiBoundary );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetMatchVolume( typename ImageType::ConstPointer matchVolume )
{
  m_CropFilter->SetMatchVolume( matchVolume );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetMatchMask( typename ImageType::Pointer maskImage )
{
  m_CropFilter->SetMatchMask( maskImage );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetSplitInput( typename ImageType::IndexType splitIndex,
  typename ImageType::IndexType roiIndex )
{
  m_CropFilter->SetSplitInput( splitIndex, roiIndex );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::SetInput( const TInputImage *inputImage )
{
  m_CropFilter->SetInput( inputImage );
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::Update()
{
  m_CropFilter->Update();
}

template< typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
CropImage< TInputImage, TOutputImage >
::GetOutput()
{
  return m_CropFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
CropImage< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << m_CropFilter << std::endl;
}

} // end namespace tube


#endif
