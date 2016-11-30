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
#ifndef __tubeConvertImagesToCSV_hxx
#define __tubeConvertImagesToCSV_hxx

#include "tubeConvertImagesToCSV.h"


namespace tube
{
template< class TInputImage, class TInputMask >
ConvertImagesToCSV< TInputImage, TInputMask >
::ConvertImagesToCSV( void )
{
  m_ConvertImagesToCSVFilter = ConvertImagesToCSVFilterType::New();
}

template< class TInputImage, class TInputMask >
typename itk::tube::ConvertImagesToCSVFilter< TInputImage, TInputMask >
::VnlMatrixType
ConvertImagesToCSV< TInputImage, TInputMask >
::GetOutput()
{
  return m_ConvertImagesToCSVFilter->GetOutput()->Get();
}

template< class TInputImage, class TInputMask >
void
ConvertImagesToCSV< TInputImage, TInputMask >
::AddImage(TInputImage* image)
{
  m_ConvertImagesToCSVFilter->AddImage( image );
  this->Modified();
}

template< class TInputImage, class TInputMask >
void
ConvertImagesToCSV< TInputImage, TInputMask >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << "Filter = " << m_ConvertImagesToCSVFilter << std::endl;
}

} // end namespace tube


#endif
