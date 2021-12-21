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
#ifndef __tubeConvertShrunkenSeedImageToList_hxx
#define __tubeConvertShrunkenSeedImageToList_hxx



namespace tube
{
template< class TImage, class TPointsImage >
ConvertShrunkenSeedImageToList< TImage, TPointsImage >
::ConvertShrunkenSeedImageToList( void )
{
  m_ConvertShrunkenSeedImageToListFilter =
    ConvertShrunkenSeedImageToListFilterType::New();
}

template< class TImage, class TPointsImage >
typename itk::tube::ConvertShrunkenSeedImageToListFilter< TImage, TPointsImage >
::VnlMatrixType
ConvertShrunkenSeedImageToList< TImage, TPointsImage >
::GetOutput()
{
  return m_ConvertShrunkenSeedImageToListFilter->GetOutput()->Get();
}

template< class TImage, class TPointsImage >
void
ConvertShrunkenSeedImageToList< TImage, TPointsImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << "Filter = " << m_ConvertShrunkenSeedImageToListFilter << std::endl;
}

} // end namespace tube


#endif
