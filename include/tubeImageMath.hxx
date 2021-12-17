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
#ifndef __tubeImageMath_hxx
#define __tubeImageMath_hxx


namespace tube
{

template< class TInputImage, class TOutputImage >
ImageMath< TInputImage, TOutputImage >
::ImageMath( void )
{
  m_HistogramBinMin = 0;
  m_HistogramBinSize = 0;
}

template< class TInputImage, class TOutputImage >
void
ImageMath< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

}

#endif
