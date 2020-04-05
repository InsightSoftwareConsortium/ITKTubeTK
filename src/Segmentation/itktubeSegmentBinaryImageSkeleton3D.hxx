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
#ifndef __itktubeSegmentBinaryImageSkeleton3D_hxx
#define __itktubeSegmentBinaryImageSkeleton3D_hxx

// TubeTK includes
#include "itktubeSegmentBinaryImageSkeleton3D.h"

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TPixel >
SegmentBinaryImageSkeleton3D< TPixel >
::SegmentBinaryImageSkeleton3D()
{
  m_Radius = 0;
}

template< class TPixel >
void
SegmentBinaryImageSkeleton3D< TPixel >
::GenerateData( void )
{
  Superclass::GenerateData();

  if( m_Radius > 0 )
    {
    typename DilateType::Pointer dilator;

    SEType binaryBall;
    binaryBall.SetRadius( m_Radius );
    binaryBall.CreateStructuringElement();
    dilator = DilateType::New();
    dilator->SetInput( this->GetOutput() );
    dilator->SetForegroundValue( 1 );
    dilator->SetKernel( binaryBall );
    dilator->Update();
    }
}

template< class TPixel >
void
SegmentBinaryImageSkeleton3D< TPixel >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << "Structuring element radius: " << m_Radius << std::endl;

}

} // End namespace tube

} // End namespace itk

#endif
