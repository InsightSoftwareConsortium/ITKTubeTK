/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeSegmentConnectedComponentsUsingParzenPDFs_hxx
#define __tubeSegmentConnectedComponentsUsingParzenPDFs_hxx


namespace tube
{

template< class TImage, class TLabelMap >
SegmentConnectedComponentsUsingParzenPDFs< TImage, TLabelMap >
::SegmentConnectedComponentsUsingParzenPDFs( void )
{
  m_Filter = FilterType::New();

  m_FVGenerator = FeatureVectorGeneratorType::New();
}

template< class TImage, class TLabelMap >
void
SegmentConnectedComponentsUsingParzenPDFs< TImage, TLabelMap >
::SetFeatureImage( InputImageType * img )
{
  m_FVGenerator->SetInput( img );
}

template< class TImage, class TLabelMap >
void
SegmentConnectedComponentsUsingParzenPDFs< TImage, TLabelMap >
::AddFeatureImage( InputImageType * img )
{
  m_FVGenerator->AddInput( img );
}

template< class TImage, class TLabelMap >
void
SegmentConnectedComponentsUsingParzenPDFs< TImage, TLabelMap >
::Update( void )
{
  m_Filter->SetFeatureVectorGenerator( m_FVGenerator );

  m_Filter->Update();
}

template< class TImage, class TLabelMap >
void
SegmentConnectedComponentsUsingParzenPDFs< TImage, TLabelMap >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << m_Filter << std::endl;
}

}

#endif
