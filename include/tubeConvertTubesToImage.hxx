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
#ifndef __tubeConvertTubesToImage_hxx
#define __tubeConvertTubesToImage_hxx


namespace tube
{

template< class TImage >
ConvertTubesToImage< TImage >
::ConvertTubesToImage( void )
{
  m_Filter = FilterType::New();
  m_Filter->SetBuildRadiusImage( false );
  m_Filter->SetBuildTangentImage( false );

  m_TemplateImage = NULL;
}

template< class TImage >
void
ConvertTubesToImage< TImage >
::SetTemplateImage( const typename ConvertTubesToImage< TImage >::
  OutputImageType * pTemplateImage )
{
  if( this->m_TemplateImage != pTemplateImage )
    {
    this->m_TemplateImage = pTemplateImage;

    m_Filter->SetSize( pTemplateImage->GetLargestPossibleRegion().GetSize() );
    m_Filter->SetSpacing( pTemplateImage->GetSpacing() );
    m_Filter->SetDirection( pTemplateImage->GetDirection() );
    m_Filter->SetOrigin( pTemplateImage->GetOrigin() );

    this->Modified();
    }
}

template< class TImage >
void
ConvertTubesToImage< TImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Filter: " << m_Filter << std::endl;
}

}

#endif
