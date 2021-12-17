/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeSubSampleSpatialObjectFilter_hxx
#define __itktubeSubSampleSpatialObjectFilter_hxx


#include "itktubeSubSampleTubeSpatialObjectFilter.h"

#include <itkSpatialObjectFactory.h>

namespace itk
{

namespace tube
{

template< unsigned int ObjectDimension >
SubSampleSpatialObjectFilter< ObjectDimension >
::SubSampleSpatialObjectFilter( void )
  : m_Sampling( 1 )
{
  SpatialObjectFactoryBase::RegisterDefaultSpatialObjects();
  SpatialObjectFactory< SpatialObjectType >::RegisterSpatialObject();
  SpatialObjectFactory< TubeSpatialObjectType >::RegisterSpatialObject();
}

template< unsigned int ObjectDimension >
SubSampleSpatialObjectFilter< ObjectDimension >
::~SubSampleSpatialObjectFilter( void )
{
}

template< unsigned int ObjectDimension >
void
SubSampleSpatialObjectFilter< ObjectDimension >
::SubSampleLevel( const SpatialObjectType * input,
  typename SpatialObjectType::Pointer output, bool graftOutput )
{
  // We make the copy and sub-sample if it is a tube.
  const TubeSpatialObjectType * inputAsTube =
    dynamic_cast< const TubeSpatialObjectType * >( input );
  typename SpatialObjectType::Pointer newParent;
  if( inputAsTube != NULL )
    {
    typedef SubSampleTubeSpatialObjectFilter< ObjectDimension >
      SubSampleTubeFilterType;
    typename SubSampleTubeFilterType::Pointer subSampleTubeFilter
      = SubSampleTubeFilterType::New();

    subSampleTubeFilter->SetSampling( this->GetSampling() );
    subSampleTubeFilter->SetInput( inputAsTube );
    if( graftOutput )
      {
      TubeSpatialObjectType * outputAsTube =
        dynamic_cast< TubeSpatialObjectType * >( output.GetPointer() );
      subSampleTubeFilter->GraftOutput( outputAsTube );
      subSampleTubeFilter->Update();
      newParent = subSampleTubeFilter->GetOutput();
      this->GraftOutput( newParent );
      }
    else
      {
      subSampleTubeFilter->Update();
      newParent = subSampleTubeFilter->GetOutput();
      output->AddChild( newParent );
      }
    }
  else
    {
    if( graftOutput )
      {
      // Output = Input in baseclass
      newParent = dynamic_cast< SpatialObjectType * >( this->GetOutput() );
      }
    else
      {
      newParent = input->Clone();
      output->AddChild( newParent );
      }
    }

  typedef typename SpatialObjectType::ChildrenListType ChildrenListType;
  ChildrenListType *children = input->GetChildren();
  typename ChildrenListType::const_iterator it = children->begin();
  while( it != children->end() )
    {
    this->SubSampleLevel( *it, newParent );
    ++it;
    }
  delete children;
}


template< unsigned int ObjectDimension >
void
SubSampleSpatialObjectFilter< ObjectDimension >
::GenerateData( void )
{
  const SpatialObjectType * input = this->GetInput();

  typename SpatialObjectType::Pointer output = this->GetOutput();

  this->SubSampleLevel( input, output, true );

  this->GraftOutput(output);
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSubSampleSpatialObjectFilter_hxx )
