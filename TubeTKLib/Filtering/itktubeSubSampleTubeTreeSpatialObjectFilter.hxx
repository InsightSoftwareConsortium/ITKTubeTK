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

#ifndef __itktubeSubSampleTubeTreeSpatialObjectFilter_hxx
#define __itktubeSubSampleTubeTreeSpatialObjectFilter_hxx

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"

#include "itktubeSubSampleTubeSpatialObjectFilter.h"

#include <itkSpatialObjectFactory.h>

namespace itk
{

namespace tube
{

template< class TSpatialObject, class TTubeSpatialObject >
SubSampleTubeTreeSpatialObjectFilter< TSpatialObject, TTubeSpatialObject >
::SubSampleTubeTreeSpatialObjectFilter( void )
  : m_Sampling( 1 )
{
  //! \todo fix this system.
  SpatialObjectFactoryBase::RegisterDefaultSpatialObjects();
  SpatialObjectFactory< SpatialObjectType >::RegisterSpatialObject();
  SpatialObjectFactory< TubeSpatialObjectType >::RegisterSpatialObject();
}

template< class TSpatialObject, class TTubeSpatialObject >
SubSampleTubeTreeSpatialObjectFilter< TSpatialObject, TTubeSpatialObject >
::~SubSampleTubeTreeSpatialObjectFilter( void )
{
}

template< class TSpatialObject, class TTubeSpatialObject >
void
SubSampleTubeTreeSpatialObjectFilter< TSpatialObject, TTubeSpatialObject >
::SubSampleLevel( const SpatialObjectBaseType * input,
  SpatialObjectBaseType * output, bool graftOutput )
{
  // We make the copy and sub-sample if it is a tube.
  const TubeSpatialObjectType * inputAsTube =
    dynamic_cast< const TubeSpatialObjectType * >( input );
  typename SpatialObjectBaseType::Pointer newParent;
  if( inputAsTube != NULL )
    {
    typedef SubSampleTubeSpatialObjectFilter< TubeSpatialObjectType >
      SubSampleTubeFilterType;
    typename SubSampleTubeFilterType::Pointer subSampleTubeFilter
      = SubSampleTubeFilterType::New();
    subSampleTubeFilter->SetSampling( this->GetSampling() );
    subSampleTubeFilter->SetInput( const_cast< TubeSpatialObjectType * >(
        inputAsTube ) );
    subSampleTubeFilter->Update();
    if( graftOutput )
      {
      output->CopyInformation( subSampleTubeFilter->GetOutput() );
      }
    else
      {
      output->AddSpatialObject( subSampleTubeFilter->GetOutput() );
      }
    newParent = dynamic_cast< SpatialObjectBaseType * >( output );
    }
  else
    {
    if( graftOutput )
      {
      output->CopyInformation( input );
      newParent = dynamic_cast< SpatialObjectBaseType * >( output );
      }
    else
      {
      const std::string spatialObjectType = input->
        GetSpatialObjectTypeAsString();
      LightObject::Pointer newSpatialObject =
        ObjectFactoryBase::CreateInstance( spatialObjectType.c_str() );

      typename SpatialObjectBaseType::Pointer newSpatialObjectBase =
        dynamic_cast< SpatialObjectBaseType * >(
        newSpatialObject.GetPointer() );
      if( newSpatialObjectBase.IsNull() )
        {
        itkExceptionMacro( << "Could not create an instance of "
          << spatialObjectType << ". The usual cause of this error is not"
          << "registering the SpatialObject with SpatialFactory." );
        }

      // Correct for extra reference count from CreateInstance()
      newSpatialObjectBase->UnRegister();

      newSpatialObjectBase->CopyInformation( input );
      output->AddSpatialObject( newSpatialObjectBase );
      newParent = newSpatialObjectBase;
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


template< class TSpatialObject, class TTubeSpatialObject >
void
SubSampleTubeTreeSpatialObjectFilter< TSpatialObject, TTubeSpatialObject >
::GenerateData( void )
{
  const SpatialObjectType * input = this->GetInput();

  typename SpatialObjectBaseType::Pointer output = this->GetOutput();
  output->Clear();

  this->SubSampleLevel( input, output, true );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSubSampleTubeTreeSpatialObjectFilter_hxx )
