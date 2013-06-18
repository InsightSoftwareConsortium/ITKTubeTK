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

#ifndef __itktubeObjectDocumentToObjectSource_hxx
#define __itktubeObjectDocumentToObjectSource_hxx

#include "itktubeObjectDocumentToObjectSource.h"

namespace itk
{

namespace tube
{

template< class TInputObjectDocument, unsigned int TDimension >
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::ObjectDocumentToObjectSource( void )
{
  m_ComposedTransformIsIdentity = true;
  this->ApplyTransforms( true );  // Set transforms to be built (default value)
  this->SetNumberOfRequiredInputs(1);
}


template< class TInputObjectDocument, unsigned int TDimension >
const typename ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>::DocumentType *
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::GetInput( void )
{
  if(this->GetNumberOfInputs() < 1)
    {
    return 0;
    }

  return static_cast<const DocumentType * > (this->ProcessObject::GetInput(0) );
}


template< class TInputObjectDocument, unsigned int TDimension >
typename ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>::TransformPointer
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::GetComposedTransform( void )
{
  ConstDocumentPointer doc = static_cast<const DocumentType * > (this->ProcessObject::GetInput(0) );
  if( m_ApplyTransforms )
    {
    return ComposeTransforms( doc, m_StartTransforms, m_EndTransforms );
    }
  else
    {
    return ComposeTransforms( doc, 0, -1 );
    }
}


template< class TInputObjectDocument, unsigned int TDimension >
void
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::ApplyTransforms( int start, int end )
{
  m_ApplyTransforms = true;
  m_StartTransforms = start;
  m_EndTransforms = end;
}


template< class TInputObjectDocument, unsigned int TDimension >
void
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::ApplyTransforms( bool b )
{
  if( b == true )
    {
    ApplyTransforms( 0, -1 );
    }
  else
    {
    m_ApplyTransforms = b;
    }
}


/** Protected method is assumed to receive valid inputs
 *  (ie. startIndex !> endIndex || none < -1
 *  startIndex included upto, but excluding endIndex.
 *  eg. all transforms : startIndex = 0, endIndex = NumberOfTransforms
 *  -1 refers to the last element.
 *  Therefore startIndex = -1 is just the last element */
template< class TInputObjectDocument, unsigned int TDimension >
typename ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>::TransformPointer
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::ComposeTransforms( ConstDocumentPointer doc, int startIndex, int endIndex ) const
{
  typename DocumentType::TransformNameListType transNames = doc->GetTransformNames();
  typename DocumentType::TransformNameListType::const_iterator it_trans = transNames.begin();

  TransformPointer transform = TransformType::New();
  transform->SetIdentity();
  m_ComposedTransformIsIdentity = true;
  if( startIndex == -1 ) { startIndex = transNames.size() - 1; }
  if( endIndex == -1 ) { endIndex = transNames.size(); }

  int i = 0;
  // Skip the starting transforms
  while( it_trans != transNames.end() && i < startIndex  )
    {
    ++it_trans;
    i++;
    }

  //Compose the transform range
  while( it_trans != transNames.end() && i < endIndex  )
    {
    m_ComposedTransformIsIdentity = false;
    transform->Compose( ReadTransform( (*it_trans) ) );
    i++;
    ++it_trans;
    }
  return transform;
}


template< class TInputObjectDocument, unsigned int TDimension >
typename ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>::TransformPointer
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
  ::ReadTransform( const std::string & file ) const
{
  typename TransformReaderType::Pointer reader = TransformReaderType::New();
  reader->SetFileName( file );
  try
    {
    reader->Update();
    }
  catch( ... )
    {
    std::cerr << "Error:: No readable Transform found " << std::endl;
    return NULL;
    }

  typename TransformReaderType::GroupPointer grp = reader->GetGroup();
  if(!grp->GetNumberOfChildren())
    {
    return grp->GetObjectToParentTransform();
    }
  else
    {
    return (* (grp->GetChildren()->begin()))->GetObjectToParentTransform();
    }
}


template< class TInputObjectDocument, unsigned int TDimension >
typename ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>::DataObjectPointer
ObjectDocumentToObjectSource<TInputObjectDocument,TDimension>
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) )
{
  std::cout << "Making output" << std::endl;
  return static_cast<DataObject*>(OutputType::New().GetPointer());
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeObjectDocumentToObjectSource_hxx)
