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

#ifndef __itktubeObjectDocumentToObjectSource_hxx
#define __itktubeObjectDocumentToObjectSource_hxx


namespace itk
{

namespace tube
{

template< class TObjectDocument, unsigned int VDimension >
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::ObjectDocumentToObjectSource( void )
{
  m_ComposedTransformIsIdentity = true;
  m_StartTransforms = 0;
  m_EndTransforms = -1;
  m_ApplyTransforms = true;

  this->SetApplyTransforms( m_ApplyTransforms );
  this->SetNumberOfRequiredInputs( 1 );
}

template< class TObjectDocument, unsigned int VDimension >
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::~ObjectDocumentToObjectSource( void )
{
}

template< class TObjectDocument, unsigned int VDimension >
typename ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::TransformPointer
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::GetComposedTransform( void )
{
  ConstDocumentPointer document
    = static_cast< const DocumentType * >( this->Superclass::GetInput( 0 ) );

  if( m_ApplyTransforms )
    {
    return this->ComposeTransforms( document, m_StartTransforms,
      m_EndTransforms );
    }

  return this->ComposeTransforms( document, 0, -1 );
}

template< class TObjectDocument, unsigned int VDimension >
const typename ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::DocumentType *
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::GetInput( void )
{
  if( this->GetNumberOfInputs() < 1 )
    {
    return NULL;
    }

  return static_cast< const DocumentType * >( this->Superclass::GetInput( 0 ) );
}

template< class TObjectDocument, unsigned int VDimension >
DataObject *
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::GetOutput( void )
{
  return this->Superclass::GetOutput( 0 );
}

template< class TObjectDocument, unsigned int VDimension >
DataObject *
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::GetOutput( unsigned int index )
{
  return this->Superclass::GetOutput( index );
}

template< class TObjectDocument, unsigned int VDimension >
void
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::SetApplyTransforms( bool applyTransforms )
{
  if( applyTransforms )
    {
    this->SetApplyTransforms( 0, -1 );
    }
  else
    {
    m_ApplyTransforms = applyTransforms;
    }
}

template< class TObjectDocument, unsigned int VDimension >
void
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::SetApplyTransforms( int start, int end )
{
  m_ApplyTransforms = true;
  m_StartTransforms = start;
  m_EndTransforms = end;
}

template< class TObjectDocument, unsigned int VDimension >
void
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::SetInput( const DocumentType * input )
{
  this->Superclass::SetNthInput( 0, const_cast< DocumentType * >( input ) );
}

template< class TObjectDocument, unsigned int VDimension >
typename ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::DataObjectPointer
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed( index ) )
{
  return static_cast< DataObject * >( OutputType::New().GetPointer() );
}

/* Function is assumed to receive valid inputs
   i.e., startIndex !> endIndex || none < -1
   startIndex included up to, but excluding, endIndex.
   e.g., all transforms: startIndex = 0, endIndex = NumberOfTransforms
   -1 refers to the last element.
   Therefore startIndex = -1 is just the last element. */
template< class TObjectDocument, unsigned int VDimension >
typename ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::TransformPointer
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::ComposeTransforms( ConstDocumentPointer document, int startIndex,
                     int endIndex ) const
{
  typedef typename DocumentType::TransformNameListType  TransformNameListType;

  TransformNameListType transformNames = document->GetTransformNames();
  typename TransformNameListType::const_iterator iter = transformNames.begin();

  TransformPointer transform = TransformType::New();
  transform->SetIdentity();
  m_ComposedTransformIsIdentity = true;

  if( startIndex == -1 )
    {
    startIndex = transformNames.size() - 1;
    }

  if( endIndex == -1 )
    {
    endIndex = transformNames.size();
    }

  int i = 0;

  // Skip the starting transforms.
  while( iter != transformNames.end() && i < startIndex )
    {
    ++i;
    ++iter;
    }

  // Compose the transform range.
  while( iter != transformNames.end() && i < endIndex )
    {
    m_ComposedTransformIsIdentity = false;
    transform->Compose( this->ReadTransform( *iter ) );

    ++i;
    ++iter;
    }

  return transform;
}

template< class TObjectDocument, unsigned int VDimension >
typename ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::TransformPointer
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
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
    itkWarningMacro( << "No readable transform found." );
    return (TransformType *)NULL;
    }

  typename TransformReaderType::GroupPointer group = reader->GetGroup();

  if( !group->GetNumberOfChildren() )
    {
    return group->GetObjectToParentTransform();
    }
  else
    {
    return ( * ( group->GetChildren()->begin() ) )->
      GetObjectToParentTransform();
    }
}

template< class TObjectDocument, unsigned int VDimension >
void
ObjectDocumentToObjectSource< TObjectDocument, VDimension >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_Input )
    {
    os << indent << "Input:                       " << m_Input << std::endl;
    }

  os << indent << "StartTransforms:             " << m_StartTransforms
     << std::endl;
  os << indent << "EndTransforms:               " << m_EndTransforms
     << std::endl;
  os << indent << "ComposedTransformIsIdentity: "
     << m_ComposedTransformIsIdentity << std::endl;
  os << indent << "ApplyTransforms:             " << m_ApplyTransforms
     << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeObjectDocumentToObjectSource_hxx )
