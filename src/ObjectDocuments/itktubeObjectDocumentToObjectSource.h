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

#ifndef __itktubeObjectDocumentToObjectSource_h
#define __itktubeObjectDocumentToObjectSource_h

#include "itktubeObjectDocument.h"

#include <itkProcessObject.h>
#include <itkSpatialObjectReader.h>

namespace itk
{

namespace tube
{

/**
 * Filter that converts an object document to an object by reading and
 * composing all transforms from an object document for a single object.
 *
 * \note  Does not hold a buffer of objects read.
 * \ingroup  ObjectDocuments
 */
template< class TObjectDocument, unsigned int VDimension = 3 >
class ObjectDocumentToObjectSource : public ProcessObject
{
public:

  typedef ObjectDocumentToObjectSource         Self;
  typedef ProcessObject                        Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  typedef typename DataObject::Pointer         DataObjectPointer;

  typedef TObjectDocument                      DocumentType;
  typedef typename DocumentType::ConstPointer  ConstDocumentPointer;

  typedef typename SpatialObject< VDimension >::TransformType
                                               TransformType;
  typedef typename TransformType::Pointer      TransformPointer;

  itkNewMacro( Self );
  itkTypeMacro( ObjectDocumentToObjectSource, ProcessObject );

  /** Return whether the transforms should be applied. */
  itkGetMacro( ApplyTransforms, bool );

  /** Return the composed transform. */
  virtual TransformPointer GetComposedTransform( void );

  /** Return whether the composed transform is the identity transform. */
  itkGetMacro( ComposedTransformIsIdentity, bool );

  /** Return the input. */
  virtual const DocumentType * GetInput( void );

  /** Return the output. */
  virtual DataObject * GetOutput( void );

  /** Set whether the transforms should be applied. */
  virtual void SetApplyTransforms( bool applyTransforms );
  itkBooleanMacro( ApplyTransforms );

  /** Set whether the transforms should be applied. */
  virtual void SetApplyTransforms( int start, int end );

  /** Set the input. */
  virtual void SetInput( const DocumentType * input );

  /* Graft the specified data object onto the specified indexed output, but note
     that this function should be implemented by derived classes. */
  virtual void GraftNthOutput( unsigned int itkNotUsed( index ),
                               DataObject * itkNotUsed( data ) )
    {
    }

  /* Graft the specified data object onto the output, but note that this
     function should be implemented by derived classes. */
  virtual void GraftOutput( DataObject * itkNotUsed( data ) )
    {
    }

  /* Make an object to be used as the specified indexed output. */
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput( DataObjectPointerArraySizeType index )
    override;

protected:

  typedef SpatialObject<>                    OutputType;
  typedef SpatialObjectReader< VDimension >  TransformReaderType;

  /** Constructor. */
  ObjectDocumentToObjectSource( void );

  /** Destructor. */
  virtual ~ObjectDocumentToObjectSource( void );

  /** Return the end index for the list of transforms. */
  itkGetMacro( EndTransforms, int );

  /** Return the specified indexed output. */
  virtual DataObject * GetOutput( unsigned int index );

  /** Return the end index for the list of transforms. */
  itkGetMacro( StartTransforms, int );

  /** Set whether the composed transform is the identity transform. */
  itkSetMacro( ComposedTransformIsIdentity, bool );
  itkBooleanMacro( ComposedTransformIsIdentity );

  /** Set the end index for the list of transforms. */
  itkSetMacro( EndTransforms, int );

  /** Set the start index for the list of transforms. */
  itkSetMacro( StartTransforms, int );

  /** Generate the output data, but note that this function should be
      implemented by derived classes. */
  virtual void GenerateData( void ) override
    {
    }

  /** Generate the information describing the output data, but note that this
      function should be implemented by derived classes. */
  virtual void GenerateOutputInformation( void ) override
    {
    }

  /** Compose the transforms of the object ranging from start to end and return
      the final transform. Note that 0 is the first transform and -1 is the last
      transform. */
  virtual TransformPointer ComposeTransforms( ConstDocumentPointer document,
                                              int startIndex = 0,
                                              int endIndex = -1 ) const;

  /** Read the transform from the specified file. */
  virtual TransformPointer ReadTransform( const std::string & file ) const;

  /** Print information about the object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Copy constructor not implemented.
  ObjectDocumentToObjectSource( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * )
    override {};

  ConstDocumentPointer  m_Input;
  int                   m_StartTransforms;
  int                   m_EndTransforms;
  mutable bool          m_ComposedTransformIsIdentity;
  bool                  m_ApplyTransforms;

}; // End class ObjectDocumentToObjectSource

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeObjectDocumentToObjectSource.hxx"
#endif

#endif // End !defined( __itktubeObjectDocumentToObjectSource_h )
