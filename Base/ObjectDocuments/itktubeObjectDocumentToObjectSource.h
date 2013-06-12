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

#ifndef __itktubeObjectDocumentToObjectSource_h
#define __itktubeObjectDocumentToObjectSource_h

#include "itktubeObjectDocument.h"

#include <itkProcessObject.h>
#include <itkSpatialObjectReader.h>

namespace itk
{

namespace tube
{

/** \class ObjectDocumentToObjectSource
 * \brief Base Class for all filters that convert an ObjectDocument to the referenced object
 *
 * This base class handles transform composition and use for all of the subclasses.
 * Reads and composes all transforms from ObjectDocument for single object.
 */
template< class TInputObjectDocument, unsigned int TDimension = 3 >
class ITK_EXPORT ObjectDocumentToObjectSource : public ProcessObject
{
public:

  typedef ObjectDocumentToObjectSource                      Self;
  typedef ProcessObject                                     Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  typedef TInputObjectDocument                              DocumentType;
  typedef typename DocumentType::Pointer                    DocumentPointer;
  typedef typename DocumentType::ConstPointer               ConstDocumentPointer;

  typedef typename SpatialObject<TDimension>::TransformType TransformType;
  typedef typename TransformType::Pointer                   TransformPointer;


protected:

  typedef SpatialObject<> OutputType;

public:

  /** Smart Pointer type to a DataObject. */
  typedef typename DataObject::Pointer DataObjectPointer;

  itkNewMacro( Self );
  itkTypeMacro( Self, Superclass );

  /** Sets the transforms that will be composed and applied to the object are applied to it */
  void ApplyTransforms( int start, int end );

  /* Sets all the transforms to be applied (equivalent to ->ApplyTransforms(0, -1) ) */
  void ApplyTransforms( bool );
  bool ApplyTransforms( void )
    {
    return m_ApplyTransforms;
    }

  /** Return the composed Transform combining start to end */
  TransformPointer GetComposedTransform( void );

  bool ComposedTransformIsIdentity( void )
    {
    return m_ComposedTransformIsIdentity;
    }

  /** Purposely not implemented -- To be implemented by deriving class */
  virtual void GenerateData( void )
    {
    std::cout << " In source filter class" << std::endl;
    }

  /** Purposely not implemented -- To be implemented by deriving class */
  virtual void GenerateOutputInformation( void ) {} // do nothing

  /** Get the output object */
  DataObject * GetOutput( void )
    {
    return this->ProcessObject::GetOutput(0);
    }

  /** Set the input document object */
  using Superclass::SetInput;
  virtual void SetInput( const DocumentType * input )
    {
    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(0, const_cast< DocumentType * >( input ) );
    }

  /** Get the input document object */
  virtual const DocumentType * GetInput( void );

  /** Must be be defined by deriving class */
  virtual void GraftOutput( DataObject * ) {}
  virtual void GraftNthOutput( unsigned int, DataObject * ) {}

  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

protected:

  /** Typedef for transform file reader */
  typedef SpatialObjectReader<TDimension>                          TransformReaderType;

  ObjectDocumentToObjectSource( void );
  ~ObjectDocumentToObjectSource( void ) {}

  DataObject * GetOutput(unsigned int idx) { return this->ProcessObject::GetOutput(idx); }

  /**
   * Composes object's transforms ranging from start to end
   * and returns the final transform
   *
   * NOTE: 0 = first transform, -1 is last transform
   */
  TransformPointer ComposeTransforms( ConstDocumentPointer doc, int startIndex=0, int endIndex=-1 ) const;

  /** Read the transform from file */
  TransformPointer ReadTransform( const std::string & file ) const;

  /** Flag to determine whether to apply transforms or not */
  bool                                    m_ApplyTransforms;

  /** Selected range of transforms to be applied  */
  int                                     m_StartTransforms;
  int                                     m_EndTransforms;
  mutable bool                            m_ComposedTransformIsIdentity;

private:

  ObjectDocumentToObjectSource(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  ConstDocumentPointer                              m_Input;

}; // End class ObjectDocumentToObjectSource

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeObjectDocumentToObjectSource.hxx"
#endif

#endif // End !defined(__itktubeObjectDocumentToObjectSource_h)
