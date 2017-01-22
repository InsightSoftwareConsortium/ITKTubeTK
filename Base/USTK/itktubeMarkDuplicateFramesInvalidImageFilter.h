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

#ifndef __itktubeMarkDuplicateFramesInvalidImageFilter_h
#define __itktubeMarkDuplicateFramesInvalidImageFilter_h

#include <itkProcessObject.h>
#include <itkDomainThreader.h>
#include <itkSimpleDataObjectDecorator.h>
#include <itkThreadedImageRegionPartitioner.h>

namespace itk
{

namespace tube
{

/** \class MarkDuplicateFramesInvalidImageFilterThreader
 *
 * \brief Threader for MarkDuplicateFramesInvalidImageFilter.
 */
template< typename TAssociate >
class MarkDuplicateFramesInvalidImageFilterThreader
  : public DomainThreader< ThreadedImageRegionPartitioner< 3 >, TAssociate >
{
public:
  const static unsigned int Dimension = 3;
  typedef ThreadedImageRegionPartitioner< Dimension > PartitionerType;

  /** Standard class typedefs. */
  typedef MarkDuplicateFramesInvalidImageFilterThreader Self;
  typedef DomainThreader< PartitionerType, TAssociate > Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  typedef TAssociate                      AssociateType;
  typedef typename Superclass::DomainType DomainType;

  itkNewMacro( Self );

protected:
  MarkDuplicateFramesInvalidImageFilterThreader() {}

private:
  MarkDuplicateFramesInvalidImageFilterThreader( const Self & );
  void operator=( const Self & ); // purposely not implemented

  virtual void BeforeThreadedExecution( void );
  virtual void ThreadedExecution( const DomainType & subDomain,
    const ThreadIdType threadId );
  virtual void AfterThreadedExecution( void );

  typedef std::list< SizeValueType >       InvalidFramesType;
  typedef std::vector< InvalidFramesType > InvalidFramesPerThreadType;
  InvalidFramesPerThreadType m_InvalidFramesPerThread;

};


/** \class MarkDuplicateFramesInvalidImageFilter
 *
 * \brief Detect duplicate frames and mark invalid in the Plus metadata.
 *
 * The screengrabber's frame rate is often higher that the ultrasounds's frame
 * rate, so detect duplicate frames ( or partial duplicates with an incomplete
 * refresh ) by subtracting subsequent frames.
 *
 * The input to this filter should be a scalar image.  It runs in place on the
 * content of the image, and leaves it unchanged.  The output image, though,
 * has a new MetaDataDictionary with modified values from the dictionary
 * specified with SetInputMetaDataDictionary.
 *
 * Pixels are considered the same if their difference is within the Tolerance
 * parameter.  Frames are considered the same if same-valued pixels exceed the
 * FractionalThreshold parameter.
 *
 */
template< typename TInputImage >
class MarkDuplicateFramesInvalidImageFilter
  : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MarkDuplicateFramesInvalidImageFilter Self;
  typedef ProcessObject                         Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( MarkDuplicateFramesInvalidImageFilter, ProcessObject );

  /** Some convenient typedefs. */
  typedef TInputImage                         InputImageType;
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename InputImageType::PixelType  InputImagePixelType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    InputImageType::ImageDimension );

  typedef SimpleDataObjectDecorator< MetaDataDictionary >
    DecoratedMetaDataDictionaryType;

  /** Set the tolerance which determines if a pixel is unchanged.  Noise with
   * analog framegrabbers requires this to be non-zero. This is considered
   * inclusive. */
  itkSetMacro( Tolerance, InputImagePixelType );
  itkGetConstMacro( Tolerance, InputImagePixelType );

  /** Set the fractional threshold after which a frame is considered a
   * duplicate. */
  itkSetMacro( FractionalThreshold, double );
  itkGetConstMacro( FractionalThreshold, double );

  virtual void SetInput( const InputImageType * image );
  virtual const InputImageType * GetInput( void ) const;

  /** Set/Get the Plus ultrasound format MetaDataDictionary that is modified
   * and added to the output image. */
  void SetInputMetaDataDictionary( const MetaDataDictionary * dictionary );
  const MetaDataDictionary * GetInputMetaDataDictionary( void ) const;

  const MetaDataDictionary & GetOutputMetaDataDictionary( void ) const;

protected:
  MarkDuplicateFramesInvalidImageFilter( void );
  virtual ~MarkDuplicateFramesInvalidImageFilter( void );

  using Superclass::MakeOutput;
  virtual DataObject::Pointer MakeOutput( DataObjectPointerArraySizeType index );

  virtual void GenerateData( void );

private:
  MarkDuplicateFramesInvalidImageFilter( const Self & );
  void operator=( const Self & ); // purposely not implemented

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  InputImagePixelType m_Tolerance;
  double m_FractionalThreshold;

  const MetaDataDictionary * m_InputMetaDataDictionary;

  friend class MarkDuplicateFramesInvalidImageFilterThreader< Self >;
  typedef MarkDuplicateFramesInvalidImageFilterThreader< Self > ThreaderType;
  typename ThreaderType::Pointer m_Threader;

}; // End class MarkDuplicateFramesInvalidImageFilter

} // End namespace tubetk

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeMarkDuplicateFramesInvalidImageFilter.hxx"
#endif

#endif // End !defined( __itktubeMarkDuplicateFramesInvalidImageFilter_h )
