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

#ifndef __itktubeVectorImageToListGenerator_h
#define __itktubeVectorImageToListGenerator_h

#include <itkDataObject.h>
#include <itkDataObjectDecorator.h>
#include <itkFixedArray.h>
#include <itkListSample.h>
#include <itkPixelTraits.h>
#include <itkProcessObject.h>

namespace itk
{

namespace tube
{

namespace Statistics
{

/** \class VectorImageToListGenerator
 *  \brief The class takes an image as input and generates a list sample as
 *  output.
 *
 *  There are differences between this class and VectorImageToListAdaptor.
 *  This class is not an adaptor. It creates a new list sample and does not
 *  provide a pseudo interface to the actual image to make it look like a
 *  list sample.
 *
 *  The class optionally allows you to specify a mask image as an input. The
 *  list sample ( if a mask is specified ) is constructed from pixels that are
 *  within the mask
 *
 * \todo
 * In future allow the filter to take a Spatial object as input so a
 * generic spatial object like an ellipse etc can be used as a mask.
 * Sure the ImageMaskSpatialObject
 * can represent image masks too, so why not make SpatialObjects the
 * default. I think the ImageMaskSpatialObject is slow in terms of
 * inefficient iteration through the image.
 *
 * \sa VectorImageToListAdaptor
 */
template< class TImage, class TMaskImage >
class VectorImageToListGenerator : public ProcessObject
{
public:
  /** Standard class typedefs */
  typedef VectorImageToListGenerator        Self;
  typedef ProcessObject                     Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( VectorImageToListGenerator, ProcessObject );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Image typedefs */
  typedef TImage                           ImageType;
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename ImageType::PixelType    PixelType;
  typedef PixelType                        MeasurementVectorType;

  /** Mask Image typedefs */
  typedef TMaskImage                           MaskImageType;
  typedef typename MaskImageType::Pointer      MaskImagePointer;
  typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
  typedef typename MaskImageType::PixelType    MaskPixelType;

   /** Type of the output list sample */
  typedef itk::Statistics::ListSample< MeasurementVectorType >  ListSampleType;

  /** Superclass typedefs for Measurement vector, measurement,
   * Instance Identifier, frequency, size, size element value */
  typedef PixelTraits< typename ImageType::PixelType > PixelTraitsType;
  typedef typename ListSampleType::MeasurementVectorSizeType
                                     MeasurementVectorSizeType;

  typedef DataObject::Pointer DataObjectPointer;

  /** ListSample is not a DataObject, we need to decorate it to push it down
   * a ProcessObject's pipeline */
  typedef DataObjectDecorator< ListSampleType >  ListSampleOutputType;

  /** the number of components in a measurement vector */
  itkStaticConstMacro( MeasurementVectorSize, unsigned int,
                      PixelTraitsType::Dimension );

  /** Standard itk::ProcessObject subclass method. */
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput( DataObjectPointerArraySizeType idx )
    override;

  virtual void SetMeasurementVectorSize( const MeasurementVectorSizeType s )
    {
    // Measurement vector size for this class is fixed as the pixel's
    // dimension. This method should throw an exception if the user tries to
    // set the dimension to a different value.
    if( s != MeasurementVectorSize )
      {
      itkExceptionMacro(
        << "Measurement vector size for the image adaptor obtained"
        << " from the pixel dimension is: "
        << MeasurementVectorSize << " but you "
        << "are setting it to " << s );
      }
    }

  unsigned int GetMeasurementVectorSize( void ) const
    {
    return MeasurementVectorSize;
    }

  /** Method to set/get the image */
  void SetInput( const ImageType* image );
  const ImageType* GetInput( void ) const;

  /** Method to set/get the mask */
  void SetMaskImage( const MaskImageType* image );
  const MaskImageType* GetMaskImage( void ) const;

  /** Method to get the list sample, the generated output. Note that this does
   * not invoke Update(). You should have called update on this class to get
   * any meaningful output. */
  const ListSampleType * GetListSample( void ) const;

  /** Set the pixel value treated as on in the mask. If a mask has been
   * specified, only pixels with this value will be added to the list sample, if
   * no mask has been specified all pixels will be added as measurement vectors
   * to the list sample. */
  void SetMaskValue( const MaskPixelType maskValue );
  itkGetMacro( MaskValue, MaskPixelType );

  itkSetMacro( UseSingleMaskValue, bool );
  itkGetMacro( UseSingleMaskValue, bool );

  /** This method causes the filter to generate its output. */
  virtual void GenerateData( void ) override;

  /** This method ensures that a mask image if specified has requested regions
   * that at least contain the input image's buffered region. */
  virtual void GenerateInputRequestedRegion( void )
    throw( InvalidRequestedRegionError ) override;

  virtual void GenerateOutputInformation( void ) override;

protected:
  VectorImageToListGenerator( void );
  virtual ~VectorImageToListGenerator( void ) {}

  void PrintSelf( std::ostream& os, Indent indent ) const override;

private:
  VectorImageToListGenerator( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * )
    override {};

  MaskPixelType       m_MaskValue;
  bool                m_UseSingleMaskValue;

}; // End class VectorImageToListGenerator

} // End namespace Statistics

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeVectorImageToListGenerator.hxx"
#endif

#endif // End !defined( VectorImageToListGenerator )
