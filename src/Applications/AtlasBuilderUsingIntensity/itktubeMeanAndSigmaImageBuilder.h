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

#ifndef __itktubeMeanAndSigmaImageBuilder_h
#define __itktubeMeanAndSigmaImageBuilder_h

#include "tubeMessage.h"

#include <itkObject.h>
#include <itkResampleImageFilter.h>

namespace itk
{

namespace tube
{

/** \class MeanAndSigmaImageBuilder
 * \brief Class builds the mean and variance from inputed images.
 *   Images are processed and discarded as they are entered,
 *   which prevents memory overload.
 *
 *  All inputed images entered to AddImage() are assumed to have the
 *  same spacing, origin.  Optionally, the output size can be changed by using
 *  AdjustOutputImageSize() and all subsequent additions to AddImage()
 *  must have that same size.  This allows for dynamic changes to the size
 *  of the output.
 *
 *  NOTE: This is not done for origin or spacing, because those factors
 *  could require interpolation, which would change the maintain base results.
 */
template< class TInputImageType, class TOutputMeanImageType,
          class TOutputSigmaImageType >
class MeanAndSigmaImageBuilder : public Object
{
public:

  typedef MeanAndSigmaImageBuilder                      Self;
  typedef Object                                        Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImageType::ImageDimension );

  itkNewMacro( Self );
  itkTypeMacro( MeanAndSigmaImageBuilder, Object );

  typedef TInputImageType                               InputImageType;
  typedef TOutputMeanImageType                          OutputMeanImageType;
  typedef TOutputSigmaImageType                         OutputSigmaImageType;
  typedef Image<
    float, itkGetStaticConstMacro( ImageDimension )>           CountImageType;

  typedef typename InputImageType::PixelType            InputPixelType;
  typedef typename OutputMeanImageType::PixelType       OutputMeanPixelType;
  typedef typename OutputSigmaImageType::PixelType      OutputSigmaPixelType;
  typedef typename CountImageType::PixelType            CountPixelType;

  typedef typename InputImageType::Pointer              InputImagePointer;
  typedef typename OutputMeanImageType::Pointer         OutputMeanImagePointer;
  typedef typename OutputSigmaImageType::Pointer        OutputSigmaImagePointer;
  typedef typename CountImageType::Pointer              CountImagePointer;

  typedef typename InputImageType::RegionType           RegionType;
  typedef typename InputImageType::SizeType             SizeType;
  typedef typename InputImageType::SpacingType          SpacingType;
  typedef typename InputImageType::PointType            PointType;

  typedef Image< float, itkGetStaticConstMacro( ImageDimension ) >
    ProcessImageType;

  /**
   * Add an image to the group being summed.
   *
   * NOTE: No check is made to insure the same image is not added twice
   *
   *  It is assumed that the input image has the same spacing,
   *  origin & size ( unless DynamicallyAdjustOutputSize() is set )
   */
  virtual void AddImage( InputImagePointer );

  /**
   * Must call this function to finish image additions and form the
   * mean and variance
   *
   * NOTE: Any images added ( AddImage() ) after the call of this
   * function will eliminate the results of this output
   * and will start a new summation.
   */
  virtual void FinalizeOutput( void );

  /** Get the output mean image */
  itkGetModifiableObjectMacro( OutputMeanImage, OutputMeanImageType );

  /** Get the output variance image */
  itkGetModifiableObjectMacro( OutputSigmaImage, OutputSigmaImageType );

  /** Get the output valid count image */
  itkGetModifiableObjectMacro( ValidCountImage, CountImageType );

  /**
   * Get the minimum number of contributing images for a given voxel
   * to count in the output images
   */
  itkGetConstMacro( ImageCountThreshold, unsigned short );

  /**
   * Set the minimum number of contributing images for a given voxel
   * to count in the output images
   */
  itkSetMacro( ImageCountThreshold, unsigned short );

  /**
   * Get the bool flag of whether the OutputSigmaImageType is
   * variance ( false ) or standard deviation ( true )
   */
  itkGetConstMacro( UseStandardDeviation, bool );

  /**
   * Set the bool flag of whether the OutputSigmaImageType is
   * variance ( false ) or standard deviation ( true ).  Default is true
   */
  itkSetMacro( UseStandardDeviation, bool );

  /**
   * Get the threshold lower value where all values less than or
   * equal to input are not counted as a valid voxel
   */
  itkGetConstMacro( ThresholdInputImageBelow, InputPixelType );

  /**
   * Set the threshold lower value where all values less than or
   * equal to input are not counted as a valid voxel - default is 0
   */
  void SetThresholdInputImageBelow( InputPixelType value )
    {
    m_ThresholdInputImageBelow = value;
    this->SetThresholdInputImageBelowOn( true );
    }

  /**
   * Turn on and off the ThresholdInputImageBelow function ( default is off ).
   * Will use last inputed value, or default otherwise
   * ( see ThresholdInputImageBelow )
   */
  itkGetConstMacro( ThresholdInputImageBelowOn, bool );

  /**
   * Turn on and off the ThresholdInputImageBelow function ( default is off ).
   * Will use last inputed value, or default otherwise
   * ( see ThresholdInputImageBelow )
   */
  itkSetMacro( ThresholdInputImageBelowOn, bool );

  /**
   * Get the bool flag of whether the output images should adjust their
   * size dynamically based on the inputed image sizes.  Default is false
   */
  itkGetConstMacro( DynamicallyAdjustOutputSize, bool );

  /**
   * Set the bool flag of whether the output images should adjust their
   * size dynamically based on the inputed image sizes.  Default is false
   */
  itkSetMacro( DynamicallyAdjustOutputSize, bool );

  /**
   * Update the output images to the inputed size.
   *
   * Can be called at any point in the mean building process AFTER
   * the first image has been added.
   *
   * NOTE: this does NOT test or warn if the current output image
   * size is decreased in any axis
   */
  virtual void UpdateOutputImageSize( SizeType );

  /** Get the current size of the output images */
  itkGetConstReferenceMacro( OutputSize, SizeType );

  /** Set the current size of the output images */
  itkSetMacro( OutputSize, SizeType );

  /** Get the current size of the output images */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  /** Set the current size of the output images */
  itkSetMacro( OutputSpacing, SpacingType );

  /** Get the current size of the output images */
  itkGetConstReferenceMacro( OutputOrigin, PointType );

  /** Set the current size of the output images */
  itkSetMacro( OutputOrigin, PointType );

  itkGetModifiableObjectMacro( SumImage, ProcessImageType );

protected:

  MeanAndSigmaImageBuilder( void );
  ~MeanAndSigmaImageBuilder( void ) {}

  /** Processing image types */
  typedef typename ProcessImageType::PixelType         ProcessPixelType;
  typedef typename ProcessImageType::Pointer           ProcessImagePointer;

  typedef ImageRegionConstIterator< InputImageType >   InputConstIteratorType;
  typedef ImageRegionConstIterator< ProcessImageType > ProcessConstIteratorType;
  typedef ImageRegionIterator< ProcessImageType >      ProcessIteratorType;
  typedef ImageRegionConstIterator< CountImageType >   CountConstIteratorType;
  typedef ImageRegionIterator< CountImageType >        CountIteratorType;
  typedef ImageRegionIterator< OutputMeanImageType >   OutputMeanIteratorType;
  typedef ImageRegionIterator< OutputSigmaImageType >  OutputSigmaIteratorType;

  itkSetObjectMacro( SumImage, ProcessImageType );

  itkGetModifiableObjectMacro( SumSquareImage, ProcessImageType );
  itkSetObjectMacro( SumSquareImage, ProcessImageType );

  itkSetObjectMacro( ValidCountImage, CountImageType );

  /** Set the output mean image */
  itkSetObjectMacro( OutputMeanImage, OutputMeanImageType );

  /** Set the output variance image */
  itkSetObjectMacro( OutputSigmaImage, OutputSigmaImageType );

  /**
   * Get processing variable flag.  Flag states whether the Processing
   * Images have been defined
   */
  itkGetConstMacro( IsProcessing, bool );

  /**
   * Set processing variable flag.  Flag states whether the Processing
   * Images have been defined
   */
  itkSetMacro( IsProcessing, bool );

  /**
   * Build new processing images, i.e.,
   * sumImage, sumSquareImage, validCountImage
   */
  virtual void BuildProcessingImages( InputImagePointer i );

private:

  ProcessImagePointer                     m_SumImage;
  ProcessImagePointer                     m_SumSquareImage;
  CountImagePointer                       m_ValidCountImage;

  OutputMeanImagePointer                  m_OutputMeanImage;
  OutputSigmaImagePointer                 m_OutputSigmaImage;

  unsigned short                          m_ImageCountThreshold;
  bool                                    m_ThresholdInputImageBelowOn;
  InputPixelType                          m_ThresholdInputImageBelow;
  bool                                    m_IsProcessing;
  bool                                    m_UseStandardDeviation;
  bool                                    m_DynamicallyAdjustOutputSize;

  SizeType                                m_OutputSize;
  SpacingType                             m_OutputSpacing;
  PointType                               m_OutputOrigin;

}; // End class MeanAndSigmaImageBuilder

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeMeanAndSigmaImageBuilder.hxx"
#endif

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMeanAndSigmaImageBuilder_h )
