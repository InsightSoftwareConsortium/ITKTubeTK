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

#ifndef __itktubeRobustMeanAndSigmaImageBuilder_h
#define __itktubeRobustMeanAndSigmaImageBuilder_h

#include "itktubeMeanAndSigmaImageBuilder.h"
#include "tubeMessage.h"

namespace itk
{

namespace tube
{


/** \class RobustMeanAndSigmaImageBuilder
 * \brief Class builds the median and robust variance from inputed images.
 * Images are processed and n/2 are kept to maintain the median calculation.
 * Therefore the number of images must be inputed prior to the start of the
 * image addition, as must be the number of outlier images to crop from the
 * ends.
 *
 * All Inputed images are assumed to have the same spacing, origin.
 * Optionally, if the tag DynamicallyAdjustOutputSize() is used, then the
 * size many vary and the images will be updated to insure the largest region
 * of all images are collected.
 *
 * This class derives from \sa MeanAndSigmaImageBuilder
 */
template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
class RobustMeanAndSigmaImageBuilder
: public MeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
{
public:

  typedef RobustMeanAndSigmaImageBuilder                    Self;
  typedef MeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
    TOutputSigmaImageType>                                  Superclass;

  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  itkNewMacro( Self );
  itkOverrideGetNameOfClassMacro( RobustMeanAndSigmaImageBuilder);

  typedef TInputImageType                              InputImageType;
  typedef TOutputMeanImageType                         OutputMeanImageType;
  typedef TOutputSigmaImageType                        OutputSigmaImageType;

  typedef typename Superclass::InputPixelType          InputPixelType;
  typedef typename Superclass::OutputMeanPixelType     OutputMeanPixelType;
  typedef typename Superclass::OutputSigmaPixelType    OutputSigmaPixelType;

  typedef typename Superclass::InputImagePointer       InputImagePointer;
  typedef typename Superclass::OutputMeanImagePointer  OutputMeanImagePointer;
  typedef typename Superclass::OutputSigmaImagePointer OutputSigmaImagePointer;

  typedef typename Superclass::RegionType              RegionType;
  typedef typename Superclass::SpacingType             SpacingType;
  typedef typename Superclass::PointType               PointType;
  typedef typename Superclass::SizeType                SizeType;

  /**
   * Add an image to the group being summed. No check is made to insure
   * the same image is not added twice
   *
   * It is assumed that the input image has the same spacing, origin & size
   * ( unless DynamicallyAdjustOutputSize() is set )
   */
  void AddImage( InputImagePointer );

  /**
   * Must call this function to finish image additions and form the mean
   * and variance. Any images added ( AddImage() ) after the call of this
   * function will eliminate the results of this output and will start a new
   * summation.
   */
  void FinalizeOutput( void );

  /**
   * Get the number of outliers on either side of the mean to remove prior
   * to statistical calculations
   */
  itkGetConstMacro( NumberOfOutlierImagesToRemove, unsigned int );

  /**
   * Set the number of outliers on either side of the mean to remove prior
   * to statistical calculations
   */
  itkSetMacro( NumberOfOutlierImagesToRemove, unsigned int );

  /**
   * Function to define class to find the median image. Requires that the
   * total number of images to be added to the form mean. Median will be
   * returned using the GetOutputMeanImage() function.
   */
  void UseMedianImage( unsigned int totalNumberOfImages )
    { m_TotalNumberOfImages = totalNumberOfImages; }

  itkGetConstMacro( TotalNumberOfImages, unsigned int );

  /**
   * Update the output images to the inputed size. Can be called at any
   * point during the mean building process.
   *
   * This does NOT test or warn if the current output image size is
   * decreased in any axis
   */
  void UpdateOutputImageSize( SizeType );


protected:

   RobustMeanAndSigmaImageBuilder( void );
  ~RobustMeanAndSigmaImageBuilder( void ) {}

  /** Processing image types */
  typedef typename Superclass::ProcessImageType         ProcessImageType;
  typedef typename Superclass::CountImageType           CountImageType;

  typedef typename Superclass::ProcessPixelType         ProcessPixelType;
  typedef typename Superclass::CountPixelType           CountPixelType;

  typedef typename Superclass::ProcessImagePointer      ProcessImagePointer;
  typedef typename Superclass::CountImagePointer        CountImagePointer;

  typedef ImageRegionIterator<InputImageType>           InputIteratorType;

  typedef typename Superclass::InputConstIteratorType
    InputConstIteratorType;
  typedef typename Superclass::ProcessConstIteratorType
    ProcessConstIteratorType;
  typedef typename Superclass::ProcessIteratorType
    ProcessIteratorType;
  typedef typename Superclass::CountConstIteratorType
    CountConstIteratorType;
  typedef typename Superclass::CountIteratorType
    CountIteratorType;
  typedef typename Superclass::OutputMeanIteratorType
    OutputMeanIteratorType;
  typedef typename Superclass::OutputSigmaIteratorType
    OutputSigmaIteratorType;

  typedef typename std::vector<InputImagePointer>       InputImageListType;

  /**
   * Build new processing images ( i.e., sumImage, sumSquareImage,
   * validCountImage )
   */
  void BuildProcessingImages( InputImagePointer i );

  /**
   * Get the ordered image list representing the lower half of voxel values.
   * Used to determine median.  Order is ascending pixel value
   */
  InputImageListType& GetLowerOutlierImages( void )
    { return m_LowerOutlierImages; }

  /**
   * Set the ordered image list representing the lower half of voxel values.
   * Used to determine median.
   */
  void  SetLowerOutlierImages( InputImageListType& list )
    { m_LowerOutlierImages = list; }

  /**
   * Get the ordered image list representing the upper outlier values.
   * Size is dependent on the number of outliers to remove for robust
   * standard deviation calculation.  Used to determine sigma calculations.
   * Order is descending starting from highest intensity
   */
  InputImageListType& GetUpperOutlierImages( void )
    { return m_UpperOutlierImages; }

  /**
   * Set the ordered image list representing the upper outlier values.
   * Size is dependent on the number of outliers to remove for robust
   * standard deviation calculation. Used to determine sigma calculations.
   * Order is descending starting from highest intensity
   */
  void SetUpperOutlierImages( InputImageListType& list )
    { m_UpperOutlierImages = list; }

  void SetUpperImages( InputImagePointer i );
  void SetLowerImages( InputImagePointer i );
  void SetMedianImages( InputImagePointer i );

  InputImagePointer GetImageCopy( InputImagePointer );
  InputImageListType& AddToUpdateImageList( InputImagePointer input,
    InputImageListType& list, bool ListIsAscending );

  bool UseMedian( void )
    { return ( m_TotalNumberOfImages > 0 ); }

  /**
   * Builds median image from the LowerOutlierImages. Can only
   * be ( reasonably ) called after at least 1 image has been added
   */
  OutputMeanImagePointer  GetMedianImage();

private:

  InputImageListType                      m_LowerOutlierImages;
  InputImageListType                      m_UpperOutlierImages;

  unsigned int                            m_NumberOfOutlierImagesToRemove;
  unsigned int                            m_TotalNumberOfImages;

}; // End class RobustMeanAndSigmaImageBuilder

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRobustMeanAndSigmaImageBuilder.hxx"
#endif

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeRobustMeanAndSigmaImageBuilder_h )
