/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __itktubeCropImageFilter_h
#define __itktubeCropImageFilter_h

#include "itkImage.h"
#include "itkExtractImageFilter.h"

// Forward declare itkTubeTK class to allow friendship
namespace tube
{
template< typename tube_TInputImage, typename tube_TOutputImage >
class CropImage;
}

namespace itk
{
namespace tube
{
/** \class CropImageFilter
 * \brief Decrease the image size by cropping the image by an itk::Size at
 * both the upper and lower bounds of the largest possible region.
 * CropImageFilter changes the image boundary of an image by removing
 * pixels outside the target region.  The target region is not specified in
 * advance, but calculated in BeforeThreadedGenerateData().
 *
 * This filter uses ExtractImageFilter to perform the cropping.
 *
 */

template< typename TInputImage, typename TOutputImage >
class CropImageFilter:
  public ExtractImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef CropImageFilter                                 Self;
  typedef ExtractImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( CropImageFilter, ExtractImageFilter );

  /** Typedef to describe the output and input image region types. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::InputImageRegionType  InputImageRegionType;

  /** Typedef to describe the type of pixel. */
  typedef typename Superclass::OutputImagePixelType OutputImagePixelType;
  typedef typename Superclass::InputImagePixelType  InputImagePixelType;

  /** Typedef to describe the output and input image index and size types. */
  typedef typename Superclass::OutputImageIndexType OutputImageIndexType;
  typedef typename Superclass::InputImageIndexType  InputImageIndexType;
  typedef typename Superclass::OutputImageSizeType  OutputImageSizeType;
  typedef typename Superclass::InputImageSizeType   InputImageSizeType;
  typedef InputImageSizeType                        SizeType;
  typedef TInputImage                               ImageType;

  /** ImageDimension constants */
  itkStaticConstMacro( InputImageDimension, unsigned int,
                      Superclass::InputImageDimension );
  itkStaticConstMacro( OutputImageDimension, unsigned int,
                      Superclass::OutputImageDimension );

  /** Set/Get the cropping sizes for the upper and lower boundaries. */
  // itkCropImageFilter parameters
  itkSetMacro( UpperBoundaryCropSize, SizeType );
  itkGetConstMacro( UpperBoundaryCropSize, SizeType );

  itkSetMacro( LowerBoundaryCropSize, SizeType );
  itkGetConstMacro( LowerBoundaryCropSize, SizeType );

  // tubeCropROI parameters
  void SetMin( typename ImageType::IndexType roiMin );
  typename ImageType::IndexType GetMin( void );
  itkSetMacro( UseROIMin, bool );

  void SetMax( typename ImageType::IndexType roiMax );
  typename ImageType::IndexType GetMax( void );
  itkSetMacro( UseROIMax, bool );

  void SetSize( typename ImageType::SizeType roiSize );
  typename ImageType::SizeType GetSize( void );
  itkSetMacro( UseROISize, bool );

  void SetCenter( typename ImageType::IndexType roiCenter );
  typename ImageType::IndexType GetCenter( void );
  itkSetMacro( UseROICenter, bool );

  void SetBoundary( typename ImageType::IndexType roiBoundary );
  typename ImageType::IndexType GetBoundary( void );
  itkSetMacro( UseROIBoundary, bool );

  void SetMatchVolume( const ImageType * matchVolume );

  void SetMatchMask( const ImageType * maskImage );

  void SetSplitInput( typename ImageType::IndexType splitIndex,
    typename ImageType::IndexType roiIndex );

  void SetBoundaryCropSize( const SizeType & s )
    {
    this->SetUpperBoundaryCropSize( s );
    this->SetLowerBoundaryCropSize( s );
    }

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputConvertibleToOutputCheck,
    ( Concept::Convertible< InputImagePixelType, OutputImagePixelType > ) );
  itkConceptMacro( SameDimensionCheck,
    ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
  // End concept checking
#endif

protected:
  CropImageFilter( void );
  ~CropImageFilter() {}

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void GenerateOutputInformation() override;

private:
  /** itkCropImageFilter parameters */
  CropImageFilter( const Self & );
  void operator=( const Self & );

  SizeType m_UpperBoundaryCropSize;
  SizeType m_LowerBoundaryCropSize;

  /** tubeCropROI parameters */
  typename ImageType::IndexType  m_ROIMin;
  bool                           m_UseROIMin;
  typename ImageType::IndexType  m_ROIMax;
  bool                           m_UseROIMax;
  typename ImageType::SizeType   m_ROISize;
  bool                           m_UseROISize;
  typename ImageType::IndexType  m_ROICenter;
  bool                           m_UseROICenter;
  typename ImageType::IndexType  m_ROIBoundary;
  bool                           m_UseROIBoundary;
  float                          m_ProgressStart;
  float                          m_ProgressRange;

  /** friendship facilitating itkTukeTK integration */
  template< typename tube_TInputImage, typename tube_TOutputImage >
  friend class ::tube::CropImage;
};
} // End namespace tube
} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeCropImageFilter.hxx"
#endif

#endif // End !defined( __itktubeCropImageFilter_h )
