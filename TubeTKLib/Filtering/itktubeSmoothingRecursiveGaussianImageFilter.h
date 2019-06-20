/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 ( Tue, 10 Jul 2007 ) $
  Version:   $Revision: 0 $

  Copyright ( c ) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itktubeSmoothingRecursiveGaussianImageFilter_h
#define __itktubeSmoothingRecursiveGaussianImageFilter_h

#include "itkRecursiveGaussianImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkPixelTraits.h"
#include "itkCommand.h"

namespace itk
{

namespace tube
{

/**
 * This is a copy of the itkSmoothingRecursiveGaussianImageFilter.
 *
 * It fixes one line to correct a zero-size array error due to
 *   order of variable declaration.
 *
 *
 */
template< typename TInputImage,
          typename TOutputImage = TInputImage >
class SmoothingRecursiveGaussianImageFilter:
  public InPlaceImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef SmoothingRecursiveGaussianImageFilter           Self;
  typedef InPlaceImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Pixel Type of the input image */
  typedef TInputImage                                     InputImageType;
  typedef TOutputImage                                    OutputImageType;
  typedef typename TInputImage::PixelType                 PixelType;
  typedef typename NumericTraits< PixelType >::RealType   RealType;

  typedef typename NumericTraits< PixelType >::ScalarRealType
    ScalarRealType;

  /** Runtime information support. */
  itkTypeMacro( SmoothingRecursiveGaussianImageFilter,
    ImageToImageFilter );

  /** Image dimension. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Define the type for the sigma array */
  typedef FixedArray< ScalarRealType,
    itkGetStaticConstMacro( ImageDimension ) > SigmaArrayType;

  /** Define the image type for internal computations
    RealType is usually 'double' in NumericTraits.
    Here we prefer float in order to save memory.  */
  typedef typename NumericTraits< PixelType >::FloatType
    InternalRealType;

  typedef typename InputImageType::template Rebind<InternalRealType>::Type
    RealImageType;

  /**  The first in the pipeline  */
  typedef RecursiveGaussianImageFilter<
    InputImageType,
    RealImageType
    >    FirstGaussianFilterType;

  /**  Smoothing filter type */
  typedef RecursiveGaussianImageFilter<
    RealImageType,
    RealImageType
    >    InternalGaussianFilterType;

  /**  The last in the pipeline  */
  typedef CastImageFilter<
    RealImageType,
    OutputImageType
    >    CastingFilterType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Set Sigma value. Sigma is measured in the units of image spacing. You
    may use the method SetSigma to set the same value across each axis or
    use the method SetSigmaArray if you need different values along each
    axis. */
  void SetSigmaArray( const SigmaArrayType & sigmas );
  void SetSigma( ScalarRealType sigma );

  SigmaArrayType GetSigmaArray() const;
  ScalarRealType GetSigma() const;

  /** This method does not effect the output of this filter.
   *
   *  \sa  RecursiveGaussianImageFilter::SetNormalizeAcrossScale
   */
  void SetNormalizeAcrossScale( bool normalizeInScaleSpace );
  itkGetConstMacro( NormalizeAcrossScale, bool );

  // See super class for doxygen documentation
  //
  void SetNumberOfWorkUnits( ThreadIdType nb ) override;

  // See super class for doxygen documentation
  //
  virtual bool CanRunInPlace( void ) const override;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  // This concept does not work with variable length vector images
  //itkConceptMacro( InputHasNumericTraitsCheck,
  //( Concept::HasNumericTraits< PixelType > ) );
  // End concept checking
#endif

protected:
  SmoothingRecursiveGaussianImageFilter();
  virtual ~SmoothingRecursiveGaussianImageFilter() {}

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Generate Data */
  void GenerateData( void ) override;

  /** SmoothingRecursiveGaussianImageFilter needs all of the input to
   * produce an output. Therefore, SmoothingRecursiveGaussianImageFilter
   * needs to provide an implementation for GenerateInputRequestedRegion
   * in order to inform * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() override;

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion( DataObject *output ) override;

private:
  SmoothingRecursiveGaussianImageFilter( const Self & ); //purposely not
                                                       // implemented
  void operator=( const Self & );                        //purposely not

  // implemented

  typedef std::vector< typename InternalGaussianFilterType::Pointer >
    SmoothingFiltersArrayType;
  SmoothingFiltersArrayType m_SmoothingFilters;

  typename FirstGaussianFilterType::Pointer    m_FirstSmoothingFilter;
  typename CastingFilterType::Pointer          m_CastingFilter;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale;

  /** Standard deviation of the gaussian used for smoothing */
  SigmaArrayType m_Sigma;
};

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSmoothingRecursiveGaussianImageFilter.hxx"
#endif

#endif
