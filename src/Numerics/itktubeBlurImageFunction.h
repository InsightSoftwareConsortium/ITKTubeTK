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

#ifndef __itktubeBlurImageFunction_h
#define __itktubeBlurImageFunction_h

#include <itkImageFunction.h>
#include <itkIndex.h>

namespace itk
{

namespace tube
{

/** \class BlurImageFunction
 * \brief Calculate the Gaussian blurred value at point
 *        given a scale and extent of the Gaussian.
 * This class is templated over the input image type.
 */
template< class TInputImage >
class BlurImageFunction
  : public ImageFunction< TInputImage, double, double >
{
public:
  /**
   * Standard "Self" typedef */
  typedef BlurImageFunction                            Self;
  typedef ImageFunction<TInputImage, double, double>   Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  itkOverrideGetNameOfClassMacro( BlurImageFunction);

  itkNewMacro( Self );

  /**
   * InputImageType typedef support. */
  typedef TInputImage                                InputImageType;
  typedef typename InputImageType::SpacingType       SpacingType;
  typedef typename InputImageType::SpacingType       SizeType;

  /**
   * IndexType typedef support. */
  typedef typename InputImageType::IndexType         IndexType;
  typedef typename Superclass::ContinuousIndexType   ContinuousIndexType;

  /**
   * Dimension of the underlying image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    InputImageType::ImageDimension );

  /**
   * Point typedef support. */
  typedef typename Superclass::PointType             PointType;

  /**
   * Set the input image. */
  virtual void SetInputImage( const InputImageType * ptr ) override;

  /**
   * Evalulate the function at specified point */
  virtual double Evaluate( const PointType& point ) const override;

  /** Evaluate the function at specified Index position. */
  virtual double EvaluateAtIndex( const IndexType & index ) const override;

  /** Evaluate the function at specified ContinousIndex position. */
  virtual double EvaluateAtContinuousIndex( const ContinuousIndexType &
    index ) const override;

  /**
   * Set the Scale */
  void SetScale( double scale );

  /**
   * Get the Scale */
  itkGetMacro( Scale, double );

  /**
   * Set the Extent */
  void SetExtent( double extent );

  /**
   * Get the Extent */
  itkGetMacro( Extent, double );

  /**
   * Get the Spacing */
  itkGetMacro( Spacing, SpacingType );

  /**
   * Interpret the sigma value to be in terms of x-spacing */
  void SetUseRelativeSpacing( bool useRelativeSpacing );

  /**
   * Get the Spacing */
  itkGetMacro( UseRelativeSpacing, bool );

protected:

  BlurImageFunction( void );
  virtual ~BlurImageFunction( void ) {}

  void PrintSelf( std::ostream& os, Indent indent ) const override;

  void RecomputeKernel( void );

private:

  BlurImageFunction( const Self& );
  void operator=( const Self& );

  typedef std::list< double > KernelWeightsListType;

  typedef std::list< typename InputImageType::IndexType > KernelXListType;

  bool                    m_UseRelativeSpacing;
  SpacingType             m_Spacing;
  SpacingType             m_OriginalSpacing;
  double                  m_Scale;
  double                  m_Extent;
  KernelWeightsListType   m_KernelWeights;
  KernelXListType         m_KernelX;
  IndexType               m_KernelMin;
  IndexType               m_KernelMax;
  SizeType                m_KernelSize;
  double                  m_KernelTotal;

  IndexType               m_ImageIndexMin;
  IndexType               m_ImageIndexMax;

}; // End class BlurImageFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeBlurImageFunction.hxx"
#endif

#endif // End !defined( __itktubeBlurImageFunction_h )
