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

#ifndef __itktubeJointHistogramImageFunction_h
#define __itktubeJointHistogramImageFunction_h

#include <itkImage.h>
#include <itkImageFunction.h>

namespace itk
{

namespace tube
{

/** \class JointHistogramImageFunction
 *
 *  Using an input image an input mask, this function allows the computation of
 *  Z-Score values at a given point in the image. This is done by calling the
 *  Precompute function to build up a mean and standard deviation for each bin
 *  in a joint-histogram. That mean and standard deviation is used to compute
 *  the Z-Score at a point when Evaluate is called. The neighborhood used in
 *  the computation of the joint histogram is determined by the feature width.
 */
template< class TInputImage, class TCoordRep = float >
class JointHistogramImageFunction
  : public ImageFunction< TInputImage, double, TCoordRep >
{
public:

  /** Class typedefs */
  typedef JointHistogramImageFunction                   Self;
  typedef ImageFunction<TInputImage, double, TCoordRep> Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;
  typedef typename Superclass::InputImageType           InputImageType;
  typedef typename TInputImage::PixelType               PixelType;
  typedef typename Superclass::PointType                PointType;
  typedef typename Superclass::IndexType                IndexType;
  typedef typename Superclass::ContinuousIndexType      ContinuousIndexType;
  typedef itk::Image<float, 2>                          HistogramType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( JointHistogramImageFunction, ImageFunction );

  /** Standard New Macro. */
  itkNewMacro( Self );

  /** Constant for fetching the dimensions of the image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

  /** Get/Set the width of a significant feature. */
  itkGetMacro( FeatureWidth, double );
  itkSetMacro( FeatureWidth, double );

  /**
   * Get/Set the factor used to compensate for standard deviations of 0 in
   * the Z-Score calculation. */
  itkGetMacro( StdevBase, double );
  itkSetMacro( StdevBase, double );

  /** Override the Set for the InputImage */
  virtual void SetInputImage( const InputImageType * ptr );

  /** Set the mask or second image used in the comparison. */
  virtual void SetInputMask( const typename InputImageType::Pointer mask );

  /** Get the mask or second image used in the comparison. */
  virtual typename InputImageType::Pointer GetInputMask( void ) const
    {
    return m_InputMask;
    }

  /** Get the size of the histogram ( It will be a size x size image ). */
  itkGetMacro( HistogramSize, unsigned int );

  /**
   * Set the size of the histogram. This will reset the bins in the
   * mean and standard deviation to zero.
   */
  virtual void SetHistogramSize( const unsigned int & size );

  /** Get the Z-score at a given point. */
  virtual double Evaluate( const PointType & point ) const
    {
    IndexType index;
    this->ConvertPointToNearestIndex( point, index );
    return ( this->EvaluateAtIndex( index ) );
    }

  /** Get the Z-score at a given continuous index. */
  virtual double EvaluateAtContinuousIndex(
    const ContinuousIndexType & index ) const
    {
    IndexType nindex;

    this->ConvertContinuousIndexToNearestIndex( index, nindex );
    return this->EvaluateAtIndex( nindex );
    }

  /** Get the Z-score at a given index. */
  virtual double EvaluateAtIndex( const IndexType & index ) const;

  /**
   * Add histograms ( based on a given point ) to the internals used to
   * calculate the mean and standard deviation histograms when needed.
   */
  virtual void Precompute( const PointType & point )
    {
    IndexType index;
    this->ConvertPointToNearestIndex( point, index );
    this->PrecomputeAtIndex( index );
    }

  /**
   * Add histograms ( based on a given continuous index ) to the internals
   * used to calculate the mean and standard deviation histograms when
   * needed.
   */
  virtual void PrecomputeAtContinuousIndex(
    const ContinuousIndexType & index )
    {
    IndexType nindex;

    this->ConvertContinuousIndexToNearestIndex( index, nindex );
    this->PrecomputeAtIndex( nindex );
    }

  itkGetObjectMacro( MeanHistogram, HistogramType );
  itkGetObjectMacro( StandardDeviationHistogram, HistogramType );

  itkSetObjectMacro( MeanHistogram, HistogramType );
  itkSetObjectMacro( StandardDeviationHistogram, HistogramType );

  itkSetMacro( ForceDiagonalHistogram, bool );
  itkGetMacro( ForceDiagonalHistogram, bool );

  // setmeanhistogram
  // setstandardeviationhistogram

  /**
   * Add histograms ( based on a given index ) to the internals used to
   * calculate the mean and standard deviation histograms when needed.
   */
  virtual void PrecomputeAtIndex( const IndexType & index );

  /**
   * Compute the mean and standard deviation histograms for use in Z-score
   * calculation.
   */
  void ComputeMeanAndStandardDeviation( void ) const;

protected:

  /** Default constructor */
  JointHistogramImageFunction( void );

  /** Default destructor */
  ~JointHistogramImageFunction( void ) {}

  /** PrintSelf function for introspection. */
  void PrintSelf( std::ostream & os, Indent indent ) const;

  typename HistogramType::Pointer & ComputeHistogramAtIndex(
    const IndexType & index, bool blur=true ) const;

  /** Get the Z-score at a given index. */
  double ComputeZScoreAtIndex( const IndexType & index ) const;

  /** Data members **/
  typename InputImageType::Pointer         m_InputMask;
  mutable typename HistogramType::Pointer  m_Histogram;
  mutable typename HistogramType::Pointer  m_SumHistogram;
  mutable typename HistogramType::Pointer  m_SumOfSquaresHistogram;
  mutable typename HistogramType::Pointer  m_MeanHistogram;
  mutable typename HistogramType::Pointer  m_StandardDeviationHistogram;
  double                                   m_FeatureWidth;
  double                                   m_StdevBase;
  unsigned int                             m_HistogramSize;
  mutable unsigned int                     m_NumberOfSamples;
  mutable unsigned int                     m_NumberOfComputedSamples;
  double                                   m_ImageMin;
  double                                   m_ImageMax;
  double                                   m_ImageStep;
  double                                   m_MaskMin;
  double                                   m_MaskMax;
  double                                   m_MaskStep;

  bool                                     m_ForceDiagonalHistogram;

private:
  JointHistogramImageFunction( const Self & ); // Purposely not implemented
  void operator=( const Self & ); // Purposely not implemented

}; // End class JointHistogramImageFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeJointHistogramImageFunction.hxx"
#endif

#endif // End !defined( __itktubeJointHistogramImageFunction_h )
