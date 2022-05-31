/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeNJetImageFunction_h
#define __itktubeNJetImageFunction_h

#include <itkArray.h>
#include <itkImageFunction.h>
#include <itkMatrix.h>
#include <itkVector.h>

#include <vnl/vnl_c_vector.h>
#include <vnl/vnl_vector.h>

namespace itk
{

namespace tube
{

/** \class NJetImageFunction
 * \brief Calculate the Gaussian blurred value, 1st derivatives, and
 *        second derivatives at point
 *        given a scale and extent of the Gaussian.
 * This class is templated over the input image type.
 */
template< class TInputImage >
class NJetImageFunction : public Object
{
public:
  /**
   * Standard "Self" typedef
   */
  typedef NJetImageFunction Self;

  /**
   * Standard "Superclass" typedef
   */
  typedef Object Superclass; // ImageFunction<TInputImage, double> Superclass;

  /**
   * Smart pointer typedef support.
   */
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro( Self );

  itkTypeMacro( NJetImageFunction, Object );
  //itkTypeMacro( NJetImageFunction, ImageFunction );

  /**
   * Dimension of the underlying image.
   */
  itkStaticConstMacro( ImageDimension,
                      unsigned int,
                      TInputImage::ImageDimension );

  typedef TInputImage                                      InputImageType;
  typedef typename InputImageType::Pointer                 InputImagePointer;

  typedef Point<double, ImageDimension >                   PointType;
  typedef Vector< double, ImageDimension >                 VectorType;
  typedef Matrix< double, ImageDimension, ImageDimension > MatrixType;

  typedef typename InputImageType::IndexType               IndexType;

  typedef ContinuousIndex<double, ImageDimension >         ContinuousIndexType;

  typedef typename InputImageType::SpacingType             SpacingType;

  typedef Size< ImageDimension >                           SizeType;

  typedef Array< VectorType >                              ArrayVectorType;

  /**
   * Set the input image.
   */
  void SetInputImage( const InputImageType * ptr );
  itkGetConstObjectMacro( InputImage, InputImageType );

  void SetInputImageMask( const InputImageType * ptr );
  itkGetConstObjectMacro( InputImageMask, InputImageType );
  itkSetMacro( UseInputImageMask, bool );
  itkGetMacro( UseInputImageMask, bool );

  void ComputeStatistics( void );

  /** Return the min over the ( possibly masked ) image.
   * Requires previous call to ComputeStatistics */
  double GetMin( void ) const;

  /** Return the max over the ( possibly masked ) image.
   * Requires previous call to ComputeStatistics */
  double GetMax( void ) const;

  itkGetConstMacro( MostRecentIntensity, double );
  itkGetConstMacro( MostRecentDerivative, VectorType );
  itkGetConstMacro( MostRecentHessian, MatrixType );
  itkGetConstMacro( MostRecentRidgeness, double );
  itkGetConstMacro( MostRecentRidgeRoundness, double );
  itkGetConstMacro( MostRecentRidgeLevelness, double );
  itkGetConstMacro( MostRecentRidgeCurvature, double );
  itkGetConstMacro( MostRecentRidgeTangent, VectorType );

  /** Evalulate the function at specified point */
  double Evaluate( const PointType & point, double scale=1 ) const;
  double Evaluate( const PointType & point, const VectorType & v1,
                      double scale=1 ) const;
  double Evaluate( const PointType & point, const VectorType & v1,
                      const VectorType & v2,
                      double scale=1 ) const;


  /** Evaluate the function at specified Index position */
  double EvaluateAtIndex( const IndexType & index, double scale=1 ) const;
  double EvaluateAtIndex( const IndexType & index, const VectorType & v1,
                      double scale=1 ) const;
  double EvaluateAtIndex( const IndexType & index, const VectorType & v1,
                      const VectorType & v2,
                      double scale=1 ) const;

  /** Evaluate the function at specified ContinousIndex position */
  double EvaluateAtContinuousIndex( const ContinuousIndexType & index,
                      double scale=1 ) const;
  double EvaluateAtContinuousIndex( const ContinuousIndexType & index,
                      const VectorType & v1,
                      double scale=1 ) const;
  double EvaluateAtContinuousIndex( const ContinuousIndexType & index,
                      const VectorType & v1,
                      const VectorType & v2,
                      double scale=1 ) const;


  double  Derivative( const PointType & point, double scale,
                       VectorType & d ) const;
  double  Derivative( const PointType & point,
                       const VectorType & v1,
                       double scale,
                       VectorType & d ) const;
  double  Derivative( const PointType & point,
                       const VectorType & v1, const VectorType & v2,
                       double scale,
                       VectorType & d ) const;

  double  DerivativeAtIndex( const IndexType & index, double scale,
                       VectorType & d ) const;
  double  DerivativeAtIndex( const IndexType & index,
                       const VectorType & v1,
                       double scale,
                       VectorType & d ) const;
  double  DerivativeAtIndex( const IndexType & index,
                       const VectorType & v1, const VectorType & v2,
                       double scale,
                       VectorType & d ) const;

  double  DerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                       double scale,
                       VectorType & d ) const;
  double  DerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                       const VectorType & v1,
                       double scale,
                       VectorType & d ) const;
  double  DerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                       const VectorType & v1, const VectorType & v2,
                       double scale,
                       VectorType & d ) const;

  double  Hessian( const PointType & point, double scale,
                       MatrixType & m ) const;
  double  Hessian( const PointType & point,
                       const VectorType & v1, double scale,
                       MatrixType & m ) const;
  double  Hessian( const PointType & point,
                       const VectorType & v1, const VectorType & v2,
                       double scale,
                       MatrixType & m ) const;

  double  HessianAtIndex( const IndexType & index, double scale,
                       MatrixType & m ) const;
  double  HessianAtIndex( const IndexType & index,
                       const VectorType & v1, double scale,
                       MatrixType & m ) const;
  double  HessianAtIndex( const IndexType & index,
                       const VectorType & v1, const VectorType & v2,
                       double scale,
                       MatrixType & m ) const;

  double  HessianAtContinuousIndex( const ContinuousIndexType & cIndex,
                       double scale,
                       MatrixType & m ) const;
  double  HessianAtContinuousIndex( const ContinuousIndexType & cIndex,
                       const VectorType & v1, double scale,
                       MatrixType & m ) const;
  double  HessianAtContinuousIndex( const ContinuousIndexType & cIndex,
                       const VectorType & v1, const VectorType & v2,
                       double scale,
                       MatrixType & m ) const;

  double  Jet( const PointType & point, VectorType & d, MatrixType & h,
                       double scale=1 ) const;

  double  JetAtIndex( const IndexType & cIndex, VectorType & d,
                       MatrixType & h,
                       double scale=1 ) const;

  double  JetAtContinuousIndex( const ContinuousIndexType & cIndex,
                       VectorType & d, MatrixType & h,
                       double scale=1 ) const;

  double  Ridgeness( const PointType & point, double scale=1 ) const;
  double  Ridgeness( const PointType & point,
                       const VectorType & v1, double scale=1 ) const;
  double  Ridgeness( const PointType & point,
                       const VectorType & v1, const VectorType & v2,
                       double scale=1 ) const;

  double  RidgenessAtIndex( const IndexType & index, double scale=1 ) const;
  double  RidgenessAtIndex( const IndexType & index,
                       const VectorType & v1, double scale=1 ) const;
  double  RidgenessAtIndex( const IndexType & index,
                       const VectorType & v1, const VectorType & v2,
                       double scale=1 ) const;

  double  RidgenessAtContinuousIndex( const ContinuousIndexType & cIndex,
                       double scale=1 ) const;
  double  RidgenessAtContinuousIndex( const ContinuousIndexType & cIndex,
                       const VectorType & v1, double scale=1 ) const;
  double  RidgenessAtContinuousIndex( const ContinuousIndexType & cIndex,
                       const VectorType & v1, const VectorType & v2,
                       double scale=1 ) const;

  InputImagePointer ScaleSubsample( double factor );

  /**
   * Set the Extent
   */
  itkSetMacro( Extent, double );

  /**
   * Get the Extent
   */
  itkGetMacro( Extent, double );

  /**
   * Get the Spacing
   */
  itkGetConstMacro( InputImageSpacing, SpacingType );

  /**
   * If set to true, then values and derivatives are computed and then
   *   projected onto the subset of directions given.   Otherwise, the
   *   computations are done explicitly within the vector/plane given.
   */
  itkSetMacro( UseProjection, bool );
  itkGetMacro( UseProjection, bool );

protected:
  NJetImageFunction( void );

  ~NJetImageFunction( void );

  void PrintSelf( std::ostream& os, Indent indent ) const override;

  typename InputImageType::ConstPointer  m_InputImage;
  typename InputImageType::ConstPointer  m_InputImageMask;
  bool                                   m_UseInputImageMask;

  IndexType               m_InputImageMinX;
  IndexType               m_InputImageMaxX;
  SizeType                m_InputImageSize;
  SpacingType             m_InputImageSpacing;
  SpacingType             m_InputImageSpacingSquared;
  double                  m_Extent;

  mutable double          m_MostRecentIntensity;
  mutable VectorType      m_MostRecentDerivative;
  mutable MatrixType      m_MostRecentHessian;
  mutable double          m_MostRecentRidgeness;
  mutable double          m_MostRecentRidgeRoundness;
  mutable double          m_MostRecentRidgeLevelness;
  mutable double          m_MostRecentRidgeCurvature;
  mutable VectorType      m_MostRecentRidgeTangent;

  double                  m_CurvatureExpectedMax;

  bool                    m_ValidStats;
  double                  m_StatsMin;
  double                  m_StatsMax;

  bool                    m_UseProjection;

private:
  NJetImageFunction( const Self& );
  void operator=( const Self& );

}; // End class NJetImageFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeNJetImageFunction.hxx"
#endif

#endif // End !defined( __itktubeNJetImageFunction_h )
