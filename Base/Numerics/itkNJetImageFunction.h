/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: main.cxx,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 (Tue, 10 Jul 2007) $
  Version:   $Revision: 48 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNJetImageFunction_h
#define __itkNJetImageFunction_h

#include "itkImageFunction.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkArray.h"


// Remove vnl warning
#ifdef VNL_CONFIG_CHECK_BOUNDS
#undef VNL_CONFIG_CHECK_BOUNDS
#endif

#include "vnl/vnl_vector.h"
#include "vnl/vnl_c_vector.txx"

namespace itk
{

/**
 * \class NJetImageFunction
 * \brief Calculate the gaussian blurred value, 1st derivatives, and
 *        second derivatives at point 
 *        given a scale and extent of the gaussian.
 * This class is templated over the input image type.
 *
 */
template <class TInputImage>
class ITK_EXPORT NJetImageFunction :
public Object // ImageFunction< TInputImage, double >
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
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  itkTypeMacro( NJetImageFunction, Object );
  //itkTypeMacro( NJetImageFunction, ImageFunction );

  /**
   * Dimension of the underlying image.
   */
  itkStaticConstMacro(ImageDimension,
                      unsigned int,
                      TInputImage::ImageDimension);
  
  typedef TInputImage                                 InputImageType;
  typedef typename InputImageType::Pointer            InputImagePointer;

  typedef Point<double, itkGetStaticConstMacro(ImageDimension) >
                                                      PointType;

  typedef Vector< double, itkGetStaticConstMacro(ImageDimension) > 
                                                      VectorType;

  typedef Matrix< double, itkGetStaticConstMacro(ImageDimension),
                  itkGetStaticConstMacro(ImageDimension) >  
                                                      MatrixType;

  typedef typename InputImageType::IndexType          IndexType;

  typedef ContinuousIndex<double, itkGetStaticConstMacro(ImageDimension) > 
                                                      ContinuousIndexType;

  typedef typename InputImageType::SpacingType              SpacingType;

  typedef Size<itkGetStaticConstMacro(ImageDimension)>      SizeType; 

  typedef Array< VectorType >                               ArrayVectorType;

  /**
   * Set the input image.
   */
  void SetInputImage( const InputImageType * ptr );
  itkGetConstObjectMacro( InputImage, InputImageType );
  
  void SetInputImageMask( const InputImageType * ptr );
  itkGetConstObjectMacro( InputImageMask, InputImageType );
  itkSetMacro( UseInputImageMask, bool );
  itkGetMacro( UseInputImageMask, bool );

  void ComputeStatistics(void);

  /** Return the min over the (possibly masked) image.
   * Requires previous call to ComputeStatistics*/
  double GetMin(void) const;

  /** Return the max over the (possibly masked) image.
   * Requires previous call to ComputeStatistics*/
  double GetMax(void) const;

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

  VectorType 
    Derivative( const PointType & point, double scale=1 ) const;
  VectorType 
    Derivative( const PointType & point, 
                const VectorType & v1, 
                double scale=1 ) const;
  VectorType 
    Derivative( const PointType & point, 
                const VectorType & v1, const VectorType & v2,
                double scale=1 ) const;

  VectorType 
    DerivativeAtIndex( const IndexType & index, double scale=1 ) const;
  VectorType 
    DerivativeAtIndex( const IndexType & index, 
                       const VectorType & v1,
                       double scale=1 ) const;
  VectorType 
    DerivativeAtIndex( const IndexType & index,
                       const VectorType & v1, const VectorType & v2,
                       double scale=1 ) const;

  VectorType
    DerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                                 double scale=1 ) const;
  VectorType
    DerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                                 const VectorType & v1,
                                 double scale=1 ) const;
  VectorType
    DerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                                 const VectorType & v1, const VectorType & v2,
                                 double scale=1 ) const;

  VectorType 
    ValueAndDerivative( const PointType & point,
                        double & val,
                        double scale=1 ) const;
  VectorType 
    ValueAndDerivative( const PointType & point,
                        double & val,
                        const VectorType & v1,
                        double scale=1 ) const;
  VectorType 
    ValueAndDerivative( const PointType & point,
                        double & val,
                        const VectorType & v1, const VectorType & v2,
                        double scale=1 ) const;

  VectorType 
    ValueAndDerivativeAtIndex( const IndexType & index,
                               double & val,
                               double scale=1 ) const;
  VectorType 
    ValueAndDerivativeAtIndex( const IndexType & index,
                               double & val,
                               const VectorType & v1,
                               double scale=1 ) const;
  VectorType 
    ValueAndDerivativeAtIndex( const IndexType & index,
                               double & val,
                               const VectorType & v1, const VectorType & v2,
                               double scale=1 ) const;
  VectorType
    ValueAndDerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                                         double & val,
                                         double scale=1 ) const;
  VectorType
    ValueAndDerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                                         double & val,
                                         const VectorType & v1,
                                         double scale=1 ) const;
  VectorType
    ValueAndDerivativeAtContinuousIndex( const ContinuousIndexType & cIndex,
                                         double & val,
                                         const VectorType & v1, 
                                         const VectorType & v2,
                                         double scale=1 ) const;

  double
    Jet( const PointType & point, VectorType & d, MatrixType & h,
         double scale=1 ) const;
  double
    Jet( const PointType & point, VectorType & d, MatrixType & h,
         const VectorType & v1, double scale=1 ) const;
  double
    Jet( const PointType & point, VectorType & d, MatrixType & h,
         const VectorType & v1, const VectorType & v2, double scale=1 ) const;

  double
    JetAtIndex( const IndexType & cIndex, VectorType & d, MatrixType & h,
                double scale=1 ) const;
  double
    JetAtIndex( const IndexType & cIndex, VectorType & d, MatrixType & h,
                const VectorType & v1, double scale=1 ) const;
  double
    JetAtIndex( const IndexType & cIndex, VectorType & d, MatrixType & h,
                const VectorType & v1, const VectorType & v2, 
                double scale=1 ) const;

  double
    JetAtContinuousIndex( const ContinuousIndexType & cIndex,
                          VectorType & d, MatrixType & h,
                          double scale=1 ) const;
  double
    JetAtContinuousIndex( const ContinuousIndexType & cIndex,
                          VectorType & d, MatrixType & h,
                          const VectorType & v1, double scale=1 ) const;
  double
    JetAtContinuousIndex( const ContinuousIndexType & cIndex,
                          VectorType & d, MatrixType & h,
                          const VectorType & v1, const VectorType & v2,
                          double scale=1 ) const;

  double
    Ridgeness( const PointType & point, double scale=1 ) const;
  double
    Ridgeness( const PointType & point, 
               const VectorType & v1, double scale=1 ) const;
  double
    Ridgeness( const PointType & point, 
               const VectorType & v1, const VectorType & v2, 
               double scale=1 ) const;

  double
    RidgenessAtIndex( const IndexType & index, double scale=1 ) const;
  double
    RidgenessAtIndex( const IndexType & index, 
                      const VectorType & v1, double scale=1 ) const;
  double
    RidgenessAtIndex( const IndexType & index, 
                      const VectorType & v1, const VectorType & v2,
                      double scale=1 ) const;

  double
    RidgenessAtContinuousIndex( const ContinuousIndexType & cIndex,
                                double scale=1 ) const;
  double
    RidgenessAtContinuousIndex( const ContinuousIndexType & cIndex,
                                const VectorType & v1, double scale=1 ) const;
  double
    RidgenessAtContinuousIndex( const ContinuousIndexType & cIndex,
                                const VectorType & v1, const VectorType & v2,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivative(const PointType &point, double & val,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivative(const PointType &point, double & val,
                                const VectorType & v1, 
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivative(const PointType &point, double & val,
                                const VectorType & v1, const VectorType & v2,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivativeAtIndex(const IndexType &index, double & val,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivativeAtIndex(const IndexType &index, double & val,
                                const VectorType & v1,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivativeAtIndex(const IndexType &index, double & val,
                                const VectorType & v1, const VectorType & v2,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivativeAtContinuousIndex(const ContinuousIndexType &cIndex,
                                double & val,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivativeAtContinuousIndex(const ContinuousIndexType &cIndex,
                                double & val,
                                const VectorType & v1,
                                double scale=1 ) const;
  VectorType
    RidgenessAndDerivativeAtContinuousIndex(const ContinuousIndexType &cIndex,
                                double & val,
                                const VectorType & v1, const VectorType & v2,
                                double scale=1 ) const;

  MatrixType
    Hessian( const PointType & point, double scale=1 ) const;
  MatrixType
    Hessian( const PointType & point, 
             const VectorType & v1, double scale=1 ) const;
  MatrixType
    Hessian( const PointType & point, 
             const VectorType & v1, const VectorType & v2, 
             double scale=1 ) const;

  MatrixType
    HessianAtIndex( const IndexType & index, double scale=1 ) const;
  MatrixType
    HessianAtIndex( const IndexType & index, 
                    const VectorType & v1, double scale=1 ) const;
  MatrixType
    HessianAtIndex( const IndexType & index, 
                    const VectorType & v1, const VectorType & v2, 
                    double scale=1 ) const;

  MatrixType
    HessianAtContinuousIndex( const ContinuousIndexType & cIndex,
                              double scale=1 ) const;
  MatrixType
    HessianAtContinuousIndex( const ContinuousIndexType & cIndex,
                              const VectorType & v1, double scale=1 ) const;
  MatrixType
    HessianAtContinuousIndex( const ContinuousIndexType & cIndex,
                              const VectorType & v1, const VectorType & v2,
                              double scale=1 ) const;

  InputImagePointer ScaleSubsample(double factor);

  /**
   * Set the Extent
   */
  itkSetMacro( Extent, double);

  /**
   * Get the Extent
   */
  itkGetMacro( Extent, double);

  /**
   * Get the Spacing
   */
  itkGetConstMacro( InputImageSpacing, SpacingType );

  /**
   * If set to true, then values and derivatives are computed and then
   *   projected onto the subset of directions given.   Otherwise, the
   *   computations are done explicitly within the vector/plane given.
   */
  itkSetMacro( UseProjection, bool);
  itkGetMacro( UseProjection, bool);

protected:
  NJetImageFunction();
  NJetImageFunction( const Self& ){};

  ~NJetImageFunction(){};

  void operator=( const Self& ){};
  void PrintSelf(std::ostream& os, Indent indent) const;

  typename InputImageType::ConstPointer  m_InputImage;
  typename InputImageType::ConstPointer  m_InputImageMask;
  bool                                   m_UseInputImageMask;

  SizeType                m_InputImageSize;
  SpacingType             m_InputImageSpacing;
  SpacingType             m_InputImageSpacingSquared;
  double                  m_Extent;

  bool                    m_ValidStats;
  double                  m_StatsMin;
  double                  m_StatsMax;

  bool                    m_UseProjection;

};
  
} // namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
  #include "itkNJetImageFunction.txx"
#endif
  
#endif
