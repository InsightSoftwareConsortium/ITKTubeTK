/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeRidgeExtractor_h
#define __itktubeRidgeExtractor_h

#include "itktubeBlurImageFunction.h"
#include "itktubeRadiusExtractor.h"
#include "tubeBrentOptimizer1D.h"
#include "tubeSplineApproximation1D.h"
#include "tubeSplineND.h"

#include <itkContinuousIndex.h>
#include <itkVesselTubeSpatialObject.h>

#include <cmath>
#include <list>

namespace itk
{

namespace tube
{

/**
 * This class extract the ridge of a tube given an image
 *
 * \sa RidgeExtractor
 */

template< class TInputImage >
class RidgeExtractor : public Object
{
public:

  typedef RidgeExtractor             Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkTypeMacro( RidgeExtractor, Object );

  itkNewMacro( RidgeExtractor );

  /**
   * Type definition for the input image. */
  typedef TInputImage                                      ImageType;

  /**
   * Standard for the number of dimension */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /**
   * Type definition for the input image. */
  typedef Image< float, TInputImage::ImageDimension >     TubeMaskImageType;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType                 PixelType;

  /** Type definition for the input image index type. */
  typedef typename TInputImage::IndexType                 IndexType;

  /** Type definition for the input image point type. */
  typedef typename TInputImage::PointType                 PointType;

  /** Type definition for the input image index type. */
  typedef ContinuousIndex< double, TInputImage::ImageDimension >
                                                ContinuousIndexType;

  /** Defines the type of vectors used */
  typedef vnl_vector< double >                  VectorType;


  /** Defines the type of matrix used */
  typedef vnl_matrix< double >                  MatrixType;

  /** Tube SpatialObject typedefs */
  typedef VesselTubeSpatialObject< TInputImage::ImageDimension >
                                                TubeType;
  typedef typename TubeType::TubePointType      TubePointType;

  /** Defines the type of vectors used */
  typedef typename TubeType::CovariantVectorType
                                                CovariantVectorType;

  /**
   * Set the input image */
  void SetInputImage( typename ImageType::Pointer inputImage );

  /**
   * Get the input image */
  typename ImageType::Pointer GetInputImage( void );

  /**
   * Get the mask image */
  itkGetObjectMacro( TubeMaskImage, TubeMaskImageType );

  /*
   * Set Data Minimum */
  void SetDataMin( double dataMin );

  /**
   * Get Data Minimum */
  itkGetMacro( DataMin, double );

  /**
   * Set Data Maximum */
  void SetDataMax( double dataMax );

  /**
   * Get Data Maximum */
  itkGetMacro( DataMax, double );

  /**
   * Set Traversal Step size */
  itkSetMacro( StepX, double );

  /**
   * Get Traversal Step size */
  itkGetMacro( StepX, double );

  /**
   * Set Tangent change threshold */
  itkSetMacro( ThreshT, double );

  /**
   * Get Tangent change threshold */
  itkGetMacro( ThreshT, double );

  /**
   * Set Spatial change threshold */
  itkSetMacro( ThreshX, double );

  /**
   * Get Spatial change threshold */
  itkGetMacro( ThreshX, double );

  /**
   * Set Ridgeness Threshold */
  itkSetMacro( ThreshRidgeness, double );

  /**
   * Get Ridgeness Threshold */
  itkGetMacro( ThreshRidgeness, double );

  /**
   * Set Ridgeness start Threshold */
  itkSetMacro( ThreshRidgenessStart, double );

  /**
   * Get Ridgeness start Threshold */
  itkGetMacro( ThreshRidgenessStart, double );

  /**
   * Set Roundness Threshold */
  itkSetMacro( ThreshRoundness, double );

  /**
   * Get Roundness Threshold */
  itkGetMacro( ThreshRoundness, double );

  /**
   * Set Roundness start Threshold */
  itkSetMacro( ThreshRoundnessStart, double );

  /**
   * Get Roundness start Threshold */
  itkGetMacro( ThreshRoundnessStart, double );

  /**
   * Set Curvature  Normalization Factor (multiplier for curvature) */
  itkSetMacro( CurvatureExpectedMax, double );

  /**
   * Get Curvature  Normalization Factor */
  itkGetMacro( CurvatureExpectedMax, double );

  /**
   * Set Curvature  Threshold */
  itkSetMacro( ThreshCurvature, double );

  /**
   * Get Curvature  Threshold */
  itkGetMacro( ThreshCurvature, double );

  /**
   * Set Curvature  Threshold */
  itkSetMacro( ThreshCurvatureStart, double );

  /**
   * Get Curvature  Threshold */
  itkGetMacro( ThreshCurvatureStart, double );

  /**
   * Set Linearity  Threshold */
  itkSetMacro( ThreshLinearity, double );

  /**
   * Get Linearity  Threshold */
  itkGetMacro( ThreshLinearity, double );

  /**
   * Set Linearity  Threshold */
  itkSetMacro( ThreshLinearityStart, double );

  /**
   * Get Linearity  Threshold */
  itkGetMacro( ThreshLinearityStart, double );


  /**
   * Set Extract Bound Minimum */
  itkSetMacro( ExtractBoundMin, IndexType );

  /**
   * get Extract Bound Minimum */
  itkGetMacro( ExtractBoundMin, IndexType );

  /**
   * Set Extract Bound Maximum */
  itkSetMacro( ExtractBoundMax, IndexType );

  /**
   * get Extract Bound Maximum */
  itkGetMacro( ExtractBoundMax, IndexType );

  /**
   * Get the data spline */
  ::tube::SplineND * GetDataSpline( void );

  /**
   * Get the data spline 1D */
  ::tube::Spline1D * GetDataSpline1D( void );

  /**
   * Get the data spline optimizer */
  ::tube::Optimizer1D * GetDataSplineOptimizer( void );

  /**
   * Set the scale */
  void SetScale( double scale );

  /**
   * Get the scale */
  double GetScale( void );

  /**
   * Set the extent */
  void SetExtent( double extent );

  /**
   * Get the extent */
  double GetExtent( void );

  /**
   * Set the dynamic Scale */
  void SetDynamicScale( bool dynamicScale );

  /**
   * Get the dynamicScale */
  itkGetMacro( DynamicScale, bool );

  itkGetMacro( DynamicScaleUsed, double );

  /**
   * Set the Recovery Maximum */
  itkSetMacro( RecoveryMax, int );

  /**
   * Get the Recovery Maximum */
  itkGetMacro( RecoveryMax, int );

  /**
   * Delete a tube */
  template< class TDrawMask >
  bool DeleteTube( TubeType * tube, TDrawMask * drawMask );
  bool DeleteTube( TubeType * tube );

  /**
   * Add a tube */
  template< class TDrawMask >
  bool AddTube( TubeType * tube, TDrawMask * drawMask );
  bool AddTube( TubeType * tube );

  /**
   * Smooth a tube */
  void SmoothTube( TubeType * tube, int h );

  /**
   * Set the radius Extractor */
  void  SetRadiusExtractor( RadiusExtractor<TInputImage> * radiusExtractor );

  /**
   * Compute the intensity at the point x */
  double  Intensity( const IndexType & x );

  /**the ridgeness at point x */
  double  Ridgeness( const ContinuousIndexType & x, double & roundness,
    double & curvature, double & linearity );

  /**
   * Compute the local Ridge */
  bool   LocalRidge( ContinuousIndexType & x );

  /**
   * Extract */
  typename TubeType::Pointer  ExtractRidge( ContinuousIndexType & x,
    int tubeID );

  /**
   * Set the idle callback */
  void   IdleCallBack( bool ( *idleCallBack )( void ) );

  /**
   * Set the status callback */
  void   StatusCallBack( void ( *statusCallBack )( const char *,
      const char *, int ) );

protected:

  RidgeExtractor( void );
  virtual ~RidgeExtractor( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

  /**
   * Traverse the ridge one way */
  bool  TraverseOneWay( ContinuousIndexType & newX, VectorType & newT,
    MatrixType & newN, int dir );

private:

  RidgeExtractor( const Self& );
  void operator=( const Self& );

  typename ImageType::Pointer                        m_InputImage;

  typename BlurImageFunction<ImageType>::Pointer     m_DataFunc;

  typename TubeMaskImageType::Pointer                m_TubeMaskImage;

  bool                                               m_DynamicScale;
  double                                             m_DynamicScaleUsed;
  RadiusExtractor<TInputImage>                     * m_RadiusExtractor;

  int                                                m_RecoveryMax;

  double                                             m_DataMin;
  double                                             m_DataMax;
  double                                             m_DataRange;

  double                                             m_StepX;
  double                                             m_ThreshT;
  double                                             m_ThreshX;

  IndexType                                          m_ExtractBoundMin;
  IndexType                                          m_ExtractBoundMax;

  ::tube::SplineApproximation1D                      m_DataSpline1D;
  ::tube::BrentOptimizer1D                           m_DataSplineOpt;
  ::tube::SplineND                                 * m_DataSpline;
  ::tube::UserFunction< vnl_vector<int>, double >  * m_SplineValueFunc;

  double                                             m_CurvatureExpectedMax;

  double                                             m_ThreshRidgeness;
  double                                             m_ThreshRidgenessStart;
  double                                             m_ThreshRoundness;
  double                                             m_ThreshRoundnessStart;
  double                                             m_ThreshCurvature;
  double                                             m_ThreshCurvatureStart;
  double                                             m_ThreshLinearity;
  double                                             m_ThreshLinearityStart;

  VectorType                                         m_X;
  VectorType                                         m_XP;
  double                                             m_XVal;

  VectorType                                         m_XD;
  MatrixType                                         m_XH;
  VectorType                                         m_XHEVal;
  MatrixType                                         m_XHEVect;

  typename TubeType::Pointer                         m_Tube;

  bool  ( *m_IdleCallBack )( void );
  void  ( *m_StatusCallBack )( const char *, const char *, int );

}; // End class RidgeExtractor

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeExtractor.hxx"
#endif

#endif // End !defined(__itktubeRidgeExtractor_h)
