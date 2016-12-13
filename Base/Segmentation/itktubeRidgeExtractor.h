/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

#ifndef __itktubeRidgeExtractor_h
#define __itktubeRidgeExtractor_h

#include "itktubeBlurImageFunction.h"
#include "itktubeRadiusExtractor2.h"
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

  /** Type definition for the input image. */
  typedef TInputImage                                      ImageType;

  /** Standard for the number of dimension */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Type definition for the input image. */
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
  typedef vnl_vector< unsigned int >            IntVectorType;

  /** Defines the type of vectors used */
  typedef vnl_vector< double >                  VectorType;

  /** Defines the type of matrix used */
  typedef vnl_matrix< double >                  MatrixType;

  /** Tube SpatialObject typedefs */
  typedef VesselTubeSpatialObject< TInputImage::ImageDimension > TubeType;

  typedef typename TubeType::TubePointType      TubePointType;

  /** Defines the type of vectors used */
  typedef typename TubeType::CovariantVectorType CovariantVectorType;

  typedef enum { SUCCESS, EXITED_IMAGE, REVISITED_VOXEL, RIDGE_FAIL,
    ROUND_FAIL, CURVE_FAIL, LEVEL_FAIL, TANGENT_FAIL, DISTANCE_FAIL,
    OTHER_FAIL }                                FailureCodeEnum;

  /** Set the input image */
  void SetInputImage( typename ImageType::Pointer inputImage );

  /** Get the input image */
  typename ImageType::Pointer GetInputImage( void );

  /** Get the mask image */
  itkGetObjectMacro( TubeMaskImage, TubeMaskImageType );

  /** Set the mask image */
  itkSetObjectMacro( TubeMaskImage, TubeMaskImageType );

  /** Set Data Minimum */
  void SetDataMin( double dataMin );

  /** Get Data Minimum */
  itkGetMacro( DataMin, double );

  /** Set Data Maximum */
  void SetDataMax( double dataMax );

  /** Get Data Maximum */
  itkGetMacro( DataMax, double );

  /** Set Traversal Step size */
  itkSetMacro( StepX, double );

  /** Get Traversal Step size */
  itkGetMacro( StepX, double );

  /** Set Tangent change threshold */
  itkSetMacro( MaxTangentChange, double );

  /** Get Tangent change threshold */
  itkGetMacro( MaxTangentChange, double );

  /** Set Spatial change threshold */
  itkSetMacro( MaxXChange, double );

  /** Get Spatial change threshold */
  itkGetMacro( MaxXChange, double );

  /** Set Ridgeness Threshold */
  itkSetMacro( MinRidgeness, double );

  /** Get Ridgeness Threshold */
  itkGetMacro( MinRidgeness, double );

  /** Set Ridgeness start Threshold */
  itkSetMacro( MinRidgenessStart, double );

  /** Get Ridgeness start Threshold */
  itkGetMacro( MinRidgenessStart, double );

  /** Set Roundness Threshold */
  itkSetMacro( MinRoundness, double );

  /** Get Roundness Threshold */
  itkGetMacro( MinRoundness, double );

  /** Set Roundness start Threshold */
  itkSetMacro( MinRoundnessStart, double );

  /** Get Roundness start Threshold */
  itkGetMacro( MinRoundnessStart, double );

  /** Set Curvature  Threshold */
  itkSetMacro( MinCurvature, double );

  /** Get Curvature  Threshold */
  itkGetMacro( MinCurvature, double );

  /** Set Curvature  Threshold */
  itkSetMacro( MinCurvatureStart, double );

  /** Get Curvature  Threshold */
  itkGetMacro( MinCurvatureStart, double );

  /** Set Levelness  Threshold */
  itkSetMacro( MinLevelness, double );

  /** Get Levelness  Threshold */
  itkGetMacro( MinLevelness, double );

  /** Set Levelness  Threshold */
  itkSetMacro( MinLevelnessStart, double );

  /** Get Levelness  Threshold */
  itkGetMacro( MinLevelnessStart, double );


  /** Set Extract Bound Minimum */
  itkSetMacro( ExtractBoundMin, IndexType );

  /** Get Extract Bound Minimum */
  itkGetMacro( ExtractBoundMin, IndexType );

  /** Set Extract Bound Maximum */
  itkSetMacro( ExtractBoundMax, IndexType );

  /** Get Extract Bound Maximum */
  itkGetMacro( ExtractBoundMax, IndexType );

  /** Get the data spline */
  ::tube::SplineND * GetDataSpline( void );

  /** Get the data spline 1D */
  ::tube::Spline1D * GetDataSpline1D( void );

  /** Get the data spline optimizer */
  ::tube::Optimizer1D * GetDataSplineOptimizer( void );

  /** Set the scale */
  void SetScale( double scale );

  /** Get the scale */
  double GetScale( void );

  /** Set the extent */
  void SetScaleKernelExtent( double extent );

  /** Get the extent */
  double GetScaleKernelExtent( void );

  /** Set to re-estimate ( based on local radius estimate ) the scale to be
   * used for image measures made during ridge extraction */
  void SetDynamicScale( bool dynamicScale );

  /** Is the scale of ridge measures optimized during ridge traversal */
  itkGetMacro( DynamicScale, bool );

  /** Get the scale currently being used for ridge traversal measures. */
  itkGetMacro( DynamicScaleUsed, double );

  /** Set to dynamic update the ridge step size to be proportional to
   * estimated tube radius during ridge traversal. */
  void SetDynamicStepSize( bool dynamicStepSize );

  /** Are step sizes made along a ridge fixed, or do they change based on
   * local tube radius estimates. */
  itkGetMacro( DynamicStepSize, bool );

  /** Set the Recovery Maximum */
  itkSetMacro( MaxRecoveryAttempts, int );

  /** Get the Recovery Maximum */
  itkGetMacro( MaxRecoveryAttempts, int );

  /** Delete a tube */
  template< class TDrawMask >
  bool DeleteTube( const TubeType * tube, TDrawMask * drawMask );
  bool DeleteTube( const TubeType * tube );

  /** Add a tube */
  template< class TDrawMask >
  bool AddTube( const TubeType * tube, TDrawMask * drawMask );
  bool AddTube( const TubeType * tube );

  /** Set the radius Extractor */
  void  SetRadiusExtractor( RadiusExtractor2<TInputImage> * radiusExtractor );

  /** Compute the intensity at the point x */
  double  Intensity( const IndexType & x );

  /** Computes ridge measures at the given point x
   *  and stores them so that they can be queried later.
   *  Returns the ridgeness at x
   *  \param x User supplied point at which ridge measures will be computed
   *  \param intensity On return equals the interpolated intensity at x
   *    Can be queried later using GetCurrentIntensity()
   *  \param roundness On return equals the roundness at x
   *    Can be queried later using GetCurrentRoundness()
   *  \param curvature On return equals the curvature at x
   *    Can be queried later using GetCurrentCurvature()
   *  \param levelness On return equals the levelness at x
   *    Can be queried later using GetCurrentLevelness()
   */
  double  Ridgeness( const ContinuousIndexType & x,
    double & intensity,
    double & roundness,
    double & curvature,
    double & levelness,
    const vnl_vector<double> & prevTangent=vnl_vector<double>() );

  /** Get current location
   *  This is location at which the Ridgness function was last called
   */
  const VectorType & GetCurrentLocation() const;

  /** Get the Hessian Eigen Basis at the current location
   *  Each column of the return matrix is an eigen vector
   *  in increasing order of the eigen values.
   *  The third column is an approximation to the ridge tangent.
   */
  const MatrixType & GetCurrentBasis() const;

  /** Get intensity at the current location */
  double GetCurrentIntensity() const;

  /** Get ridgness at the current location */
  double GetCurrentRidgeness() const;

  /** Get roundess at the current location */
  double GetCurrentRoundness() const;

  /** Get curvature at the current location */
  double GetCurrentCurvature() const;

  /** Get levelness at the current location */
  double GetCurrentLevelness() const;

  /** Compute/find the local Ridge */
  FailureCodeEnum LocalRidge( ContinuousIndexType & x,
    bool verbose=false );

  /** Extract */
  typename TubeType::Pointer ExtractRidge( const ContinuousIndexType & x,
    int tubeID,
    bool verbose=false );

  itkGetMacro( CurrentFailureCode, FailureCodeEnum );

  unsigned int      GetNumberOfFailureCodes( void ) const;
  const std::string GetFailureCodeName( FailureCodeEnum code ) const;
  unsigned int      GetFailureCodeCount( FailureCodeEnum code ) const;
  void              ResetFailureCodeCounts( void );

  /** Set the idle callback */
  void   IdleCallBack( bool ( *idleCallBack )( void ) );

  /** Set the status callback */
  void   StatusCallBack( void ( *statusCallBack )( const char *,
      const char *, int ) );

protected:

  RidgeExtractor( void );
  virtual ~RidgeExtractor( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

  /** Traverse the ridge one way */
  bool  TraverseOneWay( ContinuousIndexType & newX, VectorType & newT,
    MatrixType & newN, int dir, bool verbose=false );

private:

  RidgeExtractor( const Self& );
  void operator=( const Self& );

  typename ImageType::Pointer                        m_InputImage;

  typename BlurImageFunction<ImageType>::Pointer     m_DataFunc;

  typename TubeMaskImageType::Pointer                m_TubeMaskImage;

  bool                                               m_DynamicScale;
  double                                             m_DynamicScaleUsed;
  bool                                               m_DynamicStepSize;
  RadiusExtractor2<TInputImage>                    * m_RadiusExtractor;

  int                                                m_MaxRecoveryAttempts;

  double                                             m_DataMin;
  double                                             m_DataMax;
  double                                             m_DataRange;

  double                                             m_StepX;
  double                                             m_MaxTangentChange;
  double                                             m_MaxXChange;

  IndexType                                          m_ExtractBoundMin;
  IndexType                                          m_ExtractBoundMax;

  ::tube::SplineApproximation1D                      m_DataSpline1D;
  ::tube::BrentOptimizer1D                           m_DataSplineOpt;
  ::tube::SplineND                                 * m_DataSpline;
  ::tube::UserFunction< vnl_vector<int>, double >  * m_SplineValueFunc;

  FailureCodeEnum                                    m_CurrentFailureCode;
  IntVectorType                                      m_FailureCodeCount;

  double                                             m_MinRidgeness;
  double                                             m_MinRidgenessStart;
  double                                             m_MinRoundness;
  double                                             m_MinRoundnessStart;
  double                                             m_MinCurvature;
  double                                             m_MinCurvatureStart;
  double                                             m_MinLevelness;
  double                                             m_MinLevelnessStart;

  // current location
  VectorType                                         m_X;

  VectorType                                         m_XP;

  // current intensity
  double                                             m_XVal;

  // current gradient
  VectorType                                         m_XD;

  // current Hessian
  MatrixType                                         m_XH;

  // current Hessian Eigen Values
  VectorType                                         m_XHEVal;

  // current Hessian Eigen Vectors
  MatrixType                                         m_XHEVect;

  // current ridgeness
  double                                             m_XRidgeness;

  // current roundness
  double                                             m_XRoundness;

  // current curvature
  double                                             m_XCurvature;

  // current levelness
  double                                             m_XLevelness;

  typename TubeType::Pointer                         m_Tube;

  bool  ( *m_IdleCallBack )( void );
  void  ( *m_StatusCallBack )( const char *, const char *, int );

}; // End class RidgeExtractor

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeExtractor.hxx"
#endif

#endif // End !defined( __itktubeRidgeExtractor_h )
