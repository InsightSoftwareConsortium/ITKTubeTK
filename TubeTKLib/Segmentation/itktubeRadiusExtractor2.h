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

#ifndef __itktubeRadiusExtractor2_h
#define __itktubeRadiusExtractor2_h

#include "itktubeBlurImageFunction.h"

#include <itkVesselTubeSpatialObject.h>

#include <itkMath.h>

#include <vector>

namespace itk
{

namespace tube
{

/**
 * This class extract the radius of a tube given an image
 * \sa RidgeExtractor
 */

template< class TInputImage >
class RadiusExtractor2 : public Object
{
public:

  /**
   * Standard self typedef */
  typedef RadiusExtractor2                                   Self;
  typedef Object                                             Superclass;
  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  itkTypeMacro( RadiusExtractor2, Object );
  itkNewMacro( RadiusExtractor2 );

  /**
   * Standard for the number of dimension
   */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  typedef VesselTubeSpatialObject< TInputImage::ImageDimension > TubeType;

  typedef typename TubeType::TubePointType                   TubePointType;

  typedef typename TubeType::PointType                       ITKPointType;
  typedef typename TubeType::VectorType                      ITKVectorType;

  /**
   * Type definition for the input image. */
  typedef TInputImage                                        ImageType;

  typedef typename ImageType::IndexType                      ITKIndexType;

  /**
   * Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType                    PixelType;

  /**
   * Defines the type of vectors used
   */
  typedef vnl_vector< double >                               VectorType;

  /**
   * Defines the type of matrix used
   */
  typedef vnl_matrix< double >                               MatrixType;


  typedef enum { RADIUS_CORRECTION_NONE, RADIUS_CORRECTION_FOR_BINARY_IMAGE,
    RADIUS_CORRECTION_FOR_CTA, RADIUS_CORRECTION_FOR_MRA }
    RadiusCorrectionFunctionType;

  /**
   * Set the input image */
  void SetInputImage( typename ImageType::Pointer inputImage );

  /**
   * Get the input image */
  itkGetConstObjectMacro( Image, ImageType );

  /** Set Data Minimum */
  itkSetMacro( DataMin, double );
  itkGetMacro( DataMin, double );

  /** Set Data Maximum */
  itkSetMacro( DataMax, double );
  itkGetMacro( DataMax, double );

  /** Set Minimum Radius */
  itkSetMacro( RadiusMin, double );
  itkGetMacro( RadiusMin, double );

  /** Set Maximum Radius */
  itkSetMacro( RadiusMax, double );
  itkGetMacro( RadiusMax, double );

  /** Set Radius0 */
  itkSetMacro( RadiusStart, double );
  itkGetMacro( RadiusStart, double );

  /** Set Radius step size when searching */
  itkSetMacro( RadiusStep, double );
  itkGetMacro( RadiusStep, double );

  /** Set Radius tolerance when searching */
  itkSetMacro( RadiusTolerance, double );
  itkGetMacro( RadiusTolerance, double );

  /** Set Radius Correction Scale - multiply optimal radius to get the
   * modality-specific radius */
  itkSetMacro( RadiusCorrectionScale, double );
  itkGetMacro( RadiusCorrectionScale, double );

  itkSetMacro( RadiusCorrectionFunction, RadiusCorrectionFunctionType );
  itkGetMacro( RadiusCorrectionFunction, RadiusCorrectionFunctionType );

  /** Set ThreshMedialness */
  itkSetMacro( MinMedialness, double );
  itkGetMacro( MinMedialness, double );

  /** Set ThreshMedialness Start */
  itkSetMacro( MinMedialnessStart, double );
  itkGetMacro( MinMedialnessStart, double );

  void GetPointVectorMeasures( std::vector< TubePointType > & points,
    double pntR,
    double & mness,
    double & bness,
    bool doBNess );

  /** Calculate the optimal scale */
  bool GetPointVectorOptimalRadius( std::vector< TubePointType > & points,
    double & r0,
    double rMin,
    double rMax,
    double rStep,
    double rTolerance );

  void SetNumKernelPoints( unsigned int _numPoints );
  itkGetMacro( NumKernelPoints, unsigned int );

  itkGetMacro( KernelPointStep, unsigned int );
  itkSetMacro( KernelPointStep, unsigned int );

  itkGetMacro( KernelStep, unsigned int );
  itkSetMacro( KernelStep, unsigned int );

  itkGetMacro( KernelExtent, double );
  itkSetMacro( KernelExtent, double );

  void GenerateKernel( void );

  void SetKernelTubePoints( const std::vector< TubePointType > & tubePoints );
  itkGetMacro( KernelTubePoints, std::vector< TubePointType > );

  itkGetMacro( KernelValues, std::vector< double > );
  itkGetMacro( KernelDistances, std::vector< double > );
  itkGetMacro( KernelTangentDistances, std::vector< double > );

  double GetKernelMedialness( double r );
  double GetKernelBranchness( double r );

  bool UpdateKernelOptimalRadius( void );
  itkGetMacro( KernelOptimalRadius, double );
  itkGetMacro( KernelOptimalRadiusMedialness, double );
  itkGetMacro( KernelOptimalRadiusBranchness, double );

  /** Calculate Radii */
  bool ExtractRadii( TubeType * tube );

  void SetIdleCallBack( bool ( *idleCallBack )( void ) );
  void SetStatusCallBack( void ( *statusCallBack )( const char *,
      const char *, int ) );

protected:

  RadiusExtractor2( void );
  virtual ~RadiusExtractor2( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

  void GenerateKernelTubePoints( unsigned int tubePointNum,
    TubeType * tube );

  void RecordOptimaAtTubePoints( unsigned int tubePointNum,
    TubeType * tube );

private:

  RadiusExtractor2( const Self& );
  void operator=( const Self& );

  typename ImageType::Pointer             m_Image;
  double                                  m_DataMin;
  double                                  m_DataMax;

  double                                  m_RadiusStart;
  double                                  m_RadiusMin;
  double                                  m_RadiusMax;
  double                                  m_RadiusStep;
  double                                  m_RadiusTolerance;

  double                                  m_RadiusCorrectionScale;
  RadiusCorrectionFunctionType            m_RadiusCorrectionFunction;

  double                                  m_MinMedialness;
  double                                  m_MinMedialnessStart;

  unsigned int                            m_NumKernelPoints;
  std::vector< TubePointType >            m_KernelTubePoints;

  unsigned int                            m_KernelPointStep;
  unsigned int                            m_KernelStep;
  double                                  m_KernelExtent;

  std::vector< double >                   m_KernelValues;
  std::vector< double >                   m_KernelDistances;
  std::vector< double >                   m_KernelTangentDistances;

  double                                  m_KernelOptimalRadius;
  double                                  m_KernelOptimalRadiusMedialness;
  double                                  m_KernelOptimalRadiusBranchness;

  void ( * m_StatusCallBack )( const char *, const char *, int );
  bool ( * m_IdleCallBack )( void );

}; // End class RadiusExtractor2

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRadiusExtractor2.hxx"
#endif

#endif // End !defined( __itktubeRadiusExtractor2_h )
