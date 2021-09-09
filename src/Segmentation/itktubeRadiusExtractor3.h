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

#ifndef __itktubeRadiusExtractor3_h
#define __itktubeRadiusExtractor3_h

#include "itktubeBlurImageFunction.h"

#include <itkTubeSpatialObject.h>

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
class RadiusExtractor3 : public Object
{
public:

  /**
   * Standard self typedef */
  typedef RadiusExtractor3                                   Self;
  typedef Object                                             Superclass;
  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  itkTypeMacro( RadiusExtractor3, Object );
  itkNewMacro( RadiusExtractor3 );

  /**
   * Standard for the number of dimension
   */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  typedef TubeSpatialObject< TInputImage::ImageDimension > TubeType;

  typedef typename TubeType::TubePointType                   TubePointType;

  typedef typename TubeType::PointType                       PointType;
  typedef typename TubeType::VectorType                      VectorType;

  /**
   * Type definition for the input image. */
  typedef TInputImage                                        InputImageType;

  typedef typename InputImageType::IndexType                 IndexType;

  /**
   * Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType                    PixelType;

  /**
   * Defines the type of matrix used
   */
  typedef vnl_matrix< double >                               MatrixType;


  typedef enum { RADIUS_CORRECTION_NONE, RADIUS_CORRECTION_FOR_BINARY_IMAGE,
    RADIUS_CORRECTION_FOR_CTA, RADIUS_CORRECTION_FOR_MRA }
    RadiusCorrectionFunctionType;

  /**
   * Set the input image */
  void SetInputImage( typename InputImageType::Pointer inputImage );

  /**
   * Get the input image */
  itkGetConstObjectMacro( InputImage, InputImageType );

  /** Set Data Minimum */
  itkSetMacro( DataMin, double );
  itkGetMacro( DataMin, double );

  /** Set Data Maximum */
  itkSetMacro( DataMax, double );
  itkGetMacro( DataMax, double );

  /** Set Minimum Radius */
  itkSetMacro( RadiusMinInIndexSpace, double );
  itkGetMacro( RadiusMinInIndexSpace, double );
  void SetRadiusMin( double r )
    { this->SetRadiusMinInIndexSpace( r / m_Spacing ); }
  double GetRadiusMin()
    { return this->GetRadiusMinInIndexSpace() * m_Spacing; }

  /** Set Maximum Radius */
  itkSetMacro( RadiusMaxInIndexSpace, double );
  itkGetMacro( RadiusMaxInIndexSpace, double );
  void SetRadiusMax( double r )
    { this->SetRadiusMaxInIndexSpace( r / m_Spacing ); }
  double GetRadiusMax()
    { return this->GetRadiusMaxInIndexSpace() * m_Spacing; }

  /** Set Radius0 */
  itkSetMacro( RadiusStartInIndexSpace, double );
  itkGetMacro( RadiusStartInIndexSpace, double );
  void SetRadiusStart( double r )
    { this->SetRadiusStartInIndexSpace( r / m_Spacing ); }
  double GetRadiusStart()
    { return this->GetRadiusStartInIndexSpace() * m_Spacing; }

  /** Set Radius step size when searching */
  itkSetMacro( RadiusStepInIndexSpace, double );
  itkGetMacro( RadiusStepInIndexSpace, double );
  void SetRadiusStep( double r )
    { this->SetRadiusStepInIndexSpace( r / m_Spacing ); }
  double GetRadiusStep()
    { return this->GetRadiusStepInIndexSpace() * m_Spacing; }

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
    double rStep );

  void SetKernelNumberOfPoints( unsigned int _numPoints );
  itkGetMacro( KernelNumberOfPoints, unsigned int );

  itkGetMacro( KernelPointStep, unsigned int );
  itkSetMacro( KernelPointStep, unsigned int );

  itkGetMacro( KernelStep, unsigned int );
  itkSetMacro( KernelStep, unsigned int );

  /** Calculate Radii */
  bool ExtractRadii( TubeType * tube, bool verbose=false );

  void SetIdleCallBack( bool ( *idleCallBack )( void ) );
  void SetStatusCallBack( void ( *statusCallBack )( const char *,
      const char *, int ) );

  double GetKernelMedialness( double r );

protected:

  RadiusExtractor3( void );
  virtual ~RadiusExtractor3( void );

  void GenerateKernelProfile( void );

  void SetKernelTubePoints( const std::vector< TubePointType > & tubePoints );
  std::vector< TubePointType > & GetKernelTubePoints( void )
   { return m_KernelTube->GetPoints(); };

  double GetProfileMaxDistance();
  double GetProfileBinNumber( double x );
  double GetProfileBinMaxRadius( double i );

  itkGetMacro( ProfileBinCount, std::vector< double > );
  itkGetMacro( ProfileBinValue, std::vector< double > );

  double GetKernelBranchness( double r );

  bool UpdateKernelOptimalRadius( void );
  itkGetMacro( KernelOptimalRadius, double );
  itkGetMacro( KernelOptimalRadiusMedialness, double );
  itkGetMacro( KernelOptimalRadiusBranchness, double );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void GenerateKernelTubePoints( unsigned int tubePointNum,
    TubeType * tube );

  void RecordOptimaAtTubePoints( unsigned int tubePointNum,
    TubeType * tube );

private:

  RadiusExtractor3( const Self& );
  void operator=( const Self& );

  typename InputImageType::Pointer        m_InputImage;
  double                                  m_Spacing;
  double                                  m_DataMin;
  double                                  m_DataMax;

  double                                  m_RadiusStartInIndexSpace;
  double                                  m_RadiusMinInIndexSpace;
  double                                  m_RadiusMaxInIndexSpace;
  double                                  m_RadiusStepInIndexSpace;

  double                                  m_RadiusCorrectionScale;
  RadiusCorrectionFunctionType            m_RadiusCorrectionFunction;

  double                                  m_MinMedialness;
  double                                  m_MinMedialnessStart;

  typename TubeType::Pointer              m_KernelTube;
  unsigned int                            m_KernelNumberOfPoints;

  unsigned int                            m_KernelPointStep;
  unsigned int                            m_KernelStep;

  unsigned int                            m_ProfileNumberOfBins;
  std::vector< double >                   m_ProfileBinCount;
  std::vector< double >                   m_ProfileBinValue;

  double                                  m_KernelOptimalRadius;
  double                                  m_KernelOptimalRadiusMedialness;
  double                                  m_KernelOptimalRadiusBranchness;

  void ( * m_StatusCallBack )( const char *, const char *, int );
  bool ( * m_IdleCallBack )( void );

}; // End class RadiusExtractor3

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRadiusExtractor3.hxx"
#endif

#endif // End !defined( __itktubeRadiusExtractor3_h )
