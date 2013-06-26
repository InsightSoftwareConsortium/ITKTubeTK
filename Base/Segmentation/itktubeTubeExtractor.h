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

#ifndef __itktubeTubeExtractor_h
#define __itktubeTubeExtractor_h

#include "itktubeRadiusExtractor.h"
#include "itktubeRidgeExtractor.h"

#include <itkObject.h>

namespace itk
{

namespace tube
{

/**
 * This class extract the a tube given an image
 *
 * \sa TubeExtractor
 */

template< class TInputImage >
class TubeExtractor : public Object
{
public:

  /**
   * Standard self typedef */
  typedef TubeExtractor             Self;
  typedef Object                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /**
   * Run-time type information ( and related methods ). */
  itkTypeMacro( TubeExtractor, Object );

  itkNewMacro( TubeExtractor );

  /**
   * Type definition for the input image. */
  typedef TInputImage                                   ImageType;

  typedef typename RidgeExtractor<ImageType>::TubeMaskImageType
                                                        TubeMaskImageType;

  /**
   * Standard for the number of dimension
   */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /**
   * Type definition for the input image pixel type. */
  typedef typename ImageType::PixelType                 PixelType;

  /**  Type definition for VesselTubeSpatialObject */
  typedef VesselTubeSpatialObject< ImageDimension >     TubeType;
  typedef typename TubeType::TubePointType              TubePointType;

  typedef RidgeExtractor<ImageType>                     RidgeOpType;
  typedef RadiusExtractor<ImageType>                    RadiusOpType;

  /**
   * Type definition for the input image pixel type. */
  typedef ContinuousIndex<double, ImageDimension >      ContinuousIndexType;

  typedef typename ImageType::IndexType                 IndexType;

  /**
   * Defines the type of vectors used */
  typedef itk::Vector<double,  ImageDimension >         VectorType;


  /**
   * Set the input image */
  void SetInputImage( typename ImageType::Pointer inputImage );

  /**
   * Get the input image */
  itkGetConstObjectMacro( InputImage, ImageType );

  /**
   * Get the input image */
  typename TubeMaskImageType::Pointer GetTubeMaskImage( void );

  /**
   * Set Data Minimum */
  void SetDataMin( double dataMin );

  /**
   * Get Data Minimum */
  double GetDataMin( void );

  /**
   * Set Data Maximum */
  void SetDataMax( double dataMax );

  /**
   * Get Data Maximum */
  double GetDataMax( void );

  /**
   * Set ExtractBound Minimum */
  void SetExtractBoundMin( IndexType dataMin );

  /**
   * Get ExtractBound Minimum */
  IndexType GetExtractBoundMin( void ) const;

  /**
   * Set ExtractBound Maximum */
  void SetExtractBoundMax( IndexType dataMax );

  /**
   * Get ExtractBound Maximum */
  IndexType GetExtractBoundMax( void ) const;

  /**
   * Set the radius */
  void SetRadius( double radius );

  /**
   * Get the radius */
  double GetRadius( void );

  /**
   * Set Extract Ridge */
  void ExtractBrightTube( bool extractRidge );

  /**
   * Get Extract Ridge */
  bool ExtractBrightTube( void );

  /**
   * Get the ridge extractor */
  typename RidgeExtractor<ImageType>::Pointer GetRidgeOp( void );

  /**
   * Get the radius extractor */
  typename RadiusExtractor<ImageType>::Pointer GetRadiusOp( void );

  /**
   * Return true if a tube is found from the given seed point */
  bool LocalTube( ContinuousIndexType & x );

  /**
   * Extract the ND tube given the position of the first point
   * and the tube ID */
  typename TubeType::Pointer ExtractTube( ContinuousIndexType & x,
    unsigned int tubeID );

  /**
   * Delete a tube */
  void SmoothTube( TubeType * tube, int h=5 );

  /**
   * Delete a tube */
  bool AddTube( TubeType * tube );

  /**
   * Delete a tube */
  bool DeleteTube( TubeType * tube );

  /**
   * Set the tube color */
  void SetColor( float color[4] );

  /**
   * Get the tube color */
  itkGetMacro( Color, float * );

  /**
   * Set the idle callback */
  void   IdleCallBack( bool ( *idleCallBack )( void ) );

  /**
   * Set the status callback */
  void   StatusCallBack( void ( *statusCallBack )( const char *,
      const char *, int ) );

  /**
   * Set the tube callback */
  void   NewTubeCallBack( void ( *newTubeCallBack )( TubeType * ) );

  /**
   * Set the status callback */
  void   AbortProcess( bool ( *abortProcess )( void ) );

protected:

  TubeExtractor( void );
  virtual ~TubeExtractor( void );
  TubeExtractor( const Self& ) {}
  void operator=( const Self& ) {}

  void PrintSelf( std::ostream & os, Indent indent ) const;

  typename RidgeExtractor<ImageType>::Pointer  m_RidgeOp;
  typename RadiusExtractor<ImageType>::Pointer m_RadiusOp;

  bool ( *m_IdleCallBack )( void );
  void ( *m_StatusCallBack )( const char *, const char *, int );
  void ( *m_NewTubeCallBack )( TubeType * );
  bool ( *m_AbortProcess )( void );

private:

  typename ImageType::Pointer  m_InputImage;
  float                        m_Color[4];

}; // End class TubeExtractor

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeExtractor.hxx"
#endif

#endif // End !defined(__itktubeTubeExtractor_h)
