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

#ifndef __itktubeTubeExtractor_h
#define __itktubeTubeExtractor_h

#include "itktubeRadiusExtractor2.h"
#include "itktubeRidgeExtractor.h"

#include "itkGroupSpatialObject.h"

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
  typedef TubeExtractor               Self;
  typedef Object                      Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /**
   * Run-time type information ( and related methods ). */
  itkTypeMacro( TubeExtractor, Object );

  itkNewMacro( TubeExtractor );

  /**
   * Type definition for the input image. */
  typedef TInputImage                               ImageType;

  typedef RidgeExtractor<ImageType>                 RidgeExtractorType;
  typedef RadiusExtractor2<ImageType>               RadiusExtractorType;

  typedef typename RidgeExtractorType::TubeMaskImageType
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
  typedef TubeSpatialObject< ImageDimension >           TubeType;
  typedef typename TubeType::TubePointType              TubePointType;

  typedef itk::GroupSpatialObject< ImageDimension >     TubeGroupType;

  /**
   * Type definition for the input image pixel type. */
  typedef ContinuousIndex<double, ImageDimension >      ContinuousIndexType;
  typedef std::vector< ContinuousIndexType >            ContinuousIndexListType;
  typedef Point<double, ImageDimension >                PointType;
  typedef std::vector< PointType >                      PointListType;
  typedef double                                        RadiusType;
  typedef std::vector< RadiusType >                     RadiusListType;

  typedef typename ImageType::IndexType                 IndexType;

  /**
   * Defines the type of vectors used */
  typedef itk::Vector< double, ImageDimension >         VectorType;


  /***********/
  /***********/
  /***********/


  /**
   * Set the input image */
  void SetInputImage( ImageType * inputImage );
  const ImageType * GetInputImage( void ) const;

  /**
   * Optionally set a different input image to use for radius estimation */
  void SetRadiusInputImage( ImageType * radiusInputImage );
  const ImageType * GetRadiusInputImage( void ) const;

  /**
   * Set the tube mask image */
  void SetTubeMaskImage( TubeMaskImageType * mask );
  TubeMaskImageType * GetTubeMaskImage( void );


  /***********/
  /***********/
  /***********/

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
   * Set the border within the image edges that vessels cannot enter */
  void SetBorderInIndexSpace( int border );
  /**
   * Set ExtractBound Minimum */
  void SetExtractBoundMinInIndexSpace( const IndexType & dataMin );

  /**
   * Get ExtractBound Minimum */
  IndexType GetExtractBoundMinInIndexSpace( void ) const;

  /**
   * Set ExtractBound Maximum */
  void SetExtractBoundMaxInIndexSpace( const IndexType & dataMax );

  /**
   * Get ExtractBound Maximum */
  IndexType GetExtractBoundMaxInIndexSpace( void ) const;

  /**
   * Set the radius */
  void SetRadiusInObjectSpace( double radius );

  /**
   * Get the radius */
  double GetRadiusInObjectSpace( void );


  /***********/
  /***********/
  /***********/

  /**
   * Get the ridge extractor */
  typename RidgeExtractorType::Pointer GetRidgeExtractor( void );

  /**
   * Get the radius extractor */
  typename RadiusExtractorType::Pointer GetRadiusExtractor( void );


  /***********/
  /***********/
  /***********/

  /**
   * Return true if a tube is found from the given seed point */
  bool FindLocalTubeInObjectSpace( PointType & x );

  /**
   * Extract the ND tube given the position of the first point
   * and the tube ID */
  TubeType * ExtractTubeInObjectSpace( const PointType & x,
    unsigned int tubeID, bool verbose = false );


  /***********/
  /***********/
  /***********/

  /** Set Seed Mask Image */
  itkSetConstObjectMacro( SeedMask, TubeMaskImageType );
  itkGetConstObjectMacro( SeedMask, TubeMaskImageType );

  /** Set Seed Scale Image */
  itkSetConstObjectMacro( SeedRadiusMask, ImageType );
  itkGetConstObjectMacro( SeedRadiusMask, ImageType );

  /** Set Seed Mask Stride */
  itkSetMacro( SeedMaskStride, int );
  itkGetMacro( SeedMaskStride, int );

  /** Set Seed Index List */
  void SetSeedsInIndexSpaceList( const ContinuousIndexListType & iList );
  void SetSeedsInObjectSpaceList( const PointListType & oList );
  void SetSeedRadiiInObjectSpaceList( const RadiusListType & rList );

  /** Process seed list or seed mask */
  void ProcessSeeds( void );


  /***********/
  /***********/
  /***********/

  /** Get parameter file */
  void LoadParameterFile( const std::string & filename );


  /***********/
  /***********/
  /***********/

  /**
   * Get the list of tubes that have been extracted */
  TubeGroupType * GetTubeGroup( void );

  /**
   * Set the list of tubes that have been extracted */
  void SetTubeGroup( TubeGroupType * tubes );

  /**
   * Smooth a tube */
  void SmoothTube( TubeType * tube, int h=5 );

  /**
   * Add a tube */
  bool AddTube( TubeType * tube );

  /**
   * Delete a tube */
  bool DeleteTube( TubeType * tube );

  /**
   * Set the tube color */
  void SetTubeColor( const vnl_vector< double > & color );

  /**
   * Get the tube color */
  vnl_vector<double> & GetTubeColor( void );


  /***********/
  /***********/
  /***********/


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

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  typename RidgeExtractorType::Pointer   m_RidgeExtractor;
  typename RadiusExtractorType::Pointer  m_RadiusExtractor;

  bool ( *m_IdleCallBack )( void );
  void ( *m_StatusCallBack )( const char *, const char *, int );
  void ( *m_NewTubeCallBack )( TubeType * );
  bool ( *m_AbortProcess )( void );

private:

  TubeExtractor( const Self& );
  void operator=( const Self& );

  typename TubeGroupType::Pointer     m_TubeGroup;

  vnl_vector<double>                  m_TubeColor;

  PointListType                       m_SeedsInObjectSpaceList;
  RadiusListType                      m_SeedRadiiInObjectSpaceList;

  typename TubeMaskImageType::ConstPointer m_SeedMask;
  typename ImageType::ConstPointer         m_SeedRadiusMask;
  int                                      m_SeedMaskStride;



}; // End class TubeExtractor

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeExtractor.hxx"
#endif

#endif // End !defined( __itktubeTubeExtractor_h )
