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

#ifndef __itktubeSegmentTubes_h
#define __itktubeSegmentTubes_h

#include "itkImage.h"
#include "itkObject.h"
#include "itktubeTubeExtractor.h"
#include "itktubeRidgeExtractor.h"
#include "itktubeTubeExtractorIO.h"

namespace itk
{

namespace tube
{

/**
 * This class segments a tube
 *
 * \sa SegmentTubes
 */

template< class TInputImage >
class SegmentTubes : public Object
{
public:

  /**
   * Standard self typedef */
  typedef SegmentTubes                  Self;
  typedef Object                        Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;

  typedef TInputImage                   ImageType;
  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::IndexType IndexType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentTubes, Object );

  /** Standard for the number of dimension */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  typedef itk::tube::TubeExtractor< ImageType >        TubeExtractorFilterType;
  typedef float                                        ScaleType;
  typedef typename TubeExtractorFilterType::ContinuousIndexType
    ContinuousIndexType;
  typedef typename ImageType::PointType                PointType;
  typedef typename TubeExtractorFilterType::TubeMaskImageType
  TubeMaskImageType;
  typedef itk::Image< ScaleType, ImageDimension >      ScaleImageType;

  /**  Type definition for VesselTubeSpatialObject */
  typedef itk::GroupSpatialObject< ImageDimension >    TubeGroupType;
  typedef typename TubeExtractorFilterType::TubeType   TubeType;
  typedef typename TubeType::TransformType             TubeTransformType;
  /* Parameters file type*/
  typedef itk::tube::TubeExtractorIO< ImageType >      TubeExtractorIOType;
  typedef itk::tube::RidgeExtractor< ImageType >       RidgeExtractorFilterType;

  /** Set/Get the input image */
  itkSetObjectMacro( InputImage, ImageType );
  itkGetObjectMacro( InputImage, ImageType );
  itkSetObjectMacro( RadiusInputImage, ImageType );
  itkGetObjectMacro( RadiusInputImage, ImageType );

  /** Set/Get scale */
  itkSetMacro( Scale, ScaleType );
  itkGetMacro( Scale, ScaleType );

  /** Set Seed Index List */
  void SetSeedIndexList( std::vector< ContinuousIndexType > );
  void SetSeedPhysicalCoordinatesList( std::vector< PointType > );
  void SetSeedIndexFromFileList
    ( std::vector< ContinuousIndexType >, std::vector< ScaleType > );

  /** Set Seed Mask Image */
  itkSetObjectMacro( SeedMask, TubeMaskImageType );

  /** Set Seed Scale Image */
  itkSetObjectMacro( ScaleMask, ScaleImageType );

  /** Set Seed Mask Stride */
  itkSetMacro( SeedMaskStride, int );

  /** Set Seed Mask Image */
  itkSetObjectMacro( ExistingTubesMask, TubeMaskImageType );

  /** Set Existing tubes */
  void SetTubeGroup( TubeGroupType * tubes );

  /** Set Parameter File */
  itkSetMacro( ParameterFile, std::string );

  /** Set Border */
  itkSetMacro( Border, double );

  /* Get the list of tubes that have been extracted */
  typename TubeGroupType::Pointer GetTubeGroup( void );

  /** Get the tube mask image */
  typename TubeMaskImageType::Pointer GetTubeMaskImage( void );

  void Update();
protected:

  SegmentTubes( void );
  virtual ~SegmentTubes( void );
  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  SegmentTubes( const Self& );
  void operator=( const Self& );

  typename ImageType::Pointer               m_InputImage;
  typename ImageType::Pointer               m_RadiusInputImage;
  typename TubeExtractorFilterType::Pointer m_TubeExtractorFilter;
  typename TubeMaskImageType::Pointer       m_SeedMask;
  typename ScaleImageType::Pointer          m_ScaleMask;
  typename TubeMaskImageType::Pointer       m_ExistingTubesMask;

  typename TubeGroupType::Pointer           m_ExistingTubes;

  ScaleType                          m_Scale;
  std::vector< ContinuousIndexType > m_SeedIndexList;
  std::vector< ScaleType >           m_SeedRadiusList;
  std::vector< PointType >           m_SeedPhysicalCoordinatesList;
  std::vector< ContinuousIndexType > m_SeedIndexFromFileList;
  std::vector< ScaleType >           m_SeedScaleFromFileList;
  int                                m_SeedMaskStride;
  bool                               m_UseExistingTubes;
  std::string                        m_ParameterFile;
  double                             m_Border;
  typename TubeGroupType::Pointer    m_TubeGroup;


}; // End class SegmentTubes

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSegmentTubes.hxx"
#endif

#endif // End !defined( __itktubeSegmentTubes_h )
