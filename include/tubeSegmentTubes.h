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
#ifndef __tubeSegmentTubes_h
#define __tubeSegmentTubes_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeTubeExtractor.h"
#include "itktubeTubeExtractorIO.h"

namespace tube
{
/** \class SegmentTubes
 *
 *  \ingroup TubeTK
 */

template< class TInputImage >
class SegmentTubes:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SegmentTubes                                 Self;
  typedef itk::ProcessObject                           Superclass;
  typedef itk::SmartPointer< Self >                    Pointer;
  typedef itk::SmartPointer< const Self >              ConstPointer;

  typedef TInputImage                                  ImageType;
  typedef typename ImageType::PixelType                PixelType;
  typedef typename ImageType::IndexType                IndexType;

  typedef itk::tube::TubeExtractor< ImageType >        FilterType;
  typedef itk::tube::RidgeExtractor< ImageType >       RidgeFilterType;
  typedef itk::tube::RadiusExtractor3< ImageType >     RadiusFilterType;

  typedef typename FilterType::TubeMaskImageType       TubeMaskImageType;
  typedef typename FilterType::TubeType                TubeType;
  typedef typename FilterType::TubeGroupType           TubeGroupType;

  typedef typename FilterType::ContinuousIndexType     ContinuousIndexType;
  typedef std::vector< ContinuousIndexType >           ContinuousIndexListType;
  typedef typename FilterType::RadiusListType          RadiusListType;
  typedef typename FilterType::PointType               PointType;
  typedef std::vector< PointType >                     PointListType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentTubes, ProcessObject );

  /***/
  /***/
  /***/

  /** Set the source image. */
  void SetInput( ImageType * img )
  { m_Filter->SetInputImage( img ); };

  const ImageType * GetInput( void )
  { return m_Filter->GetInputImage(); };

  /** Set the radius image. */
  tubeWrapSetObjectMacro( RadiusInputImage, ImageType, Filter );
  tubeWrapGetConstObjectMacro( RadiusInputImage, ImageType, Filter );

  /** The image that records extracted vessels */
  tubeWrapSetObjectMacro( TubeMaskImage, TubeMaskImageType, Filter );
  tubeWrapGetObjectMacro( TubeMaskImage, TubeMaskImageType, Filter );

  /***/
  /***/
  /***/

  itkSetMacro( Verbose, bool );
  itkGetMacro( Verbose, bool );

  tubeWrapSetMacro( DataMin, double, Filter );
  tubeWrapGetMacro( DataMin, double, Filter );

  tubeWrapSetMacro( DataMax, double, Filter );
  tubeWrapGetMacro( DataMax, double, Filter );

  void SetDataMinMaxLimits( double minLimit, double maxLimit )
  { return this->m_Filter->SetDataMinMaxLimits( minLimit, maxLimit ); };

  tubeWrapForceSetMacro( BorderInIndexSpace, int, Filter );

  tubeWrapSetConstReferenceMacro( ExtractBoundMinInIndexSpace, IndexType,
    Filter );
  tubeWrapGetMacro( ExtractBoundMinInIndexSpace, IndexType, Filter );

  tubeWrapSetConstReferenceMacro( ExtractBoundMaxInIndexSpace, IndexType,
    Filter );
  tubeWrapGetMacro( ExtractBoundMaxInIndexSpace, IndexType, Filter );

  tubeWrapSetMacro( RadiusInObjectSpace, double, Filter );
  tubeWrapGetMacro( RadiusInObjectSpace, double, Filter );

  tubeWrapSetMacro( OptimizeRadius, bool, Filter );
  tubeWrapGetMacro( OptimizeRadius, bool, Filter );

  /***/
  /***/
  /***/

  bool FindLocalTubeInObjectSpace( PointType & x )
  { return m_Filter->FindLocalTubeInObjectSpace( x ); }

  TubeType * ExtractTubeInObjectSpace( const PointType & x,
    unsigned int tubeId )
  { return m_Filter->ExtractTubeInObjectSpace( x, tubeId, m_Verbose ); }

  /***/
  /***/
  /***/

  tubeWrapForceSetMacro( SeedsInIndexSpaceList, ContinuousIndexListType, Filter );
  tubeWrapForceSetMacro( SeedsInObjectSpaceList, PointListType, Filter );
  tubeWrapForceSetMacro( SeedRadiiInObjectSpaceList, RadiusListType, Filter );

  tubeWrapSetObjectMacro( SeedMask, TubeMaskImageType, Filter );
  tubeWrapGetConstObjectMacro( SeedMask, TubeMaskImageType, Filter );

  tubeWrapSetObjectMacro( SeedRadiusMask, ImageType, Filter );
  tubeWrapGetConstObjectMacro( SeedRadiusMask, ImageType, Filter );

  tubeWrapSetMacro( UseSeedMaskAsProbabilities, bool, Filter );
  tubeWrapGetMacro( UseSeedMaskAsProbabilities, bool, Filter );

  tubeWrapSetMacro( SeedExtractionMinimumSuccessRatio, double, Filter );
  tubeWrapGetMacro( SeedExtractionMinimumSuccessRatio, double, Filter );

  tubeWrapSetMacro( SeedExtractionMinimumProbability, double, Filter );
  tubeWrapGetMacro( SeedExtractionMinimumProbability, double, Filter );

  tubeWrapSetMacro( SeedMaskMaximumNumberOfPoints, unsigned int, Filter );
  tubeWrapGetMacro( SeedMaskMaximumNumberOfPoints, unsigned int, Filter );

  tubeWrapSetMacro( SeedMaskStride, int, Filter );
  tubeWrapGetMacro( SeedMaskStride, int, Filter );

  void ProcessSeeds( void )
  { this->m_Filter->ProcessSeeds( m_Verbose ); };

  /** Load parameters of tube extraction from a file */
  void LoadParameterFile( const std::string & filename )
  { ::itk::tube::TubeExtractorIO< ImageType > teReader;
    teReader.SetTubeExtractor( this->m_Filter );
    teReader.Read( filename.c_str() ); };

  /** Get the list of tubes that have been extracted */
  tubeWrapSetObjectMacro( TubeGroup, TubeGroupType, Filter );
  tubeWrapGetObjectMacro( TubeGroup, TubeGroupType, Filter );

  /** Add a tube */
  bool AddTube( TubeType * tube )
  { return this->m_Filter->AddTube( tube ); };

  /** Delete a tube */
  bool DeleteTube( TubeType * tube )
  { return this->m_Filter->DeleteTube( tube ); };

  /** Set/Get debug status */
  tubeWrapSetMacro( Debug, bool, Filter );
  tubeWrapGetMacro( Debug, bool, Filter );

  /***/
  /***/
  /***/

  tubeWrapSetMacro( StepX, double, RidgeFilter );
  tubeWrapGetMacro( StepX, double, RidgeFilter );

  tubeWrapSetMacro( MaxTangentChange, double, RidgeFilter );
  tubeWrapGetMacro( MaxTangentChange, double, RidgeFilter );

  tubeWrapSetMacro( MaxXChange, double, RidgeFilter );
  tubeWrapGetMacro( MaxXChange, double, RidgeFilter );

  tubeWrapSetMacro( MinRidgeness, double, RidgeFilter );
  tubeWrapGetMacro( MinRidgeness, double, RidgeFilter );
  tubeWrapSetMacro( MinRidgenessStart, double, RidgeFilter );
  tubeWrapGetMacro( MinRidgenessStart, double, RidgeFilter );

  tubeWrapSetMacro( MinRoundness, double, RidgeFilter );
  tubeWrapGetMacro( MinRoundness, double, RidgeFilter );
  tubeWrapSetMacro( MinRoundnessStart, double, RidgeFilter );
  tubeWrapGetMacro( MinRoundnessStart, double, RidgeFilter );

  tubeWrapSetMacro( MinCurvature, double, RidgeFilter );
  tubeWrapGetMacro( MinCurvature, double, RidgeFilter );
  tubeWrapSetMacro( MinCurvatureStart, double, RidgeFilter );
  tubeWrapGetMacro( MinCurvatureStart, double, RidgeFilter );

  tubeWrapSetMacro( MinLevelness, double, RidgeFilter );
  tubeWrapGetMacro( MinLevelness, double, RidgeFilter );
  tubeWrapSetMacro( MinLevelnessStart, double, RidgeFilter );
  tubeWrapGetMacro( MinLevelnessStart, double, RidgeFilter );

  tubeWrapSetMacro( Scale, double, RidgeFilter );
  tubeWrapGetMacro( Scale, double, RidgeFilter );

  tubeWrapSetMacro( ScaleKernelExtent, double, RidgeFilter );
  tubeWrapGetMacro( ScaleKernelExtent, double, RidgeFilter );

  tubeWrapSetMacro( DynamicScale, bool, RidgeFilter );
  tubeWrapGetMacro( DynamicScale, bool, RidgeFilter );

  double RidgenessInObjectSpace( const PointType & x )
  { m_Ridgeness = m_RidgeFilter->Ridgeness( x, m_Intensity, m_Roundness,
      m_Curvature, m_Levelness );  return m_Ridgeness; };

  tubeWrapGetMacro( CurrentRidgeness, double, RidgeFilter );
  tubeWrapGetMacro( CurrentRoundness, double, RidgeFilter );
  tubeWrapGetMacro( CurrentCurvature, double, RidgeFilter );
  tubeWrapGetMacro( CurrentLevelness, double, RidgeFilter );

  PointType FindLocalRidgeInObjectSpace( const PointType & x )
  { PointType localX = x; m_RidgeFilter->LocalRidge( localX, m_Verbose );
    return localX; };

  TubeType * ExtractRidgeInObjectSpace( const PointType & x, int tubeId )
  { return m_RidgeFilter->ExtractRidge( x, tubeId, m_Verbose ); };


  /***/
  /***/
  /***/

  tubeWrapSetMacro( RadiusMin, double, RadiusFilter );
  tubeWrapGetMacro( RadiusMin, double, RadiusFilter );

  tubeWrapSetMacro( RadiusMax, double, RadiusFilter );
  tubeWrapGetMacro( RadiusMax, double, RadiusFilter );

  tubeWrapSetMacro( RadiusStart, double, RadiusFilter );
  tubeWrapGetMacro( RadiusStart, double, RadiusFilter );

  tubeWrapSetMacro( MinMedialness, double, RadiusFilter );
  tubeWrapGetMacro( MinMedialness, double, RadiusFilter );
  tubeWrapSetMacro( MinMedialnessStart, double, RadiusFilter );
  tubeWrapGetMacro( MinMedialnessStart, double, RadiusFilter );

  tubeWrapSetMacro( KernelNumberOfPoints, unsigned int, RadiusFilter );
  tubeWrapGetMacro( KernelNumberOfPoints, unsigned int, RadiusFilter );

  tubeWrapSetMacro( KernelPointStep, unsigned int, RadiusFilter );
  tubeWrapGetMacro( KernelPointStep, unsigned int, RadiusFilter );

  bool ExtractRadii( TubeType * tube )
  { return m_RadiusFilter->ExtractRadii( tube ); };

protected:
  SegmentTubes( void );
  ~SegmentTubes() {};
  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itktubeTubeExtractor parameters **/
  SegmentTubes( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) 
    override {};

  typename FilterType::Pointer           m_Filter;
  typename RidgeFilterType::Pointer      m_RidgeFilter;
  typename RadiusFilterType::Pointer     m_RadiusFilter;

  bool m_Verbose;

  double m_Ridgeness;
  double m_Intensity;
  double m_Roundness;
  double m_Curvature;
  double m_Levelness;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentTubes.hxx"
#endif

#endif // End !defined( __tubeSegmentTubes_h )
