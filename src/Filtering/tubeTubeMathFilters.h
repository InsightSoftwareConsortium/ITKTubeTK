/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __tubeTubeMathFilters_h
#define __tubeTubeMathFilters_h

#include <itkImageFileReader.h>
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"

namespace tube
{
template< unsigned int DimensionT, class ImagePixelT=float >
class TubeMathFilters
{
public:
  //typedefs
  typedef itk::GroupSpatialObject< DimensionT >         TubeGroupType;
  typedef typename TubeGroupType::ChildrenListPointer   TubeListPointerType;

  typedef itk::TubeSpatialObject< DimensionT >          TubeType;
  typedef typename TubeType::TubePointType              TubePointType;

  typedef itk::Image< ImagePixelT, DimensionT >         ImageType;
  typedef itk::Image< float, DimensionT >               FloatImageType;
  typedef typename FloatImageType::OffsetType           VectorPixelType;
  typedef itk::Image< VectorPixelType, DimensionT >     VectorImageType;

  typedef typename TubeType::PointType                  PositionType;
  typedef itk::IndexValueType                           TubeIdType;
  typedef typename TubeType::TubePointListType          TubePointListType;

  TubeMathFilters();
  ~TubeMathFilters();

  void SetInputTubeGroup( TubeGroupType * inputTubeGroup );
  void SetInputTube( TubeType * inputTube );
  
  typename TubeGroupType::Pointer & GetOutputTubeGroup( void );
  typename TubeType::Pointer &      GetOutputTube( void );

  /** SetCurrentTubeId
   *  Specifies which tube to process.  Use -1 to specify all tubes. */
  void SetCurrentTubeId( int currentTubeId );
  void SetUseAllTubes( void );

  void SetPointValues( std::string propertyId, double val, double blend=1 );

  void SetPointValuesFromImage( const ImageType * inputImage,
    std::string propertyId, double blend=1 );

  void SetPointValuesFromImageMean( const ImageType * inputImage,
    std::string propertyId );

  void ComputeTubeRegions( const ImageType * referenceImage );

  void SetPointValuesFromTubeRegions(
    const ImageType * inputImage,
    const std::string & propertyId,
    double minRFactor=1, double maxRFactor=3 );

  void SetPointValuesFromTubeRadius(
    const ImageType * inputImage,
    const std::string & propertyId,
    double minRFactor=1, double maxRFactor=3 );

  /** Run Fill Gap on the tube-tree. */
  void FillGapToParent( double stepSize = 0.1 );

  /** Smooth a tube
   * The parameter h has different meanings when using different smoothing
   * functions:
   *
   * smoothFunction = SMOOTH_TUBE_USING_INDEX_AVERAGE:
   *    h is half of the window size
   * smoothFunction = SMOOTH_TUBE_USING_INDEX_GAUSSIAN:
   *    h is the gaussian's standard deviation
   */
  enum SmoothTubeFunctionEnum { SMOOTH_TUBE_USING_INDEX_AVERAGE,
    SMOOTH_TUBE_USING_INDEX_GAUSSIAN };
  void SmoothTube( double h = 2,
    SmoothTubeFunctionEnum smoothFunction = SMOOTH_TUBE_USING_INDEX_AVERAGE );
  void SmoothTubeProperty( const std::string & propertyId, double h = 2,
    SmoothTubeFunctionEnum smoothFunction = SMOOTH_TUBE_USING_INDEX_AVERAGE );

  void RenumberTubes( void );
  void RenumberPoints( void );

  void SubsampleTube( int N = 2 );
  
  double ComputeTubeLength( void );

protected:

  static void InterpolatePath(
    typename TubeType::TubePointType * parentNearestPoint,
    typename TubeType::TubePointType * childEndPoint,
    float stepSize,
    typename TubeType::TubePointListType & newTubePoints );

private:
  typename TubeType::Pointer      m_InputTube;
  typename TubeGroupType::Pointer m_InputTubeGroup;

  int m_CurrentTubeId;

  typename FloatImageType::Pointer         m_TubePointIdImage;
  typename FloatImageType::Pointer         m_TubeRadiusImage;
  typename FloatImageType::Pointer         m_TubeDistanceImage;
  typename VectorImageType::Pointer        m_TubeDirectionImage;


}; // End class ImageFilters

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeTubeMathFilters.hxx"
#endif

#endif // End !defined( __tubeTubeMathFilters_h )
