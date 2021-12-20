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
#ifndef __tubeTubeMath_h
#define __tubeTubeMath_h

// ITK includes
#include "itkProcessObject.h"
#include "itkCastImageFilter.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "tubeTubeMathFilters.h"

namespace tube
{
/** \class TubeMath
 *  \brief Common image processing functions.
 *
 *  \ingroup TubeTK
 */

template< unsigned int DimensionT, class ImagePixelT=float >
class TubeMath:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef TubeMath                                   Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef tube::TubeMathFilters< DimensionT, ImagePixelT >  FilterType;

  typedef typename FilterType::TubeGroupType         TubeGroupType;
  typedef typename FilterType::TubeType              TubeType;
  typedef typename FilterType::ImageType             ImageType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubeMath, ProcessObject );

  void SetInputTubeGroup( TubeGroupType * tubeGroup )
  { m_Filter.SetInputTubeGroup( tubeGroup ); this->Modified(); };

  void SetInputTube( TubeType * tube )
  { m_Filter.SetInputTube(tube); this->Modified(); };

  TubeGroupType * GetOutputTubeGroup( void )
  { return m_Filter.GetOutputTubeGroup(); };

  TubeType * GetOutputTube()
  { return m_Filter.GetOutputTube(); };

  void SetCurrentTubeId( int tubeId )
  { m_Filter.SetCurrentTubeId(tubeId); this->Modified(); };

  void SetUseAllTubes( void )
  { m_Filter.SetUseAllTubes(); this->Modified(); };

  void SetPointValues( std::string propertyId, double val, double blend=1 )
  { m_Filter.SetPointValues(propertyId, val, blend); this->Modified(); };

  void SetPointValuesFromImage( const ImageType * inputImage,
    std::string propertyId )
  { m_Filter.SetPointValuesFromImage(inputImage, propertyId);
    this->Modified(); };

  void SetPointValuesFromImageMean( const ImageType * inputImage,
    std::string propertyId )
  { m_Filter.SetPointValuesFromImageMean(inputImage, propertyId);
    this->Modified(); };

  void ComputeTubeRegions( const ImageType * referenceImage )
  { m_Filter.ComputeTubeRegions(referenceImage); this->Modified(); };

  void SetPointValuesFromTubeRegions( const ImageType * inputImage,
    const std::string & propertyId, double minRFactor=1, double maxRFactor=3 )
  { m_Filter.SetPointValuesFromTubeRegions(inputImage, propertyId,
    minRFactor, maxRFactor); this->Modified(); };

  void SetPointValuesFromTubeRadius( const ImageType * inputImage,
    const std::string & propertyId, double minRFactor=1, double maxRFactor=3 )
  { m_Filter.SetPointValuesFromTubeRadius(inputImage, propertyId,
    minRFactor, maxRFactor); this->Modified(); };

  void SmoothTube( double h=2, const std::string & smoothTubeFunction =
    "SMOOTH_TUBE_USING_INDEX_AVERAGE" )
  { if( smoothTubeFunction == "SMOOTH_TUBE_USING_INDEX_AVERAGE" )
      { m_Filter.SmoothTube(h, FilterType::SMOOTH_TUBE_USING_INDEX_AVERAGE ); }
    else
      { m_Filter.SmoothTube(h, FilterType::SMOOTH_TUBE_USING_INDEX_GAUSSIAN ); }
    this->Modified(); };

  void SmoothTubeProperty( const std::string & propertyId, 
    double h=2, const std::string & smoothTubeFunction =
    "SMOOTH_TUBE_USING_INDEX_AVERAGE" )
  { if( smoothTubeFunction == "SMOOTH_TUBE_USING_INDEX_AVERAGE" )
      { m_Filter.SmoothTubeProperty(propertyId, h,
        FilterType::SMOOTH_TUBE_USING_INDEX_AVERAGE ); }
    else
      { m_Filter.SmoothTubeProperty(propertyId, h,
        FilterType::SMOOTH_TUBE_USING_INDEX_GAUSSIAN ); }
    this->Modified(); };

  void RenumberPoints( void )
  { m_Filter.RenumberPoints(); this->Modified(); };

  void SubsampleTube( int N=2 )
  { m_Filter.SubsampleTube(N); this->Modified(); };

  double ComputeTubeLength( void )
  { return m_Filter.ComputeTubeLength(); this->Modified(); };

protected:
  TubeMath( void );
  ~TubeMath() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** TubeMath parameters **/
  TubeMath( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  FilterType m_Filter;
};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeTubeMath.hxx"
#endif

#endif // End !defined( __tubeTubeMath_h )
