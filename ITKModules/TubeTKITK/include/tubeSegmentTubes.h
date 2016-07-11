/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeSegmentTubes_h
#define __tubeSegmentTubes_h

#include "itktubeTubeExtractor.h"
#include "itktubeRidgeExtractor.h"
#include "itktubeRadiusExtractor2.h"
#include "itkObject.h"

#include "tubeWrappingMacros.h"

namespace tube
{
/** \class SegmentTubes
 *
 *  \ingroup TubeTKITK
 */

template< class TInputImage >
class SegmentTubes:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef SegmentTubes                                 Self;
  typedef itk::SmartPointer< Self >                    Pointer;
  typedef itk::SmartPointer< const Self >              ConstPointer;
  typedef TInputImage                                  ImageType;
  typedef typename ImageType::PixelType                PixelType;
  typedef typename ImageType::IndexType                IndexType;
  typedef itk::tube::TubeExtractor< ImageType >        TubeExtractorFilterType;
  typedef itk::tube::RidgeExtractor< ImageType >       RidgeExtractorFilterType;
  typedef itk::tube::RadiusExtractor2< ImageType >     RidgeExtractor2FilterType;
  typedef typename TubeExtractorFilterType::TubeMaskImageType
    TubeMaskImageType;
  typedef typename TubeExtractorFilterType::TubeType   TubeType;
  typedef typename TubeType::Pointer                   TubePointerType;

  typedef typename TubeExtractorFilterType::ContinuousIndexType
    ContinuousIndexType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SegmentTubes, Object );

  /** Set the source image. */
  tubeWrapSetObjectMacro( InputImage, ImageType, TubeExtractorFilter );

  /** Set the radius image. */
  tubeWrapSetObjectMacro( RadiusInputImage, ImageType, TubeExtractorFilter );

  /** Set/Get radius value */
  tubeWrapSetMacro( Radius, double, TubeExtractorFilter );
  tubeWrapGetMacro( Radius, double, TubeExtractorFilter );

  /** Set the tube mask image. */
  void SetTubeMaskImage( typename TubeMaskImageType::Pointer & mask );
  typename itk::tube::TubeExtractor< TInputImage >::TubeMaskImageType::Pointer
    GetTubeMaskImage( void );

  /** Add a tube */
  bool AddTube( TubeType * tube );

  /** Set/Get debug status */
  tubeWrapSetMacro( Debug, bool, TubeExtractorFilter );
  tubeWrapGetMacro( Debug, bool, TubeExtractorFilter );

  /** Set ExtractBound Minimum */
  void SetExtractBoundMin( const IndexType & dataMin );

  /** Set ExtractBound Maximum */
  void SetExtractBoundMax( const IndexType & dataMax );

  /*** Extract the ND tube given the position of the first point
   * and the tube ID */
  typename itk::tube::TubeExtractor< TInputImage >::TubeType::Pointer
    ExtractTube( const ContinuousIndexType & x,
    unsigned int tubeID,
    bool verbose = false );

  /** Get the list of tubes that have been extracted */
  typename itk::tube::TubeExtractor< TInputImage >::TubeGroupType::Pointer
    GetTubeGroup( void );

  /* Get the ridge extractor */
  typename itk::tube::RidgeExtractor< ImageType >::Pointer GetRidgeOp( void );

    /*Get the radius extractor */
  typename itk::tube::RadiusExtractor2< ImageType >::Pointer GetRadiusOp( void );

    /*Get the tube extractor */
  typename itk::tube::TubeExtractor< ImageType >::Pointer GetTubeOp( void );

protected:
  SegmentTubes( void );
  ~SegmentTubes() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeTubeExtractor parameters **/
  SegmentTubes( const Self & );
  void operator=( const Self & );

  typename TubeExtractorFilterType::Pointer m_TubeExtractorFilter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentTubes.hxx"
#endif

#endif // End !defined( __tubeSegmentTubes_h )
