/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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
  tubeWrapSetObjectMacro( InputImage, ImageType, Filter );
  tubeWrapGetConstObjectMacro( InputImage, ImageType, Filter );

  /** Set the radius image. */
  tubeWrapSetObjectMacro( RadiusInputImage, ImageType, Filter );
  tubeWrapGetConstObjectMacro( RadiusInputImage, ImageType, Filter );

  /** The image that records extracted vessels */
  tubeWrapSetObjectMacro( TubeMaskImage, TubeMaskImageType, Filter );
  tubeWrapGetObjectMacro( TubeMaskImage, TubeMaskImageType, Filter );

  /***/
  /***/
  /***/

  tubeWrapSetMacro( DataMin, double, Filter );
  tubeWrapGetMacro( DataMin, double, Filter );

  tubeWrapSetMacro( DataMax, double, Filter );
  tubeWrapGetMacro( DataMax, double, Filter );

  tubeWrapForceSetMacro( BorderInIndexSpace, int, Filter );

  tubeWrapSetConstReferenceMacro( ExtractBoundMinInIndexSpace, IndexType,
    Filter );
  tubeWrapGetMacro( ExtractBoundMinInIndexSpace, IndexType, Filter );

  tubeWrapSetConstReferenceMacro( ExtractBoundMaxInIndexSpace, IndexType,
    Filter );
  tubeWrapGetMacro( ExtractBoundMaxInIndexSpace, IndexType, Filter );

  tubeWrapSetMacro( RadiusInObjectSpace, double, Filter );
  tubeWrapGetMacro( RadiusInObjectSpace, double, Filter );

  /***/
  /***/
  /***/

  bool FindLocalTubeInObjectSpace( PointType & x )
  { return m_Filter->FindLocalTubeInObjectSpace( x ); }

  TubeType * ExtractTubeInObjectSpace( const PointType & x,
    unsigned int tubeId, bool verbose = false )
  { return m_Filter->ExtractTubeInObjectSpace( x, tubeId, verbose ); }

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
  tubeWrapSetMacro( SeedMaskStride, int, Filter );
  tubeWrapGetMacro( SeedMaskStride, int, Filter );

  void ProcessSeeds( void )
  { this->m_Filter->ProcessSeeds(); }

  /** Load parameters of tube extraction from a file */
  void LoadParameterFile( const std::string & filename )
  { this->m_Filter->LoadParameterFile( filename ); }

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

protected:
  SegmentTubes( void );
  ~SegmentTubes() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeTubeExtractor parameters **/
  SegmentTubes( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename FilterType::Pointer m_Filter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentTubes.hxx"
#endif

#endif // End !defined( __tubeSegmentTubes_h )
