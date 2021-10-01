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

#ifndef __itktubePointBasedSpatialObjectTransformFilter_h
#define __itktubePointBasedSpatialObjectTransformFilter_h

#include "itktubeSpatialObjectFilter.h"

#include <itkSpatialObject.h>

#include <itkObject.h>
#include <itkObjectFactory.h>

#include <itkPointBasedSpatialObject.h>
#include <itkSpatialObjectPoint.h>

#include <itkContourSpatialObject.h>
#include <itkContourSpatialObjectPoint.h>
#include <itkDTITubeSpatialObject.h>
#include <itkDTITubeSpatialObjectPoint.h>
#include <itkLineSpatialObject.h>
#include <itkLineSpatialObjectPoint.h>
#include <itkSurfaceSpatialObject.h>
#include <itkSurfaceSpatialObjectPoint.h>
#include <itkTubeSpatialObject.h>
#include <itkTubeSpatialObjectPoint.h>

namespace itk
{

namespace tube
{

/**
 *  This class applies a transformation to tubes in a group and returns
 *  a the group with transformed tubes.
 *
 *  \warning Transform Class MUST have a proper implementation of
 *    ::TransformCovariantVector( void )
 *
 *  \warning The scale is applied before computing the transformation.
 *
 *  The resulting tube could be cropped and/or a narrow band could be
 *  defined.
 */
template< class TTransformType, unsigned int TDimension >
class PointBasedSpatialObjectTransformFilter :
  public SpatialObjectFilter< TDimension >
{
public:

  typedef TTransformType                                TransformType;

  typedef SpatialObject< TDimension >                   SpatialObjectType;

  /** Standard class typedefs. */
  typedef PointBasedSpatialObjectTransformFilter< TTransformType, TDimension >
    Self;

  typedef SpatialObjectFilter< TDimension >                  Superclass;

  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  typedef PointBasedSpatialObject< TDimension >              PointBasedType;
  typedef TubeSpatialObject< TDimension >                    TubeType;
  typedef SurfaceSpatialObject< TDimension >                 SurfaceType;
  typedef LineSpatialObject< TDimension >                    LineType;
  typedef DTITubeSpatialObject< TDimension >                 DTITubeType;
  typedef ContourSpatialObject< TDimension >                 ContourType;

  typedef typename SpatialObject< TDimension >::TransformType
    SpatialObjectTransformType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( PointBasedSpatialObjectTransformFilter, SpatialObjectFilter );

  /** Set the Transformation */
  itkSetConstObjectMacro( Transform, TransformType );

  /** Set the Object to Parent transform for the output tubes */
  itkSetConstObjectMacro( OutputObjectToParentTransform, SpatialObjectTransformType );

protected:

  PointBasedSpatialObjectTransformFilter( void );
  virtual ~PointBasedSpatialObjectTransformFilter( void ) {}

  void PrintSelf( std::ostream& os, Indent indent ) const override;

  /** Apply the transformation to the tube */
  void GenerateData( void ) override;


private:

  //purposely not implemented
  PointBasedSpatialObjectTransformFilter( const Self& );
  void operator=( const Self& );

  void UpdateLevel( const SpatialObject< TDimension > * inputSO,
    SpatialObject< TDimension > * parentSO );

  bool Transform( const SpatialObject< TDimension > * inputSO,
    SpatialObject< TDimension > * outputSO );

  void TransformPointBased( const PointBasedSpatialObject< TDimension > * inputSO,
    PointBasedSpatialObject< TDimension > * outputSO );

  void TransformTube( const TubeSpatialObject< TDimension > * inputSO,
    TubeSpatialObject< TDimension > * outputSO );

  void TransformSurface( const SurfaceSpatialObject< TDimension > * inputSO,
    SurfaceSpatialObject< TDimension > * outputSO );

  typename TransformType::ConstPointer              m_Transform;

  typename SpatialObjectTransformType::ConstPointer m_OutputObjectToParentTransform;

}; // End class PointBasedSpatialObjectTransformFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePointBasedSpatialObjectTransformFilter.hxx"
#endif

#endif // End !defined( __itktubePointBasedSpatialObjectTransformFilter_h )
