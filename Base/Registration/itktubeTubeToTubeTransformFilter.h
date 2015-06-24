/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeTubeToTubeTransformFilter_h
#define __itktubeTubeToTubeTransformFilter_h

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

#include <itkGroupSpatialObject.h>
#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkVesselTubeSpatialObject.h>
#include <itkVesselTubeSpatialObjectPoint.h>

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
template< class TTransformType, unsigned int VDimension >
class TubeToTubeTransformFilter :
  public SpatialObjectToSpatialObjectFilter<
  GroupSpatialObject< VDimension >,
  GroupSpatialObject< VDimension > >
{
public:

  typedef GroupSpatialObject<VDimension>                     GroupType;

  /** Standard class typedefs. */
  typedef TubeToTubeTransformFilter<TTransformType,VDimension>
    Self;

  typedef SpatialObjectToSpatialObjectFilter<GroupType, GroupType>
    Superclass;

  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  typedef VesselTubeSpatialObject<VDimension>                TubeType;

  typedef typename VesselTubeSpatialObject<VDimension>::TransformType
    TubeTransformType;

  typedef typename VesselTubeSpatialObject<VDimension>::
    AffineGeometryFrameType                                  TubeFrameType;

  /** Typedef for the transformations */
  typedef TTransformType                                     TransformType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TubeToTubeTransformFilter,
    SpatialObjectToSpatialObjectFilter );

  /** Apply the transformation to the tube */
  void Update( void );

  /** Set the Transformation */
  itkSetObjectMacro( Transform, TransformType );

  /** Set the Transformation */
  itkSetObjectMacro( OutputIndexToObjectFrame, TubeFrameType );

  /** Get the output tubenet */
  itkGetObjectMacro( Output, GroupType );

protected:

  TubeToTubeTransformFilter( void );
  virtual ~TubeToTubeTransformFilter( void ) {}
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  TubeToTubeTransformFilter(const Self&); //purposely not implemented
  void operator=(const Self&);            //purposely not implemented

  typename TransformType::Pointer            m_Transform;

  typename TubeType::AffineGeometryFrameType::Pointer
    m_OutputIndexToObjectFrame;

  typename GroupType::Pointer                m_Output;

}; // End class TubeToTubeTransformFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeToTubeTransformFilter.hxx"
#endif

#endif // End !defined(__itktubeTubeToTubeTransformFilter_h)
