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

#ifndef __itktubeSubSampleTubeTreeSpatialObjectFilter_h
#define __itktubeSubSampleTubeTreeSpatialObjectFilter_h

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

namespace itk
{

namespace tube
{

/** \class SubSampleTubeTreeSpatialObjectFilter
 * \brief Sub-sample tubes within a SpatialObject hierarchy.
 *
 * The input to this SpatialObjectFilter can be a single TubeSpatialObject
 * or a hierarchy of SpatialObject's that contain tubes to be sub-sampled.
 * All TubeSpatialObjects in the output hierarchy will be sub-sampled by
 * the \c Sampling factor.  Non-Tube spatial objects are passed to the
 * output unchanged.
 *
 * \sa SubSampleTubeSpatialObjectFilter
 */
template< class TSpatialObject, class TTubeSpatialObject >
class SubSampleTubeTreeSpatialObjectFilter
  : public SpatialObjectToSpatialObjectFilter< TSpatialObject,
    TSpatialObject >
{
public:
  /** Standard class typedefs. */
  typedef SubSampleTubeTreeSpatialObjectFilter      Self;
  typedef SpatialObjectToSpatialObjectFilter< TSpatialObject,
    TSpatialObject >                                Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  typedef TSpatialObject     SpatialObjectType;
  typedef TTubeSpatialObject TubeSpatialObjectType;

  /** Run-time type information ( and related methods ).   */
  itkTypeMacro( SubSampleTubeTreeSpatialObjectFilter,
    SpatialObjectToSpatialObjectFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  itkStaticConstMacro( ObjectDimension, unsigned int,
    SpatialObjectType::ObjectDimension );

  /** Set the sampling factor.  The output points taken every sampling
  * factor from the input points. */
  itkSetClampMacro( Sampling, SizeValueType, 1, NumericTraits<
    SizeValueType >::max() );
  itkGetConstMacro( Sampling, SizeValueType );

protected:
  typedef SpatialObject< ObjectDimension > SpatialObjectBaseType;

  SubSampleTubeTreeSpatialObjectFilter( void );
  virtual ~SubSampleTubeTreeSpatialObjectFilter( void );

  virtual void GenerateData( void );

  /** Sub-sample at the tubes at a given level, then sub-sample their
  * children. */
  virtual void SubSampleLevel( const SpatialObjectBaseType * input,
    SpatialObjectBaseType * output, bool graftOutput=false );

private:
  // purposely not implemented
  SubSampleTubeTreeSpatialObjectFilter( const Self & );

  // purposely not implemented
  void operator=( const Self & );

  SizeValueType m_Sampling;

}; // End class SubSampleTubeTreeSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSubSampleTubeTreeSpatialObjectFilter.hxx"
#endif

#endif // End !defined( __itktubeSubSampleTubeTreeSpatialObjectFilter_h )
