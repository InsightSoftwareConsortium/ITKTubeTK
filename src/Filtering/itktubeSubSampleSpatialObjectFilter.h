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

#ifndef __itktubeSubSampleSpatialObjectFilter_h
#define __itktubeSubSampleSpatialObjectFilter_h

#include "itktubeSpatialObjectFilter.h"
#include "itkTubeSpatialObject.h"

namespace itk
{

namespace tube
{

/** \class SubSampleSpatialObjectFilter
 * \brief Sub-sample objects within a SpatialObject hierarchy.
 *
 * The input to this SpatialObjectFilter can be a single SpatialObject
 * or a hierarchy of SpatialObject's that contain objects to be sub-sampled.
 * All supported SpatialObjects in the output hierarchy will be sub-sampled by
 * the \c Sampling factor.  Non-supported spatial objects are passed to the
 * output unchanged.
 */
template< unsigned int ObjectDimension=3 >
class SubSampleSpatialObjectFilter
  : public SpatialObjectFilter< ObjectDimension >
{
public:
  /** Standard class typedefs. */
  typedef SubSampleSpatialObjectFilter              Self;
  typedef SpatialObjectFilter< ObjectDimension >    Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  typedef SpatialObject<ObjectDimension>            SpatialObjectType;

  /** Run-time type information ( and related methods ).   */
  itkTypeMacro( SubSampleSpatialObjectFilter,
    SpatialObjectFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Set the sampling factor.  The output points taken every sampling
  * factor from the input points. */
  itkSetClampMacro( Sampling, SizeValueType, 1, NumericTraits<
    SizeValueType >::max() );
  itkGetConstMacro( Sampling, SizeValueType );

protected:
  typedef TubeSpatialObject< ObjectDimension > TubeSpatialObjectType;

  SubSampleSpatialObjectFilter( void );
  virtual ~SubSampleSpatialObjectFilter( void );

  virtual void GenerateData( void ) override;

  /** Sub-sample at the at a given level, then sub-sample their
  * children. */
  virtual void SubSampleLevel( const SpatialObjectType * input,
    typename SpatialObjectType::Pointer output, bool graftOutput=false );

private:
  // purposely not implemented
  SubSampleSpatialObjectFilter( const Self & );

  // purposely not implemented
  void operator=( const Self & );

  SizeValueType m_Sampling;

}; // End class SubSampleSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSubSampleSpatialObjectFilter.hxx"
#endif

#endif // End !defined( __itktubeSubSampleSpatialObjectFilter_h )