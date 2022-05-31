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

#ifndef __itktubeSubSampleTubeSpatialObjectFilter_h
#define __itktubeSubSampleTubeSpatialObjectFilter_h

#include "itkTubeSpatialObject.h"
#include "itktubeSpatialObjectFilter.h"

namespace itk
{

namespace tube
{

/** \class SubSampleTubeSpatialObjectFilter
 * \brief Sub-sample points from a tube.
 *
 * The input tube are sub-sampled by the \c Sampling
 * ( an integer greater or equal to one ).  The beginning and end
 * points of the tube are always included.
 *
 * \sa SubSampleTubeTreeSpatialObjectFilter
 */
template< unsigned int ObjectDimension >
class SubSampleTubeSpatialObjectFilter
  : public SpatialObjectFilter< ObjectDimension >
{
public:
  /** Standard class typedefs. */
  typedef SubSampleTubeSpatialObjectFilter          Self;
  typedef SpatialObjectFilter< ObjectDimension >    Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  typedef TubeSpatialObject<ObjectDimension>        TubeSpatialObjectType;

  /** Run-time type information ( and related methods ).   */
  itkTypeMacro( SubSampleTubeSpatialObjectFilter,
    SpatialObjectFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Set the sampling factor.  The output points taken every sampling
   * factor from the input points. */
  itkSetClampMacro( Sampling, SizeValueType, 1,
    NumericTraits< SizeValueType >::max() );
  itkGetConstMacro( Sampling, SizeValueType );

protected:
  SubSampleTubeSpatialObjectFilter( void );
  virtual ~SubSampleTubeSpatialObjectFilter( void );

  virtual void GenerateData( void ) override;

private:
  // purposely not implemented
  SubSampleTubeSpatialObjectFilter( const Self & );
  // purposely not implemented
  void operator=( const Self & );

  SizeValueType m_Sampling;

}; // End class SubSampleTubeSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSubSampleTubeSpatialObjectFilter.hxx"
#endif

#endif // End !defined( __itktubeSubSampleTubeSpatialObjectFilter_h )
