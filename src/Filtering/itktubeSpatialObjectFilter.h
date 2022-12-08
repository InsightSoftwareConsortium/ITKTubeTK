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

#ifndef __itktubeSpatialObjectFilter_h
#define __itktubeSpatialObjectFilter_h

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

namespace itk
{

namespace tube
{
/** \class SpatialObjectFilter
 *
 * \brief Base class for filters that take a SpatialObject as input
 * and produce a SpatialObject as output.
 */
template< unsigned int ObjectDimension >
class SpatialObjectFilter
  : public SpatialObjectToSpatialObjectFilter< SpatialObject<ObjectDimension>,
    SpatialObject<ObjectDimension> >
{
public:
  /** Standard class typedefs */
  typedef SpatialObjectFilter                                  Self;
  typedef SpatialObjectToSpatialObjectFilter< SpatialObject<ObjectDimension >,
          SpatialObject<ObjectDimension> >                     Superclass;
  typedef SmartPointer< Self >                                 Pointer;
  typedef SmartPointer< const Self >                           ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SpatialObjectFilter, SpatialObjectSource );

  typedef SpatialObject<ObjectDimension>  SpatialObjectType;

  virtual void SetInput( const SpatialObjectType * spatialObject ) override;

  using Superclass::MakeOutput;
  virtual ProcessObject::DataObjectPointer
    MakeOutput( ProcessObject::DataObjectPointerArraySizeType idx ) override;

protected:
  SpatialObjectFilter( void );
  virtual ~SpatialObjectFilter( void ) {}

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const typename Superclass::DataObjectIdentifierType &,
    itk::DataObject * ) override {};

private:
  // purposely not implemented
  SpatialObjectFilter( const Self & );

  // purposely not implemented
  void operator=( const Self & );

}; // End class SpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSpatialObjectFilter.hxx"
#endif

#endif // End !defined( __itktubeSpatialObjectFilter_h )
