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

#ifndef __itktubeSpatialObjectToSpatialObjectFilter_h
#define __itktubeSpatialObjectToSpatialObjectFilter_h

#include "itktubeSpatialObjectSource.h"

namespace itk
{

namespace tube
{
/** \class SpatialObjectToSpatialObjectFilter
 *
 * \brief Base class for filters that take a SpatialObject as input
 * and produce a SpatialObject as output.
 */
template< class TInputSpatialObject, class TOutputSpatialObject >
class SpatialObjectToSpatialObjectFilter
  : public SpatialObjectSource< TOutputSpatialObject >
{
public:
  /** Standard class typedefs */
  typedef SpatialObjectToSpatialObjectFilter          Self;
  typedef SpatialObjectSource< TOutputSpatialObject > Superclass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SpatialObjectToSpatialObjectFilter, SpatialObjectSource );

  typedef TInputSpatialObject  InputSpatialObjectType;
  typedef TOutputSpatialObject OutputSpatialObjectType;

  virtual void SetInput( const InputSpatialObjectType * spatialObject );

  virtual void SetInput( unsigned int, const InputSpatialObjectType *
    spatialObject );

  // from itkProcessObject
  virtual void SetInput( const
    itk::ProcessObject::DataObjectIdentifierType & index, itk::DataObject *
    input ) override;

  const InputSpatialObjectType * GetInput( void ) const;

  const InputSpatialObjectType * GetInput( unsigned int idx ) const;

protected:
  SpatialObjectToSpatialObjectFilter( void );
  virtual ~SpatialObjectToSpatialObjectFilter( void ) {}

  // To remove warning "was hidden [-Woverloaded-virtual]"
  //void SetInput( const typename Superclass::DataObjectIdentifierType &,
    //itk::DataObject * ) override {};

private:
  // purposely not implemented
  SpatialObjectToSpatialObjectFilter( const Self & );

  // purposely not implemented
  void operator=( const Self & );

}; // End class SpatialObjectToSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSpatialObjectToSpatialObjectFilter.hxx"
#endif

#endif // End !defined( __itktubeSpatialObjectToSpatialObjectFilter_h )
