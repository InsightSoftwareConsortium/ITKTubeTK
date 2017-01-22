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

#ifndef __itktubeExtractTubePointsSpatialObjectFilter_h
#define __itktubeExtractTubePointsSpatialObjectFilter_h

#include <itkDataObjectDecorator.h>
#include <itkGroupSpatialObject.h>
#include <itkProcessObject.h>
#include <itkSpatialObject.h>

namespace itk
{

namespace tube
{
/** \class ExtractTubePointsSpatialObjectFilter
 *
 * \brief Extract all points in a tree of TubeSpatialObject's into an Array.
 *
 * ExtractTubePointsSpatialObjectFilter iterates through a tree of
 * TubeSpatialObject's ( or classes derived from
 * itk::TubeSpatialObject ), and places all the TubeSpatialObjectPointType
 * points
 * into an itk::VectorContainer< TubeSpatialObjectType::TubePointType >.
 * It is assumed that
 * This filter is templated over the type of the TubeSpatialObject.
 *
 * \sa GroupSpatialObject
 * \sa TubeSpatialObject
 * \sa TubeSpatialObjectPoint
 * \sa VesselTubeSpatialObject
 * \sa VesselTubeSpatialObjectPoint
 */
template< class TTubeSpatialObject >
class ExtractTubePointsSpatialObjectFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ExtractTubePointsSpatialObjectFilter Self;
  typedef ProcessObject                        Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  typedef TTubeSpatialObject                   TubeSpatialObjectType;

  typedef typename TTubeSpatialObject::TubePointType   TubePointType;

  typedef VectorContainer< IdentifierType, TubePointType > PointsContainerType;

  typedef DataObjectDecorator< PointsContainerType >
                                               PointsContainerDecoratorType;

  typedef GroupSpatialObject< TubeSpatialObjectType::ObjectDimension >
                                               GroupSpatialObjectType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ExtractTubePointsSpatialObjectFilter, ProcessObject );

  /** Standard New method. */
  itkNewMacro( Self );

  /** Set/Get the input GroupSpatialObject. */
  virtual void SetInput( const GroupSpatialObjectType * group );
  const GroupSpatialObjectType * GetInput( void ) const;

  /** Get the output PointsContainerType decorated as an itk::DataObject. */
  PointsContainerDecoratorType *       GetPointsContainerOutput( void );
  const PointsContainerDecoratorType * GetPointsContainerOutput( void ) const;

  /** Get the output PointsContainerType. */
  const PointsContainerType * GetPointsContainer( void ) const;

  using Superclass::MakeOutput;
  virtual ProcessObject::DataObjectPointer
    MakeOutput( ProcessObject::DataObjectPointerArraySizeType idx );

protected:
  ExtractTubePointsSpatialObjectFilter( void );
  virtual ~ExtractTubePointsSpatialObjectFilter( void );

  virtual void GenerateData();

private:
  // purposely not implemented
  ExtractTubePointsSpatialObjectFilter( const Self & );

  // purposely not implemented
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  typename PointsContainerDecoratorType::Pointer m_PointsContainerDecorator;
  typename PointsContainerType::Pointer          m_PointsContainer;

}; // End class ExtractTubePointsSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeExtractTubePointsSpatialObjectFilter.hxx"
#endif

#endif // End !defined( __itktubeExtractTubePointsSpatialObjectFilter_h )
