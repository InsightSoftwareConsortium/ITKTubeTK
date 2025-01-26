/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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
#ifndef __tubeConvertTubesToTubeGraph_h
#define __tubeConvertTubesToTubeGraph_h

// ITK includes
#include "itkMacro.h"

#include "itkProcessObject.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeTubeSpatialObjectToTubeGraphFilter.h"

namespace tube
{
/** \class ConvertTubesToTubeGraph
 *
 *  \ingroup TubeTK
 */

template< class TPixel, unsigned int Dimension >
class ConvertTubesToTubeGraph:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertTubesToTubeGraph                    Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::TubeSpatialObjectToTubeGraphFilter
    < TPixel, Dimension > FilterType;

  typedef typename FilterType::InputImageType     InputImageType;
  typedef typename FilterType::TubeGroupType      TubeGroupType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro( ConvertTubesToTubeGraph);

  /** Set Number of Centroids */
  tubeWrapSetMacro( NumberOfCenteroids, int, Filter );
  tubeWrapGetMacro( NumberOfCenteroids, int, Filter );

  /* Set Central Voronoi Tesselation image */
  tubeWrapSetObjectMacro( CVTImage, InputImageType, Filter );
  tubeWrapGetObjectMacro( CVTImage, InputImageType, Filter );

  /* Set input tubes */
  tubeWrapSetObjectMacro( InputTubeGroup, TubeGroupType, Filter );
  tubeWrapGetObjectMacro( InputTubeGroup, TubeGroupType, Filter );

  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro( Filter );

  /** Get Adjacency Matrix */
  vnl_matrix< double > GetAdjacencyMatrix( void );

  /** Get Root Nodes Vector */
  vnl_vector< int > GetRootNodes( void );

  /** Get Branch Nodes Vector */
  vnl_vector< double > GetBranchNodes( void );

protected:
  ConvertTubesToTubeGraph( void );
  ~ConvertTubesToTubeGraph() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itktubeTubeSpatialObjectToTubeGraphFilter parameters **/
  ConvertTubesToTubeGraph( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typename FilterType::Pointer m_Filter;
};

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertTubesToTubeGraph.hxx"
#endif

#endif // End !defined( __tubeConvertTubesToTubeGraph_h )
