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

#ifndef __itktubeTubeSpatialObjectToTubeGraphFilter_h
#define __itktubeTubeSpatialObjectToTubeGraphFilter_h

#include <itkMinimumMaximumImageFilter.h>
#include <itkGroupSpatialObject.h>
#include <itkVesselTubeSpatialObject.h>
#include <itkImage.h>
#include <metaScene.h>
#include <metaTubeGraph.h>

namespace itk
{

namespace tube
{

/** \class TubeSpatialObjectToTubeGraphFilter
 */

template< class TPixel, unsigned int Dimension >
class TubeSpatialObjectToTubeGraphFilter
  : public Object
{
public:

  /** Standard class typedefs. */
  typedef TubeSpatialObjectToTubeGraphFilter             Self;
  typedef Object                                         SuperClass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  /** Tube class typedef */
  typedef Image< TPixel, Dimension >                     InputImageType;
  typedef typename InputImageType::Pointer               InputImagePointer;
  typedef GroupSpatialObject< Dimension >                TubeGroupType;
  typedef typename TubeGroupType::Pointer                TubeGroupPointer;
  typedef VesselTubeSpatialObject< Dimension >           TubeSpatialObjectType;
  typedef typename TubeSpatialObjectType::TubePointType  TubePointType;
  typedef typename TubeSpatialObjectType::TransformType  TubeTransformType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TubeSpatialObjectToTubeGraphFilter,
                Object );

  /** Set Number of Centroids */
  itkSetMacro( NumberOfCenteroids, int );
  itkGetMacro( NumberOfCenteroids, int );

  /** Set Central Voronoi Tesselation Image*/
  itkSetObjectMacro( CVTImage, InputImageType );
  itkGetObjectMacro( CVTImage, InputImageType );

  /** Sets the input tubes */
  itkSetMacro( InputTubeGroup, TubeGroupPointer );
  itkGetMacro( InputTubeGroup, TubeGroupPointer );

  /** Get Adjacency Matrix */
  vnl_matrix< double > GetAdjacencyMatrix( void );

  /** Get Root Nodes Vector */
  vnl_vector< int > GetRootNodes( void );

  /** Get Branch Nodes Vector */
  vnl_vector< double > GetBranchNodes( void );

  void Update( void );
protected:

  TubeSpatialObjectToTubeGraphFilter( void );
  ~TubeSpatialObjectToTubeGraphFilter( void );

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  int                    m_NumberOfCenteroids;
  InputImagePointer      m_CVTImage;
  TubeGroupPointer       m_InputTubeGroup;
  vnl_matrix< double >   m_AdjacencyMatrix;
  vnl_vector< int >      m_RootNodes;
  vnl_vector< double >   m_BranchNodes;

}; // End class TubeSpatialObjectToTubeGraphFilter

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeSpatialObjectToTubeGraphFilter.hxx"
#endif

} // End namespace tube

} // End namespace itk

#endif // End !defined(__TubeSpatialObjectToTubeGraphFilter_h)
