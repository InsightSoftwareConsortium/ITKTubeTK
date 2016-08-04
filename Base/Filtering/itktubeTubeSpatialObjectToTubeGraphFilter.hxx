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

#ifndef __itktubeTubeSpatialObjectToTubeGraphFilter_hxx
#define __itktubeTubeSpatialObjectToTubeGraphFilter_hxx

#include "itktubeTubeSpatialObjectToTubeGraphFilter.h"

/** Constructor */
template< class TPixel, unsigned int Dimension >
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::TubeSpatialObjectToTubeGraphFilter( void )
{
  m_CVTImage = NULL;
}

/** Destructor */
template< class TPixel, unsigned int Dimension >
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::~TubeSpatialObjectToTubeGraphFilter( void )
{
}

/** Get Adjacency Matrix */
template< class TPixel, unsigned int Dimension >
vnl_matrix< double >
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::GetAdjacencyMatrix( void )
{
  return m_AdjacencyMatrix;
}

/** Get Root Nodes Vector */
template< class TPixel, unsigned int Dimension >
vnl_vector< int >
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::GetRootNodes( void )
{
  return m_RootNodes;
}

/** Get Branch Nodes Vector */
template< class TPixel, unsigned int Dimension >
vnl_vector< double >
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::GetBranchNodes( void )
{
  return m_BranchNodes;
}

/** Update */
template< class TPixel, unsigned int Dimension >
void
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::Update( void )
{
  itkDebugMacro( << "TubeSpatialObjectToTubeGraphFilter::Update() called." );

  typedef itk::MinimumMaximumImageFilter< InputImageType > MinMaxFilterType;
  MinMaxFilterType::Pointer mmFilter = MinMaxFilterType::New();
  mmFilter->SetInput( m_CVTImage );
  mmFilter->Update();

  m_NumberOfCenteroids = mmFilter->GetMaximum();

  m_AdjacencyMatrix.set_size( m_NumberOfCenteroids, m_NumberOfCenteroids );
  m_AdjacencyMatrix.fill( 0 );
  m_RootNodes.set_size( m_NumberOfCenteroids );
  m_RootNodes.fill( 0 );
  m_BranchNodes.set_size( m_NumberOfCenteroids );
  m_BranchNodes.fill( 0 );

  vnl_matrix<double> cMat( Dimension, Dimension );
  vnl_vector<double> cVect( Dimension );

  char tubeName[] = "Tube";
  TubeSpatialObjectType::ChildrenListType * tubeList =
    m_InputTubeGroup->GetChildren
      ( m_InputTubeGroup->GetMaximumDepth(), tubeName );
  TubeSpatialObjectType::ChildrenListType::const_iterator
      tubeIt = tubeList->begin();
  int numTubes = tubeList->size();
  TubePointType tubePoint;
  MetaScene scene( Dimension );
  MetaTubeGraph * graph;
  TubeTransformType::Pointer tubeTransform;
  while( tubeIt != tubeList->end() )
    {
    TubeSpatialObjectType::Pointer tube =
        dynamic_cast< TubeSpatialObjectType * >(( *tubeIt ).GetPointer() );

    tube->RemoveDuplicatePoints();
    tube->ComputeTangentAndNormals();

    int numberOfPoints = tube->GetNumberOfPoints();

    graph = new MetaTubeGraph( Dimension );

    itk::Point< double, Dimension > pnt;
    itk::Index< Dimension > indx;
    tubePoint = static_cast< TubePointType >( tube->GetPoints()[0] );
    pnt = tubePoint.GetPosition();
    tube->ComputeObjectToWorldTransform();
    tubeTransform = tube->GetIndexToWorldTransform();
    pnt = tubeTransform->TransformPoint( pnt );
    m_CVTImage->TransformPhysicalPointToIndex( pnt, indx );
    double cCount = 1;
    int cNode = m_CVTImage->GetPixel( indx );
    double cRadius = tubePoint.GetRadius();
    for( int i = 0; i < Dimension; i++ )
      {
      cVect[i] = tubePoint.GetTangent()[i];
      }
    cMat = outer_product( cVect, cVect );
    if( tube->GetRoot() )
      {
      m_RootNodes[cNode - 1] = m_RootNodes[cNode - 1] + 1;
      }
    m_BranchNodes[cNode - 1] = m_BranchNodes[cNode - 1] + 1.0 / numTubes;
    int numberOfNodesCrossed = 0;
    for( int p = 1; p < numberOfPoints; p++ )
      {
      tubePoint = static_cast< TubePointType >( tube->GetPoints()[p] );
      pnt = tubePoint.GetPosition();
      pnt = tubeTransform->TransformPoint( pnt );
      m_CVTImage->TransformPhysicalPointToIndex( pnt, indx );
      int tNode = m_CVTImage->GetPixel( indx );
      if( tNode == cNode )
        {
        cCount++;
        cRadius += tubePoint.GetRadius();
        for( int i = 0; i < Dimension; i++ )
          {
          cVect[i] = tubePoint.GetTangent()[i];
          }
        cMat = cMat + outer_product( cVect, cVect );
        }
      else
        {
        int len = graph->GetPoints().size();
        if( graph->GetPoints().size() > 3
          && graph->GetPoints().at( len - 1 )->m_GraphNode == tNode
          && graph->GetPoints().at( len - 2 )->m_GraphNode == cNode )
          {
          itkDebugMacro( << " " );
          itkDebugMacro(  << "Oscillation detected"
                  << " : tube = " << cNode
                  << " : seq = " << graph->GetPoints().at( len - 3 )->m_GraphNode
                  << " " << graph->GetPoints().at( len - 2 )->m_GraphNode
                  << " " << graph->GetPoints().at( len - 1 )->m_GraphNode
                  << " " << cNode << " " << tNode );

          TubeGraphPnt * tgP = graph->GetPoints().back();
          cNode = tNode;
          cRadius = tgP->m_R;
          for( int i = 0; i < Dimension; i++ )
            {
            for( int j = 0; j < Dimension; j++ )
              {
              cMat[i][j] = tgP->m_T[i * Dimension + j];
              }
            }
          cCount = tgP->m_P;
          graph->GetPoints().pop_back();
          /* Memory allocated for each element of list returned by
          graph->GetPoints() usually released when destructor of graph called,
          but since tgP is popped off back of list, memory would not be
          released without explicit delete. */
          delete tgP;
          }
        else
          {
          numberOfNodesCrossed++;
          m_AdjacencyMatrix[cNode - 1][tNode - 1] =
            m_AdjacencyMatrix[cNode - 1][tNode - 1] + 1;
          TubeGraphPnt * tgP = new TubeGraphPnt( Dimension );
          tgP->m_GraphNode = cNode;
          tgP->m_R = cRadius / cCount;
          tgP->m_P = cCount;
          for( int i = 0; i < Dimension; i++ )
            {
            for( int j = 0; j < Dimension; j++ )
              {
              tgP->m_T[i * Dimension + j] = cMat[i][j] / cCount;
              }
            }
          graph->GetPoints().push_back( tgP );
          cNode = tNode;
          cRadius = tubePoint.GetRadius();
          for( int i = 0; i < Dimension; i++ )
            {
            cVect[i] = tubePoint.GetTangent()[i];
            }
          cMat = outer_product( cVect, cVect );
          cCount = 1;
          }
        }
      }
    if( numberOfNodesCrossed > 0 )
      {
      TubeGraphPnt * tgP = new TubeGraphPnt( Dimension );
      tgP->m_GraphNode = cNode;
      tgP->m_R = cRadius / cCount;
      for( int i = 0; i < Dimension; i++ )
        {
        for( int j = 0; j < Dimension; j++ )
          {
          tgP->m_T[i * Dimension + j] = cMat[i][j] / cCount;
          }
        }
      graph->GetPoints().push_back( tgP );
      scene.AddObject( graph );
      }
    else
      {
      delete graph;
      }
    ++tubeIt;
    }

  delete tubeList;

  itkDebugMacro( << "TubeSpatialObjectToTubeGraphFilter::Update() finished." );

} // End update function

template< class TPixel, unsigned int Dimension >
void
TubeSpatialObjectToTubeGraphFilter< TPixel, Dimension >
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Number of Centroids: "  << m_NumberOfCenteroids << std::endl;
}

#endif // End !defined(__itktubeTubeSpatialObjectToTubeGraphFilter_hxx)
