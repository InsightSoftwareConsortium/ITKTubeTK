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

#ifndef __itktubeMinimumSpanningTreeVesselConnectivityFilter_hxx
#define __itktubeMinimumSpanningTreeVesselConnectivityFilter_hxx

#include "itktubeMinimumSpanningTreeVesselConnectivityFilter.h"
#include "tubeMatrixMath.h"
#include "itkMath.h"
#include "tubeMacro.h"
#include "tubeTubeMath.h"

#include <utility>
#include <algorithm>
#include <iomanip>

namespace itk
{
namespace tube
{

//----------------------------------------------------------------------------
template< unsigned int VDimension >
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::MinimumSpanningTreeVesselConnectivityFilter( void )
{
  m_MaxTubeDistanceToRadiusRatio = 2.0;
  m_MaxContinuityAngleError = 180.0; // degrees
  m_RemoveOrphanTubes = false;
  m_numOutputConnectedComponents = 0;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::~MinimumSpanningTreeVesselConnectivityFilter( void )
{
  m_RootTubeIdList.clear();
  m_TubeGraph.clear();
  m_TubeIdToObjectMap.clear();
}

template< unsigned int VDimension >
bool
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::GraphEdgeType
::operator>( const GraphEdgeType & rhs ) const
{
  return weight > rhs.weight;
}

template< unsigned int VDimension >
bool
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::ConnectionPointType
::operator>( const ConnectionPointType & rhs ) const
{
  if( dist != rhs.dist )
    {
    return dist > rhs.dist;
    }
  else
    {
    return angle > rhs.angle;
    }
}

template< unsigned int VDimension >
bool
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::TubePQElementType
::operator<( const TubePQElementType & rhs ) const
{
  if( outDegree != rhs.outDegree )
    {
    return outDegree < rhs.outDegree;
    }
  else
    {
    return tubeLength < rhs.tubeLength;
    }
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::SetRootTubeIdList( const TubeIdListType & rootTubeIdList )
{
  m_RootTubeIdList = rootTubeIdList;
}

template< unsigned int VDimension >
const typename MinimumSpanningTreeVesselConnectivityFilter< VDimension >::
TubeIdListType &
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::GetRootTubeIdList( void ) const
{
  return m_RootTubeIdList;
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::BuildTubeGraph( void )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();

  // build a map between tube id and tube object
  tubeDebugMacro( << "Computing Tube ID to Object Map" );

  typedef typename TubeGroupType::ChildrenListPointer TubeListPointerType;

  char tubeName[] = "Tube";
  TubeListPointerType pTubeList
    = inputTubeGroup->GetChildren(
    inputTubeGroup->GetMaximumDepth(), tubeName );

  m_TubeIdToObjectMap.clear();

  for( typename TubeGroupType::ChildrenListType::iterator
    itTubes = pTubeList->begin();
    itTubes != pTubeList->end(); ++itTubes )
    {
    m_TubeIdToObjectMap[( *itTubes )->GetId()]
      = dynamic_cast< TubeType * >( itTubes->GetPointer() );
    }

  // build graph
  tubeDebugMacro( << "Building tube graph" );

  typedef typename TubeType::TubePointListType TubePointListType;
  typedef typename TubeType::TubePointType     TubePointType;
  typedef typename TubeType::PointType         PositionType;
  typedef typename PositionType::VectorType    PositionVectorType;

  m_TubeGraph.clear();

  for( typename TubeGroupType::ChildrenListType::iterator
    itSourceTubes = pTubeList->begin();
    itSourceTubes != pTubeList->end(); ++itSourceTubes )
    {
    TubePointerType pCurSourceTube
      = dynamic_cast< TubeType * >( itSourceTubes->GetPointer() );
    TubeIdType curSourceTubeId = pCurSourceTube->GetId();
    TubePointListType sourcePointList = pCurSourceTube->GetPoints();

    m_TubeGraph[curSourceTubeId].clear();

    int curSourceTubePointId = 0;
    for( typename TubePointListType::const_iterator
      itSourcePoints = sourcePointList.begin();
      itSourcePoints != sourcePointList.end(); ++itSourcePoints )
      {
      TubePointType ptSource = *itSourcePoints;
      PositionVectorType ptSourcePos
        = ptSource.GetPositionInObjectSpace().GetVectorFromOrigin();

      for( typename TubeGroupType::ChildrenListType::iterator
        itTargetTubes = pTubeList->begin();
        itTargetTubes != pTubeList->end(); ++itTargetTubes )
        {
        TubePointerType curTargetTube
          = dynamic_cast< TubeType * >( itTargetTubes->GetPointer() );
        TubeIdType curTargetTubeId = curTargetTube->GetId();

        if( curSourceTubeId == curTargetTubeId )
          {
          continue;
          }

        TubePointListType targetPointList = curTargetTube->GetPoints();

        if( targetPointList.size() <= 1 )
          {
          continue;
          }

        int ptCandidateIdList[] = {0, ( int ) targetPointList.size() - 1};

        std::priority_queue< ConnectionPointType,
          std::vector< ConnectionPointType >,
          std::greater< ConnectionPointType > >
        minpqConnPoint;
        ConnectionPointType ePtConn;

        for( unsigned int i = 0; i < 2; i++ )
          {
          int curPtId = ptCandidateIdList[i];
          TubePointType ptCur = targetPointList[curPtId];

          PositionVectorType ptCurPos
            = ptCur.GetPositionInObjectSpace().GetVectorFromOrigin();

          PositionVectorType vecToCurPt = ptCurPos - ptSourcePos;

          // compute and check distance
          double curDist = vecToCurPt.GetNorm();

          if( curDist > m_MaxTubeDistanceToRadiusRatio
            * ptSource.GetRadiusInObjectSpace() )
            {
            continue;
            }

          // compute and check angular continuity
          PositionVectorType ptNextPos;

          if( curPtId == 0 )
            {
            ptNextPos = targetPointList[ curPtId + 1 ].GetPositionInObjectSpace()
              .GetVectorFromOrigin();
            }
            else
            {
            ptNextPos = targetPointList[ curPtId - 1 ].GetPositionInObjectSpace()
              .GetVectorFromOrigin();
            }

          PositionVectorType curVecToNextPt = ptNextPos - ptCurPos;

          vecToCurPt.Normalize();
          curVecToNextPt.Normalize();

          double curAngle = std::acos( vecToCurPt * curVecToNextPt );
          curAngle *= 180.0 / itk::Math::pi;

          if( curAngle > m_MaxContinuityAngleError )
            {
            continue;
            }

          // add to queue
          ePtConn.dist = curDist;
          ePtConn.angle = curAngle;
          ePtConn.pointId = curPtId;
          minpqConnPoint.push( ePtConn );
          }

        if( minpqConnPoint.empty() )
          {
          continue;
          }

        ePtConn = minpqConnPoint.top();

        while( !minpqConnPoint.empty() )
          {
          minpqConnPoint.pop();
          }

        GraphEdgeType e;
        e.sourceTube       = pCurSourceTube;
        e.sourceTubeId      = curSourceTubeId;
        e.sourceTubePointId = curSourceTubePointId;

        e.targetTube       = curTargetTube;
        e.targetTubeId      = curTargetTubeId;
        e.targetTubePointId = ePtConn.pointId;

        e.weight = ePtConn.dist;
        e.distToRadRatio = ePtConn.dist / ptSource.GetRadiusInObjectSpace();
        e.continuityAngleError = ePtConn.angle;

        // if edge to current target is present then update it, else add it
        if( m_TubeGraph[curSourceTubeId].find( curTargetTubeId )
          != m_TubeGraph[curSourceTubeId].end() )
          {
          GraphEdgeType eOld = m_TubeGraph[curSourceTubeId][curTargetTubeId];

          // add only if current weight is better
          if( e.weight < eOld.weight )
            {
            m_TubeGraph[curSourceTubeId][curTargetTubeId] = e;
            }
          }
        else
          {
          m_TubeGraph[curSourceTubeId][curTargetTubeId] = e;
          }

        }
      ++curSourceTubePointId;
      }
    }

  // print graph
  tubeDebugMacro( << "\nNumber of tubes = "
    << m_TubeIdToObjectMap.size() );

  for( typename TubeAdjacencyListGraphType::const_iterator
    itV = m_TubeGraph.begin(); itV != m_TubeGraph.end(); ++itV )
    {
    tubeDebugMacro( << "Id = " << itV->first
      << ", numPoints = "
      << m_TubeIdToObjectMap[itV->first]->GetNumberOfPoints()
      << ", Out-degree = "
      <<  m_TubeGraph[itV->first].size() );

    for( typename GraphEdgeListType::const_iterator
      itE = itV->second.begin(); itE != itV->second.end(); ++itE )
      {
      GraphEdgeType e = itE->second;
      tubeDebugMacro( << "  Id = " << e.targetTubeId
        << ", sourcePointId = " << e.sourceTubePointId
        << ", targetPointId = " << e.targetTubePointId
        << ", weight = " << std::setprecision( 3 ) << e.weight
        << ", distanceToRadiusRatio = "
        << std::setprecision( 3 ) << e.distToRadRatio
        << ", continuityAngleError = "
        << std::setprecision( 3 ) << e.continuityAngleError );
      }
    }

  pTubeList->clear();
  delete pTubeList;
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::VisitTube( TubePointerType pTube )
{
  TubeIdType tubeId = pTube->GetId();

  m_SetTubesVisited.insert( tubeId );

  // add neighbors to priority queue
  for( typename GraphEdgeListType::iterator
    itRootNeighbors = m_TubeGraph[tubeId].begin();
    itRootNeighbors != m_TubeGraph[tubeId].end(); ++itRootNeighbors )
    {
    TubeIdType curTargetTubId = itRootNeighbors->second.targetTubeId;

    // push neighbor to priority queue if not already visited
    if( m_SetTubesVisited.find( curTargetTubId ) == m_SetTubesVisited.end() )
      {
      itRootNeighbors->second.sourceTube = pTube;
      m_minpqGraphEdge.push( itRootNeighbors->second );
      }
    }
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::RunMinimumSpanningTree( TubeIdType rootTubeId )
{
  TubeGroupType * outputTubeGroup = this->GetOutput();

  // make sure if given root tube was not visited already
  if( m_SetTubesVisited.find( rootTubeId ) != m_SetTubesVisited.end() )
    {
    return;
    }

  // Add root tube to output
  TubePointerType inputRootTube = m_TubeIdToObjectMap[rootTubeId];
  TubePointerType rootTube = inputRootTube->Clone();

  // TODO: make CopyInformation of itk::SpatialObject do this
  rootTube->GetObjectToParentTransform()->SetFixedParameters(
    inputRootTube->GetObjectToParentTransform()->GetFixedParameters() );
  rootTube->GetObjectToParentTransform()->SetParameters(
    inputRootTube->GetObjectToParentTransform()->GetParameters() );
  rootTube->Update();

  rootTube->ComputeTangentAndNormals();
  rootTube->SetRoot( true );

  // visit root tube
  VisitTube( rootTube );
  // add root to output
  if( m_RemoveOrphanTubes && m_minpqGraphEdge.empty() )
    {
    bool isActualRoot = std::find( m_RootTubeIdList.begin(),
      m_RootTubeIdList.end(), rootTubeId ) != m_RootTubeIdList.end();
    if( !isActualRoot )
      {
      return;
      }
    }
  outputTubeGroup->AddChild( rootTube );

  // recusrively process all children in increasing order of connection weight
  int numChildren = 0;

  tubeDebugMacro( << "\nRunning MST on tube " << rootTubeId
    <<  ", with out-degree = "
    << m_TubeGraph[rootTubeId].size()
    <<  ", and numPoints = "
    << rootTube->GetNumberOfPoints() );

  while( !m_minpqGraphEdge.empty() )
    {
    // pop smallest weighted edge
    GraphEdgeType eTop = m_minpqGraphEdge.top();
    m_minpqGraphEdge.pop();

    // if target was already visited, do nothing
    if( m_SetTubesVisited.find( eTop.targetTubeId ) != m_SetTubesVisited.end() )
      {
      continue;
      }

    // get tube object
    TubePointerType curTube = eTop.targetTube->Clone();

    // TODO: make CopyInformation of itk::SpatialObject do this
    curTube->GetObjectToParentTransform()->SetFixedParameters(
      eTop.targetTube->GetObjectToParentTransform()->GetFixedParameters() );
    curTube->GetObjectToParentTransform()->SetParameters(
      eTop.targetTube->GetObjectToParentTransform()->GetParameters() );
    curTube->Update();

    curTube->SetRoot( false );

    // set parent tube id
    curTube->SetParentId( eTop.sourceTubeId );

    // correct parent point id if parent tube points were reversed
    if( m_SetTubesReversed.find( eTop.sourceTubeId )
      != m_SetTubesReversed.end() )
      {
      int numParentTubePoints = ( int ) eTop.sourceTube->GetNumberOfPoints();
      eTop.sourceTubePointId = numParentTubePoints - eTop.sourceTubePointId - 1;
      }

    // set parent tube point id
    curTube->SetParentPoint( eTop.sourceTubePointId );

    // reverse point list if parent wants to connect to last point
    if( eTop.targetTubePointId > 0 )
      {
      std::reverse( curTube->GetPoints().begin(), curTube->GetPoints().end() );
      m_SetTubesReversed.insert( eTop.targetTubeId );
      }

    curTube->ComputeTangentAndNormals();

    // add tube to output
    eTop.sourceTube->AddChild( curTube );
    // print some info
    tubeDebugMacro(
      << "  sourceTubeId = "    << eTop.sourceTubeId
      << ", targetTubeId = "    << eTop.targetTubeId
      << ", reversed points = " << ( eTop.targetTubePointId > 0 )
      << std::setprecision( 3 )
      << ", weight = "          << eTop.weight
      << std::setprecision( 3 )
      << ", continuity angle error = " << eTop.continuityAngleError );

    // visit tube
    VisitTube( curTube );
    ++numChildren;
    }

  tubeDebugMacro( << "  Number of children = " << numChildren );

  ++m_numOutputConnectedComponents;
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::ComputeTubeConnectivity( void )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();
  TubeGroupType * outputTubeGroup = this->GetOutput();

  // initialize output and copy metadata from input
  outputTubeGroup->Clear();

  outputTubeGroup->CopyInformation( inputTubeGroup );

  // TODO: make CopyInformation of itk::SpatialObject do this
  outputTubeGroup->GetObjectToParentTransform()->SetFixedParameters(
    inputTubeGroup->GetObjectToParentTransform()->GetFixedParameters() );
  outputTubeGroup->GetObjectToParentTransform()->SetParameters(
    inputTubeGroup->GetObjectToParentTransform()->GetParameters() );
  outputTubeGroup->Update();

  // initialize
  m_SetTubesVisited.clear();
  m_SetTubesReversed.clear();

  // prepare a list of root tube ids to seed from
  std::priority_queue< TubePQElementType > maxpqVOutDegree;

  if( !m_RootTubeIdList.empty() )
    {
    // push given root tubes into a priority queue with outDegree as priority
    for( typename TubeIdListType::const_iterator
      itRootTubeId = m_RootTubeIdList.begin();
      itRootTubeId != m_RootTubeIdList.end(); ++itRootTubeId )
      {
      // throw exception if current root tube id is not present
      if( m_TubeIdToObjectMap.find( *itRootTubeId )
        == m_TubeIdToObjectMap.end() )
        {
        itkExceptionMacro( << "Could not find root tube id: "
          << *itRootTubeId );
        }

      TubePQElementType epTube;
      epTube.tubeId = *itRootTubeId;
      epTube.outDegree = m_TubeGraph[epTube.tubeId].size();
      epTube.tubeLength = ::tube::ComputeTubeLength< TubeType >(
        m_TubeIdToObjectMap[epTube.tubeId] );
      maxpqVOutDegree.push( epTube );
      }
    }
  else
    {
    // push all tubes into a priority queue with outDegree as priotiry
    for( typename TubeAdjacencyListGraphType::const_iterator
      itV = m_TubeGraph.begin(); itV != m_TubeGraph.end(); ++itV )
      {
      TubePQElementType epTube;
      epTube.tubeId = itV->first;
      epTube.outDegree = itV->second.size();
      epTube.tubeLength = ::tube::ComputeTubeLength< TubeType >(
        m_TubeIdToObjectMap[epTube.tubeId] );
      if( m_RemoveOrphanTubes && epTube.outDegree == 0 )
        {
        continue;
        }
      maxpqVOutDegree.push( epTube );
      }
    }

  // Run Minimum spanning tree like algorith to find subtree rooted at
  // each tube
  // Note: the root tubes are explored in decreasing order of out-degree
  m_numOutputConnectedComponents = 0;

  while( !maxpqVOutDegree.empty() )
    {
    // Get next tube with largest out degree
    TubePQElementType eTop = maxpqVOutDegree.top();
    maxpqVOutDegree.pop();

    // Run MST
    RunMinimumSpanningTree( eTop.tubeId );
    }

  tubeDebugMacro( << "Number of output connected components = "
    <<  m_numOutputConnectedComponents );
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::AddRemainingTubes( void )
{
  const TubeGroupType * inputTubeGroup = this->GetInput();
  TubeGroupType * outputTubeGroup = this->GetOutput();

  tubeDebugMacro( << "Adding remaining tubes to the output spatial group" );

  typedef typename TubeGroupType::ChildrenListPointer TubeListPointerType;

  char tubeName[] = "Tube";
  TubeListPointerType pTubeList
    = inputTubeGroup->GetChildren(
    inputTubeGroup->GetMaximumDepth(), tubeName );

  for( typename TubeGroupType::ChildrenListType::iterator
       itSourceTubes = pTubeList->begin();
       itSourceTubes != pTubeList->end(); ++itSourceTubes )
    {
    TubePointerType pCurSourceTube
      = dynamic_cast< TubeType * >( itSourceTubes->GetPointer() );
    TubeIdType curSourceTubeId = pCurSourceTube->GetId();

    if( m_SetTubesVisited.find( curSourceTubeId ) == m_SetTubesVisited.end() )
      {
      TubePointerType curTube = pCurSourceTube->Clone();

      // TODO: make CopyInformation of itk::SpatialObject do this
      curTube->GetObjectToParentTransform()->SetFixedParameters(
        pCurSourceTube->GetObjectToParentTransform()->GetFixedParameters() );
      curTube->GetObjectToParentTransform()->SetParameters(
        pCurSourceTube->GetObjectToParentTransform()->GetParameters() );
      curTube->Update();
      curTube->ComputeTangentAndNormals();
      curTube->SetRoot( false );

      outputTubeGroup->AddChild( curTube );
      }
    }
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::GenerateData( void )
{
  BuildTubeGraph();
  ComputeTubeConnectivity();
  if( !m_RemoveOrphanTubes )
    {
    AddRemainingTubes();
    }
}

template< unsigned int VDimension >
void
MinimumSpanningTreeVesselConnectivityFilter< VDimension >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  // print parameters
  tubeStandardOutputMacro( << "\nMaxTubeDistanceToRadiusRatio = "
                           << m_MaxTubeDistanceToRadiusRatio );
  tubeStandardOutputMacro( << "\nMaxContinuityAngleError = "
                           << m_MaxContinuityAngleError );
  tubeStandardOutputMacro( << "\nRemove Orphan Tubes = "
                           << m_RemoveOrphanTubes );
}

} // End namespace tube

} // End namespace itk

// End !defined( __itktubeMinimumSpanningTreeVesselConnectivityFilter_hxx )
#endif
