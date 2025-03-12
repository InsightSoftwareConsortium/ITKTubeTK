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

#ifndef itktubeMinimumSpanningTreeVesselConnectivityFilter_h
#define itktubeMinimumSpanningTreeVesselConnectivityFilter_h

#include "itkGroupSpatialObject.h"
#include "itkImage.h"
#include "itktubeSpatialObjectToSpatialObjectFilter.h"
#include "itkTubeSpatialObject.h"
#include <unordered_map>

#include <map>
#include <queue>
#include <math.h>
#include <functional>
#include <vector>

namespace itk
{
namespace tube
{

/** \class MinimumSpanningTreeVesselConnectivityFilter
 * \brief Computes connectivity between tubes in a TRE file
 *
 */
template <unsigned int ObjectDimension>
class MinimumSpanningTreeVesselConnectivityFilter
  : public SpatialObjectToSpatialObjectFilter<GroupSpatialObject<ObjectDimension>, GroupSpatialObject<ObjectDimension>>
{
public:
  /** Standard class type alias. */
  using TubeGroupType = GroupSpatialObject<ObjectDimension>;

  typedef MinimumSpanningTreeVesselConnectivityFilter Self;
  using Superclass = SpatialObjectToSpatialObjectFilter<TubeGroupType, TubeGroupType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using TubeType = TubeSpatialObject<ObjectDimension>;
  using TubePointerType = typename TubeType::Pointer;
  using TubeConstPointerType = typename TubeType::ConstPointer;
  using TubeIdType = itk::IndexValueType;
  using TubeIdListType = std::vector<TubeIdType>;

  /** Run-time type information ( and related methods ).   */
  itkOverrideGetNameOfClassMacro(MinimumSpanningTreeVesselConnectivityFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get a list of root tube ids */
  void
  SetRootTubeIdList(const TubeIdListType & rootTubeIdList);
  const TubeIdListType &
  GetRootTubeIdList(void) const;

  /** Set/Get max tube distance */
  itkSetMacro(MaxTubeDistanceToRadiusRatio, double);
  itkGetConstMacro(MaxTubeDistanceToRadiusRatio, double);

  /** Set/Get bifurcation angle continuity */
  itkSetMacro(MaxContinuityAngleError, double);
  itkGetConstMacro(MaxContinuityAngleError, double);

  /** Set/Get whether or not to remove orphan tubes */
  itkSetMacro(RemoveOrphanTubes, bool);
  itkGetMacro(RemoveOrphanTubes, bool);
  itkBooleanMacro(RemoveOrphanTubes);

protected:
  MinimumSpanningTreeVesselConnectivityFilter(void);
  virtual ~MinimumSpanningTreeVesselConnectivityFilter(void);

  virtual void
  GenerateData(void) override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  // purposely not implemented
  MinimumSpanningTreeVesselConnectivityFilter(const Self &);
  // purposely not implemented
  void
  operator=(const Self &);

  double
  ComputeTubeLength();
  void
  BuildTubeGraph(void);
  void
  ComputeTubeConnectivity(void);
  void
  VisitTube(TubePointerType pTube);
  void
  RunMinimumSpanningTree(TubeIdType rootTubeId);
  void
  AddRemainingTubes();

  struct GraphEdgeType
  {
    double weight;
    double distToRadRatio;
    double continuityAngleError;

    TubePointerType sourceTube;
    TubeIdType      sourceTubeId;
    int             sourceTubePointId;

    TubePointerType targetTube;
    TubeIdType      targetTubeId;
    int             targetTubePointId;

    bool
    operator>(const GraphEdgeType & rhs) const;
  };

  struct ConnectionPointType
  {
    int    pointId;
    double dist;
    double angle;

    bool
    operator>(const ConnectionPointType & rhs) const;
  };

  using GraphEdgeListType = std::unordered_map<TubeIdType, GraphEdgeType>;
  using TubeAdjacencyListGraphType = std::unordered_map<TubeIdType, GraphEdgeListType>;
  using TubeIdToPointerMapType = std::unordered_map<TubeIdType, TubePointerType>;

  struct TubePQElementType
  {
    TubeIdType                            tubeId;
    typename GraphEdgeListType::size_type outDegree;
    double                                tubeLength;

    bool
    operator<(const TubePQElementType & rhs) const;
  };

  double         m_MaxTubeDistanceToRadiusRatio;
  double         m_MaxContinuityAngleError;
  TubeIdListType m_RootTubeIdList;
  bool           m_RemoveOrphanTubes;

  TubeAdjacencyListGraphType                                                                  m_TubeGraph;
  TubeIdToPointerMapType                                                                      m_TubeIdToObjectMap;
  std::set<TubeIdType>                                                                        m_SetTubesVisited;
  std::set<TubeIdType>                                                                        m_SetTubesReversed;
  std::priority_queue<GraphEdgeType, std::vector<GraphEdgeType>, std::greater<GraphEdgeType>> m_minpqGraphEdge;
  int m_numOutputConnectedComponents;

}; // End class MinimumSpanningTreeVesselConnectivityFilter

} // End namespace tube
} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeMinimumSpanningTreeVesselConnectivityFilter.hxx"
#endif

#endif // End !defined( __itktubeMinimumSpanningTreeVesselConnectivityFilter_h )
