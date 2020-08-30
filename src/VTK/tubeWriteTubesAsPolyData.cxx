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

#include "tubeWriteTubesAsPolyData.h"

#include "tubeMessage.h"

#include <itkGroupSpatialObject.h>
#include <itkTubeSpatialObject.h>
#include <itksys/SystemTools.hxx>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyLine.h>
#include <vtkSTLWriter.h>
#include <vtkTubeFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

namespace tube
{

WriteTubesAsPolyData::
WriteTubesAsPolyData()
{
  m_GroupSpatialObject = nullptr;

  m_FileName = "";
  m_CenterlineFileName = "";

  m_NumberOfSides = 5;
}

/** Main work happens here */
void WriteTubesAsPolyData::
Update()
{
  const unsigned int Dimension = 3;

  typedef itk::TubeSpatialObject< Dimension >       TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >      GroupSpatialObjectType;

  m_GroupSpatialObject->Update();

  char childName[] = "Tube";
  typedef TubeSpatialObjectType::ChildrenListType ChildrenListType;
  ChildrenListType * tubeList =
    m_GroupSpatialObject->GetChildren( m_GroupSpatialObject->GetMaximumDepth(),
      childName );

  itk::SizeValueType totalNumberOfPoints = 0;
  for( ChildrenListType::iterator tubeIt = tubeList->begin();
      tubeIt != tubeList->end();
      ++tubeIt )
    {
    TubeSpatialObjectType * tube =
      static_cast< TubeSpatialObjectType * >( tubeIt->GetPointer() );

    tube->RemoveDuplicatePointsInObjectSpace();

    const itk::SizeValueType numberOfPoints = tube->GetNumberOfPoints();
    if( numberOfPoints < 2 )
      {
      continue;
      }

    totalNumberOfPoints += numberOfPoints;
    }

  // Create the points
  vtkNew< vtkPoints > tubeSpatialPoints;
  tubeSpatialPoints->SetNumberOfPoints( totalNumberOfPoints );

  // Create the Lines
  vtkNew< vtkCellArray > tubeLines;

  // Create scalar array that indicates the radius at each
  // centerline point.
  vtkNew< vtkDoubleArray > tubeRadius;
  tubeRadius->SetName( "Radius" );
  tubeRadius->SetNumberOfTuples( totalNumberOfPoints );

  // Create scalar array that indicates TubeId.
  vtkNew< vtkDoubleArray > tubeIds;
  tubeIds->SetName( "Id" );
  tubeIds->SetNumberOfTuples( totalNumberOfPoints );

  // Create scalar array that indicates both tangents at each
  // centerline point.
  vtkNew< vtkDoubleArray > normal1;
  normal1->SetName( "Normal1" );
  normal1->SetNumberOfTuples( Dimension * totalNumberOfPoints );
  normal1->SetNumberOfComponents( Dimension );

  vtkNew< vtkDoubleArray > normal2;
  normal2->SetName( "Normal2" );
  normal2->SetNumberOfTuples( Dimension * totalNumberOfPoints );
  normal2->SetNumberOfComponents( Dimension );

  bool containsAlpha1Info = false;
  vtkNew< vtkDoubleArray > alpha1;
  alpha1->SetName( "Alpha1" );
  alpha1->SetNumberOfTuples( totalNumberOfPoints );

  bool containsAlpha2Info = false;
  vtkNew< vtkDoubleArray > alpha2;
  alpha2->SetName( "Alpha2" );
  alpha2->SetNumberOfTuples( totalNumberOfPoints );

  bool containsAlpha3Info = false;
  vtkNew< vtkDoubleArray > alpha3;
  alpha3->SetName( "Alpha3" );
  alpha3->SetNumberOfTuples( totalNumberOfPoints );

  bool containsMedialnessInfo = false;
  vtkNew< vtkDoubleArray > medialness;
  medialness->SetName( "Medialness" );
  medialness->SetNumberOfTuples( totalNumberOfPoints );

  bool containsRidgenessInfo = false;
  vtkNew< vtkDoubleArray > ridgeness;
  ridgeness->SetName( "Ridgeness" );
  ridgeness->SetNumberOfTuples( totalNumberOfPoints );

  bool containsBranchnessInfo = false;
  vtkNew< vtkDoubleArray > branchness;
  branchness->SetName( "Branchness" );
  branchness->SetNumberOfTuples( totalNumberOfPoints );

  bool containsIntensityInfo = false;
  vtkNew< vtkDoubleArray > intensity;
  intensity->SetName( "Intensity" );
  intensity->SetNumberOfTuples( totalNumberOfPoints );

  bool containsCurvatureInfo = false;
  vtkNew< vtkDoubleArray > curvature;
  curvature->SetName( "Curvature" );
  curvature->SetNumberOfTuples( totalNumberOfPoints );

  bool containsRoundnessInfo = false;
  vtkNew< vtkDoubleArray > roundness;
  roundness->SetName( "Roundness" );
  roundness->SetNumberOfTuples( totalNumberOfPoints );

  bool containsLevelnessInfo = false;
  vtkNew< vtkDoubleArray > levelness;
  levelness->SetName( "Levelness" );
  levelness->SetNumberOfTuples( totalNumberOfPoints );

  std::list<vtkDoubleArray *> extraPointData;
  TubeSpatialObjectType * tube =
    static_cast< TubeSpatialObjectType * >( tubeList->begin()->GetPointer() );
  unsigned int numExtra =
    tube->GetPoints().begin()->GetTagScalarDictionary().size();
  auto pntDictIt = tube->GetPoints().begin()->GetTagScalarDictionary().begin();
  auto pntDictItEnd = tube->GetPoints().begin()->GetTagScalarDictionary().end();
  vtkDoubleArray * extraPointField;
  while( pntDictIt != pntDictItEnd )
    {
    extraPointField = vtkDoubleArray::New();
    extraPointField->SetName( pntDictIt->first.c_str() );
    extraPointField->SetNumberOfTuples( totalNumberOfPoints );
    extraPointData.push_back( extraPointField );
    ++pntDictIt;
    }

  itk::SizeValueType pointId = 0;
  for( ChildrenListType::iterator tubeIt = tubeList->begin();
      tubeIt != tubeList->end();
      ++tubeIt )
    {
    tube = static_cast< TubeSpatialObjectType * >( tubeIt->GetPointer() );

    const itk::SizeValueType numberOfPoints = tube->GetNumberOfPoints();
    if( numberOfPoints < 2 )
      {
      continue;
      }

    tube->RemoveDuplicatePointsInObjectSpace();

    tube->ComputeTangentsAndNormals();

    // Create a pointID list [linear for a polyline]
    vtkIdType * pointIds = new vtkIdType[numberOfPoints];
    vtkNew<vtkPolyLine> tubeLine;

    typedef TubeSpatialObjectType::TubePointListType::const_iterator
      TubePointIteratorType;
    const TubePointIteratorType tubePointsEnd = tube->GetPoints().end();
    itk::SizeValueType index = 0;
    for( TubePointIteratorType pointIt = tube->GetPoints().begin();
      pointIt != tubePointsEnd;
      ++pointIt, ++pointId, ++index )
      {
      TubeSpatialObjectType::PointType point =
        pointIt->GetPositionInWorldSpace();
      pointIds[index] = pointId;

      // Insert points using the element spacing information.
      tubeSpatialPoints->SetPoint( pointId, -1 * point[0],
        -1 * point[1], point[2] );
      // TubeId
      tubeIds->SetTuple1( pointId, tube->GetId() );

      // Radius
      tubeRadius->SetTuple1( pointId, pointIt->GetRadiusInWorldSpace() );

      // Tangeantes
      normal1->SetTuple3( pointId, pointIt->GetNormal1InWorldSpace()[0],
        pointIt->GetNormal1InWorldSpace()[1],
        pointIt->GetNormal1InWorldSpace()[2] );

      normal2->SetTuple3( pointId, pointIt->GetNormal2InWorldSpace()[0],
        pointIt->GetNormal2InWorldSpace()[1],
        pointIt->GetNormal2InWorldSpace()[2] );

      double tf = pointIt->GetRidgeness();
      if( tf != 0.0 )
        {
        containsRidgenessInfo = true;
        }
      ridgeness->SetTuple1( pointId, tf );

      tf = pointIt->GetBranchness();
      if( tf != 0.0 )
        {
        containsBranchnessInfo = true;
        }
      branchness->SetTuple1( pointId, tf );

      tf = pointIt->GetMedialness();
      if( tf != 0.0 )
        {
        containsMedialnessInfo = true;
        }
      medialness->SetTuple1( pointId, tf );

      tf = pointIt->GetCurvature();
      if( tf != 0.0 )
        {
        containsCurvatureInfo = true;
        }
      curvature->SetTuple1( pointId, tf );

      tf = pointIt->GetIntensity();
      if( tf != 0.0 )
        {
        containsIntensityInfo = true;
        }
      intensity->SetTuple1( pointId, tf );

      tf = pointIt->GetRoundness();
      if( tf != 0.0 )
        {
        containsRoundnessInfo = true;
        }
      roundness->SetTuple1( pointId, tf );

      tf = pointIt->GetLevelness();
      if( tf != 0.0 )
        {
        containsLevelnessInfo = true;
        }
      levelness->SetTuple1( pointId, tf );

      tf = pointIt->GetAlpha1();
      if( tf != 0.0 )
        {
        containsAlpha1Info = true;
        }
      alpha1->SetTuple1( pointId, tf );

      tf = pointIt->GetAlpha2();
      if( tf != 0.0 )
        {
        containsAlpha2Info = true;
        }
      alpha2->SetTuple1( pointId, tf );

      tf = pointIt->GetAlpha3();
      if( tf != 0.0 )
        {
        containsAlpha3Info = true;
        }
      alpha3->SetTuple1( pointId, tf );

      auto extraIt = extraPointData.begin();
      auto pntExtraIt = pointIt->GetTagScalarDictionary().begin();
      auto pntExtraItEnd = pointIt->GetTagScalarDictionary().end();
      while( pntExtraIt != pntExtraItEnd )
        {
        if( (*extraIt)->GetName() != pntExtraIt->first )
          {
          std::cerr << "Error: point tagscalar dictionary don't match: point = "
            << pntExtraIt->first << " != " << (*extraIt)->GetName()
            << std::endl;
          }
        (*extraIt)->SetTuple1( pointId, pntExtraIt->second );
        ++extraIt;
        ++pntExtraIt;
        }
      }


    tubeLine->Initialize( numberOfPoints,
                          pointIds,
                          tubeSpatialPoints.GetPointer() );

    tubeLines->InsertNextCell( tubeLine.GetPointer() );
    delete[] pointIds;
    }
  delete tubeList;

  // Convert spatial objects to a PolyData
  vtkNew< vtkPolyData > tubesCenterlineData;
  tubesCenterlineData->SetLines( tubeLines.GetPointer() );
  tubesCenterlineData->SetPoints( tubeSpatialPoints.GetPointer() );

  // Add the Radius information
  tubesCenterlineData->GetPointData()->AddArray( tubeRadius.GetPointer() );
  tubesCenterlineData->GetPointData()->SetActiveScalars( "Radius" );

  // Add the TudeId information
  tubesCenterlineData->GetPointData()->AddArray( tubeIds.GetPointer() );

  // Add Tangeantes information
  tubesCenterlineData->GetPointData()->AddArray( normal1.GetPointer() );
  tubesCenterlineData->GetPointData()->AddArray( normal2.GetPointer() );

  // Add Medialness & Ridgness if contains information
  if( containsIntensityInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( intensity.GetPointer() );
    }

  if( containsRidgenessInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( ridgeness.GetPointer() );
    }

  if( containsBranchnessInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( branchness.GetPointer() );
    }

  if( containsMedialnessInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( medialness.GetPointer() );
    }

  if( containsRoundnessInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( roundness.GetPointer() );
    }

  if( containsCurvatureInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( curvature.GetPointer() );
    }

  if( containsLevelnessInfo == true )
    {
    tubesCenterlineData->GetPointData()->AddArray( levelness.GetPointer() );
    }

  auto extraIt = extraPointData.begin();
  auto extraItEnd = extraPointData.end();
  while( extraIt != extraItEnd )
    {
    tubesCenterlineData->GetPointData()->AddArray( (*extraIt) );
    ++extraIt;
    }

  if( !m_CenterlineFileName.empty() )
    {
    vtkSmartPointer<vtkXMLPolyDataWriter> centerlineVTKwriter =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    centerlineVTKwriter->SetInputData( tubesCenterlineData.GetPointer() );
    centerlineVTKwriter->SetFileName( m_CenterlineFileName.c_str() );
    centerlineVTKwriter->Write();
    }

  vtkNew< vtkTubeFilter > tubesSurfaceFilter;
  tubesSurfaceFilter->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
  tubesSurfaceFilter->CappingOn();
  tubesSurfaceFilter->SetNumberOfSides( m_NumberOfSides );
  tubesSurfaceFilter->SetInputData( tubesCenterlineData.GetPointer() );
  tubesSurfaceFilter->Update();

  const std::string extension =
    itksys::SystemTools::GetFilenameLastExtension( m_FileName );

  if( extension == ".vtp" )
    {
    vtkNew< vtkXMLPolyDataWriter > polyDataWriter;
    polyDataWriter->SetInputConnection( tubesSurfaceFilter->GetOutputPort() );
    polyDataWriter->SetFileName( m_FileName.c_str() );
    polyDataWriter->SetDataModeToBinary();
    polyDataWriter->Write();
    }
  else if( extension == ".stl" )
    {
    // STL files only write triangles
    vtkNew< vtkTriangleFilter > triangleFilter;
    triangleFilter->SetInputConnection( tubesSurfaceFilter->GetOutputPort() );

    vtkNew< vtkSTLWriter > stlWriter;
    stlWriter->SetInputConnection( triangleFilter->GetOutputPort() );
    stlWriter->SetFileName( m_FileName.c_str() );
    stlWriter->Write();
    }
  else
    {
    std::cerr << "Unrecognized output file extension: " << extension
      << std::endl;
    std::cerr << "Appending .vtp instead" << std::endl;
    std::string filenamevtp = m_FileName + ".vtp";
    vtkNew< vtkXMLPolyDataWriter > polyDataWriter;
    polyDataWriter->SetInputConnection( tubesSurfaceFilter->GetOutputPort() );
    polyDataWriter->SetFileName( filenamevtp.c_str() );
    polyDataWriter->SetDataModeToBinary();
    polyDataWriter->Write();
    }

}

void WriteTubesAsPolyData::
PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Input group = " << m_GroupSpatialObject << std::endl;
}

}; //namespace
