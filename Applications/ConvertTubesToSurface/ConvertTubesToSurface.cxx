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

#include "tubeMessage.h"
#include "tubeCLIProgressReporter.h"

#include <itkSpatialObjectReader.h>
#include <itkGroupSpatialObject.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkVesselTubeSpatialObject.h>

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

#include "ConvertTubesToSurfaceCLP.h"

/** Forward declarations */
int DoIt( int argc, char * argv[] );

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}

/** Main work happens here */
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter(
    "ConvertTubesToSurface",
    CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  const unsigned int Dimension = 3;

  typedef itk::VesselTubeSpatialObject< Dimension > TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >      GroupSpatialObjectType;

  timeCollector.Start( "Read tubes" );
  typedef itk::SpatialObjectReader< Dimension >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTubeFile );
  try
    {
    reader->Update();
    }
  catch( const std::exception & error )
    {
    tube::ErrorMessage( error.what() );
    return EXIT_FAILURE;
    }
  GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  std::ostringstream ostrm;
  ostrm << "Number of children = "
    << groupSpatialObject->GetNumberOfChildren();
  tube::InformationMessage( ostrm.str() );
  timeCollector.Stop( "Read tubes" );

  progress = 0.3;
  progressReporter.Report( progress );

  timeCollector.Start( "Convert to surface" );
  char childName[] = "Tube";
  typedef TubeSpatialObjectType::ChildrenListType ChildrenListType;
  ChildrenListType * tubeList =
    groupSpatialObject->GetChildren( groupSpatialObject->GetMaximumDepth(),
                                     childName );

  itk::SizeValueType totalNumberOfPoints = 0;
  for(ChildrenListType::iterator tubeIt = tubeList->begin();
      tubeIt != tubeList->end();
      ++tubeIt )
    {
    TubeSpatialObjectType * tube =
      static_cast< TubeSpatialObjectType * >( tubeIt->GetPointer() );

    tube->RemoveDuplicatePoints();

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
  tubeRadius->SetName( "TubeRadius" );
  tubeRadius->SetNumberOfTuples( totalNumberOfPoints );

  // Create scalar array that indicates TubeId.
  vtkNew< vtkDoubleArray > tubeIds;
  tubeIds->SetName( "TubeIds" );
  tubeIds->SetNumberOfTuples( totalNumberOfPoints );

  // Create scalar array that indicates both tangents at each
  // centerline point.
  vtkNew< vtkDoubleArray > tan1;
  tan1->SetName("Tan1");
  tan1->SetNumberOfTuples( Dimension * totalNumberOfPoints );
  tan1->SetNumberOfComponents( Dimension );

  vtkNew< vtkDoubleArray > tan2;
  tan2->SetName( "Tan2" );
  tan2->SetNumberOfTuples( Dimension * totalNumberOfPoints);
  tan2->SetNumberOfComponents( Dimension );

  // Create scalar array that indicates Ridgeness and medialness at each
  // centerline point.
  bool containsMedialnessInfo = false;
  vtkNew< vtkDoubleArray > medialness;
  medialness->SetName( "Medialness" );
  medialness->SetNumberOfTuples( totalNumberOfPoints );

  bool containsRidgenessInfo = false;
  vtkNew< vtkDoubleArray > ridgeness;
  ridgeness->SetName( "Ridgeness" );
  ridgeness->SetNumberOfTuples( totalNumberOfPoints );

  itk::SizeValueType pointId = 0;
  for(ChildrenListType::iterator tubeIt = tubeList->begin();
      tubeIt != tubeList->end();
      ++tubeIt )
    {
    TubeSpatialObjectType * tube =
      static_cast< TubeSpatialObjectType * >( tubeIt->GetPointer() );

    const itk::SizeValueType numberOfPoints = tube->GetNumberOfPoints();
    if( numberOfPoints < 2 )
      {
      continue;
      }

    tube->ComputeTangentAndNormals();

    // Create a pointID list [linear for a polyline]
    vtkIdType * pointIds = new vtkIdType[numberOfPoints];
    vtkNew<vtkPolyLine> tubeLine;

    // Get the tube element spacing information.
    const double* axesRatio = tube->GetSpacing();

    const TubeSpatialObjectType::PointListType & tubePoints =
      tube->GetPoints();
    typedef TubeSpatialObjectType::PointListType::const_iterator TubePointIteratorType;
    const TubePointIteratorType tubePointsEnd = tubePoints.end();
    itk::SizeValueType index = 0;
    for( TubePointIteratorType pointIt = tubePoints.begin();
         pointIt != tubePointsEnd;
         ++pointIt, ++pointId, ++index )
      {
      TubeSpatialObjectType::PointType point = pointIt->GetPosition();
      pointIds[index] = pointId;

      // Insert points using the element spacing information.
      tubeSpatialPoints->SetPoint( pointId,
                              point[0],
                              point[1] * axesRatio[1] / axesRatio[0],
                              point[2] * axesRatio[2] / axesRatio[0] );

      // TubeId
      tubeIds->SetTuple1( pointId, tube->GetId() );

      // Radius
      tubeRadius->SetTuple1( pointId, pointIt->GetRadius() );

      // Tangeantes
      tan1->SetTuple3( pointId,
                       pointIt->GetNormal1()[0],
                       pointIt->GetNormal1()[1],
                       pointIt->GetNormal1()[2] );

      tan2->SetTuple3( pointId,
                       pointIt->GetNormal2()[0],
                       pointIt->GetNormal2()[1],
                       pointIt->GetNormal2()[2] );

      // Medialness & Ridgness
      if( pointIt->GetMedialness() != 0.0 )
        {
        containsMedialnessInfo = true;
        }
      medialness->SetTuple1( pointId, pointIt->GetMedialness() );

      if( pointIt->GetRidgeness() != 0.0 )
        {
        containsRidgenessInfo = true;
        }
      ridgeness->SetTuple1( pointId, pointIt->GetRidgeness() );
      }

    tubeLine->Initialize( numberOfPoints,
                            pointIds,
                            tubeSpatialPoints.GetPointer() );

    tubeLines->InsertNextCell( tubeLine.GetPointer() );
    delete[] pointIds;
    }
  delete tubeList;

  // Convert spatial objects to a PolyData
  vtkNew< vtkPolyData > tubesPolyData;
  tubesPolyData->SetLines( tubeLines.GetPointer() );
  tubesPolyData->SetPoints( tubeSpatialPoints.GetPointer() );

  // Add the Radius information
  tubesPolyData->GetPointData()->AddArray( tubeRadius.GetPointer() );
  tubesPolyData->GetPointData()->SetActiveScalars( "TubeRadius" );

  // Add the TudeId information
  tubesPolyData->GetPointData()->AddArray( tubeIds.GetPointer() );

  // Add Tangeantes information
  tubesPolyData->GetPointData()->AddArray( tan1.GetPointer() );
  tubesPolyData->GetPointData()->AddArray( tan2.GetPointer() );

  // Add Medialness & Ridgness if contains information
  if( containsMedialnessInfo == true )
    {
    tubesPolyData->GetPointData()->AddArray( medialness.GetPointer() );
    }

  if( containsRidgenessInfo == true )
    {
    tubesPolyData->GetPointData()->AddArray( ridgeness.GetPointer() );
    }

  vtkNew< vtkTubeFilter > tubeSurfaceFilter;
  tubeSurfaceFilter->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
  tubeSurfaceFilter->CappingOn();
  tubeSurfaceFilter->SetNumberOfSides( 32 );
#if VTK_MAJOR_VERSION <= 5
  tubeSurfaceFilter->SetInput( tubesPolyData.GetPointer() );
#else
  tubeSurfaceFilter->SetInputData( tubesPolyData.GetPointer() );
#endif
  tubeSurfaceFilter->Update();

  timeCollector.Stop( "Convert to surface" );

  progress = 0.6;
  timeCollector.Start( "Write polydata file" );

  const std::string extension =
    itksys::SystemTools::GetFilenameLastExtension( outputSurfaceFile );

  if( extension == ".vtk" )
    {
    vtkNew< vtkPolyDataWriter > polyDataWriter;
    polyDataWriter->SetInputConnection( tubeSurfaceFilter->GetOutputPort() );
    polyDataWriter->SetFileName( outputSurfaceFile.c_str() );
    polyDataWriter->Write();
    }
  else if( extension == ".stl" )
    {
    // STL files only write triangles
    vtkNew< vtkTriangleFilter > triangleFilter;
    triangleFilter->SetInputConnection( tubeSurfaceFilter->GetOutputPort() );

    vtkNew< vtkSTLWriter > stlWriter;
    stlWriter->SetInputConnection( triangleFilter->GetOutputPort() );
    stlWriter->SetFileName( outputSurfaceFile.c_str() );
    stlWriter->Write();
    }
  else
    {
    tube::ErrorMessage( "Unrecognized output file extension: " + extension );
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Write polydata file" );
  progressReporter.Report( progress );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}
