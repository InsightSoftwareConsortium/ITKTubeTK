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
#include <iostream>
#include <sstream>

#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "metaScene.h"

#include "itkGroupSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkTimeProbesCollectorBase.h"

#include "ClipTubesCLP.h"

#include <vtkCubeSource.h>
#include <vtkNew.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

template< unsigned int DimensionT >
bool IsInside( itk::Point< double, DimensionT > pointPos, double tubeRadius,
  itk::Point< double, DimensionT > boxPos,
  itk::Vector< double, DimensionT > boxSize,
  std::vector<  typename itk::VesselTubeSpatialObjectPoint
    < DimensionT >::CovariantVectorType > normalList )
{
  // Return a boolean indicating if any slice of the tube
  // is included in the box.
  // A slice is considered as a point and an associated radius
  for( unsigned int i = 0; i < normalList.size(); i++ )
    {
    bool hasXInside = false;
    bool hasYInside = false;
    bool hasZInside = false;
    if( pointPos[0] + tubeRadius * normalList[i][0] >= boxPos[0] &&
      pointPos[0] - tubeRadius * normalList[i][0] <= boxPos[0] + boxSize[0] )
      {
      hasXInside = true;
      }
    if( pointPos[1] + tubeRadius * normalList[i][1] >= boxPos[1] &&
      pointPos[1] - tubeRadius * normalList[i][1] <= boxPos[1] + boxSize[1] )
      {
      hasYInside = true;
      }
    switch( DimensionT )
      {
      case 2:
        {
        hasZInside = true;
        break;
        }
      case 3:
        {
        if( pointPos[2] + tubeRadius  * normalList[i][2] >= boxPos[2] &&
          pointPos[2] - tubeRadius * normalList[i][2] <=
            boxPos[2] + boxSize[2] )
          {
          hasZInside = true;
          }
        break;
        }
      default:
        {
        tubeErrorMacro(
          << "Error: Only 2D and 3D data is currently supported." );
        return EXIT_FAILURE;
        }
      }
    if( hasXInside && hasYInside && hasZInside )
      {
      return true;
      }
    }
  return false;
}

template< unsigned int DimensionT >
void WriteBox( std::vector< double > boxPos, std::vector< double > boxSize,
  std::vector< double > spacing, std::string boxFileName )
{
  if( DimensionT == 3 )
    {
    vtkSmartPointer<vtkCubeSource> cubeSource =
      vtkSmartPointer<vtkCubeSource>::New();
    cubeSource->SetBounds(
      -spacing[0] * ( boxPos[0] + boxSize[0] ), -spacing[0] * boxPos[0],
      -spacing[1] * ( boxPos[1] + boxSize[1] ), -spacing[1] * boxPos[1],
      spacing[2] * boxPos[2], spacing[2] * ( boxPos[2] + boxSize[2] ) );
    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName( boxFileName.c_str() );
    writer->SetInputConnection( cubeSource->GetOutputPort() );
    writer->Write();
    }
  else
    {
    vtkSmartPointer< vtkCubeSource > cubeSource =
      vtkSmartPointer< vtkCubeSource >::New();
    cubeSource->SetBounds(
      -spacing[0] * ( boxPos[0] + boxSize[0] ), -spacing[0] * boxPos[0],
      -spacing[1] * ( boxPos[1] + boxSize[1] ), -spacing[1] * boxPos[1],
      0, 0 );
    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName( boxFileName.c_str() );
    writer->SetInputConnection( cubeSource->GetOutputPort() );
    writer->Write();
    }
}

template< unsigned int DimensionT >
int DoIt (int argc, char * argv[])
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if( DimensionT != 2 && DimensionT != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< DimensionT >          TubesReaderType;
  typedef itk::GroupSpatialObject< DimensionT >           TubeGroupType;
  typedef itk::VesselTubeSpatialObject< DimensionT >      TubeType;
  typedef itk::VesselTubeSpatialObjectPoint< DimensionT > TubePointType;
  typedef double                                          PixelType;
  typedef itk::Image< PixelType, DimensionT >             ImageType;
  typedef itk::ImageFileReader< ImageType >               ImageReaderType;
  typedef itk::Vector< double, DimensionT >               VectorType;
  typedef itk::Point< double, DimensionT >                PointType;

  timeCollector.Start( "Loading Input TRE File" );

  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();

  try
    {
    tubeFileReader->SetFileName( inputTREFile.c_str() );
    tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE File: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  //Cast XML vector parameters
  PointType boxPositionVector;
  VectorType boxSizeVector;
  if ( !boxCorner.empty() )
    {
    for ( unsigned int i = 0; i < DimensionT; i++ )
      {
      boxPositionVector[i] = boxCorner[i];
      boxSizeVector[i] = boxSize[i];
      }
    }
  else
    {
    for ( unsigned int i = 0; i < DimensionT; i++ )
      {
      boxPositionVector[i] = -1;
      boxSizeVector[i] = -1;
      }
    }

  typename TubeGroupType::Pointer pSourceTubeGroup =
    tubeFileReader->GetGroup();
  typename TubeGroupType::ChildrenListPointer pSourceTubeList =
    pSourceTubeGroup->GetChildren();

  timeCollector.Stop( "Loading Input TRE File" );

  //loading Volume mask if its there
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  ImageType::Pointer image;
  if ( !volumeMask.empty() )
    {
    tube::InfoMessage( "Reading volume mask..." );
    timeCollector.Start( "Load Volume Mask" );
    imReader->SetFileName( volumeMask.c_str() );
    try
      {
      imReader->Update();
      image = imReader->GetOutput();
      }
    catch ( itk::ExceptionObject & err )
      {
      tube::FmtErrorMessage( "Cannot read volume mask file: %s",
        err.what() );
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Load Volume Mask" );
    }

  // Compute clipping
  tubeStandardOutputMacro( << "\n>> Finding Tubes for Clipping" );

  timeCollector.Start( "Selecting Tubes" );
  //Target Group to save desired tubes
  typename TubeGroupType::Pointer pTargetTubeGroup = TubeGroupType::New();

  pTargetTubeGroup->CopyInformation( pSourceTubeGroup );
  // TODO: make CopyInformation of itk::SpatialObject do this
  pTargetTubeGroup->GetObjectToParentTransform()->SetScale(
    pSourceTubeGroup->GetObjectToParentTransform()->GetScale() );
  pTargetTubeGroup->GetObjectToParentTransform()->SetOffset(
    pSourceTubeGroup->GetObjectToParentTransform()->GetOffset() );
  pTargetTubeGroup->GetObjectToParentTransform()->SetMatrix(
    pSourceTubeGroup->GetObjectToParentTransform()->GetMatrix() );
  pTargetTubeGroup->SetSpacing( pSourceTubeGroup->GetSpacing() );
  pTargetTubeGroup->ComputeObjectToWorldTransform();
  typename TubePointType::PointType curSourceOffset =
    pSourceTubeGroup->GetObjectToParentTransform()->GetOffset();
  int targetTubeId=0;
  std::vector< double > spacing( DimensionT, 0 );

  for( typename TubeGroupType::ChildrenListType::iterator
    tubeList_it = pSourceTubeList->begin();
    tubeList_it != pSourceTubeList->end(); ++tubeList_it )
    {
    //**** Source Tube **** :
    typename TubeType::Pointer pCurSourceTube =
      dynamic_cast< TubeType* >( tubeList_it->GetPointer() );
    //dynamic_cast verification
    if( !pCurSourceTube )
      {
      return EXIT_FAILURE;
      }
    //Compute Tangent and Normals
    pCurSourceTube->ComputeTangentAndNormals();
    //pCurSourceTube->ComputeObjectToWorldTransform();//BUG
    //Point List for TargetTube
    typename TubeType::PointListType TargetPointList;
    //Get points in current source tube
    typename TubeType::PointListType pointList =
      pCurSourceTube->GetPoints();

    //Get Index to World Transformation
    typename TubeType::TransformType * pTubeIndexPhysTransform =
      pCurSourceTube->GetIndexToWorldTransform();

    for( typename TubeType::PointListType::const_iterator
      pointList_it = pointList.begin();
      pointList_it != pointList.end(); ++pointList_it )
      {
      TubePointType curSourcePoint = *pointList_it;
      //Transform parameters in physical space
      typename TubePointType::PointType curSourcePos =
        pTubeIndexPhysTransform->TransformPoint(
          curSourcePoint.GetPosition() );
      typename TubePointType::PointType worldBoxposition =
        pTubeIndexPhysTransform->TransformPoint( boxPositionVector );
      typename TubePointType::VectorType worldBoxSize =
        pTubeIndexPhysTransform->TransformVector( boxSizeVector );
      typename TubePointType::CovariantVectorType curTubeNormal1 =
        pTubeIndexPhysTransform->TransformCovariantVector(
          curSourcePoint.GetNormal1() );
      typename TubePointType::CovariantVectorType curTubeNormal2 =
        pTubeIndexPhysTransform->TransformCovariantVector(
          curSourcePoint.GetNormal2() );
      VectorType curRadius;
      for ( unsigned int i = 0; i < DimensionT; i++ )
        {
        curRadius[i] = curSourcePoint.GetRadius();
        }
      typename TubePointType::VectorType curRadiusVector =
        pTubeIndexPhysTransform->TransformVector( curRadius );
      //Save Normals in a vector to pass it as an argument for IsIside()
      std::vector<typename TubePointType::CovariantVectorType> normalList;
      normalList.push_back( curTubeNormal1 );
      if( DimensionT == 3 )
        {
        normalList.push_back( curTubeNormal2 );
        }

      bool volumeMaskFlag = false;
      if ( !volumeMask.empty() )
        {
        typename ImageType::IndexType imageIndex;
        if ( image->TransformPhysicalPointToIndex( curSourcePos, imageIndex ) )
          {
          double val = 0;
          val = image->GetPixel( imageIndex );
          if ( val != 0 )
            {
            volumeMaskFlag = true;
            }
          }
        }
      //Save point in target tube if it belongs to the box
      if ( volumeMaskFlag || IsInside( curSourcePos, curRadiusVector[0],
        worldBoxposition, worldBoxSize, normalList ) )
        {
        if( ClipTubes )
          {
          TargetPointList.push_back( curSourcePoint );
          }
        else
          {
          pCurSourceTube->SetId( targetTubeId );
          ++targetTubeId;
          pTargetTubeGroup->AddSpatialObject( pCurSourceTube );
          break;
          }
        }
      else
        {
        if( TargetPointList.size() > 0 )
          {
          //**** Target Tube **** :
          typename TubeType::Pointer pTargetTube = TubeType::New();

          pTargetTube->CopyInformation( pCurSourceTube );

          // TODO: make CopyInformation of itk::SpatialObject do this
          pTargetTube->GetObjectToParentTransform()->SetScale(
            pCurSourceTube->GetObjectToParentTransform()->GetScale() );
          pTargetTube->GetObjectToParentTransform()->SetOffset(
            pCurSourceTube->GetObjectToParentTransform()->GetOffset() );
          pTargetTube->GetObjectToParentTransform()->SetMatrix(
            pCurSourceTube->GetObjectToParentTransform()->GetMatrix() );
          pTargetTube->SetSpacing( pCurSourceTube->GetSpacing() );
          pTargetTube->ComputeObjectToWorldTransform();

          pTargetTube->ComputeTangentAndNormals();

          pTargetTube->SetId( targetTubeId );
          ++targetTubeId;
          //Save clipped tube
          pTargetTube->SetPoints( TargetPointList );
          pTargetTubeGroup->AddSpatialObject( pTargetTube );

          TargetPointList.clear();
          }
        }
      }
    if( TargetPointList.size() > 0 )
      {
      //**** Target Tube **** :
      typename TubeType::Pointer pTargetTube = TubeType::New();

      pTargetTube->CopyInformation( pCurSourceTube );

      // TODO: make CopyInformation of itk::SpatialObject do this
      pTargetTube->GetObjectToParentTransform()->SetScale(
        pCurSourceTube->GetObjectToParentTransform()->GetScale() );
      pTargetTube->GetObjectToParentTransform()->SetOffset(
        pCurSourceTube->GetObjectToParentTransform()->GetOffset() );
      pTargetTube->GetObjectToParentTransform()->SetMatrix(
        pCurSourceTube->GetObjectToParentTransform()->GetMatrix() );
      pTargetTube->SetSpacing( pCurSourceTube->GetSpacing() );
      pTargetTube->ComputeObjectToWorldTransform();

      pTargetTube->ComputeTangentAndNormals();

      pTargetTube->SetId( targetTubeId );
      ++targetTubeId;
      //Save clipped tube
      pTargetTube->SetPoints( TargetPointList );
      pTargetTubeGroup->AddSpatialObject( pTargetTube );

      TargetPointList.clear();
      }
    for( unsigned int d=0; d<DimensionT; ++d )
      {
      spacing[d] = pCurSourceTube->GetSpacing()[d];
      }
    }
  if( !boxCorner.empty() && !outputBoxFile.empty() )
    {
    for ( unsigned int i = 0; i < DimensionT; i++ )
      {
      boxCorner[i] = boxCorner[i] + curSourceOffset[i];
      }
    WriteBox<DimensionT>( boxCorner, boxSize, spacing, outputBoxFile );
    }


  timeCollector.Stop( "Selecting Tubes" );

  // Write output TRE file
  tubeStandardOutputMacro(
    << "\n>> Writing TRE file" );

  timeCollector.Start( "Writing output TRE file" );

  typedef itk::SpatialObjectWriter< DimensionT > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();

  try
    {
    tubeWriter->SetFileName( outputTREFile.c_str() );
    tubeWriter->SetInput( pTargetTubeGroup );
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing output TRE file" );
  timeCollector.Report();
  return EXIT_SUCCESS;
}


// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;

  if( (boxCorner.empty() || boxSize.empty()) && volumeMask.empty() )
    {
    tube::ErrorMessage(
      "Error: Either both longflags --boxCorner and --boxSize "
      "or the flag --volumeMask is required." );
    return EXIT_FAILURE;
    }

  MetaScene *mScene = new MetaScene;
  mScene->Read( inputTREFile.c_str() );

  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input TRE file has no spatial objects" );
    return EXIT_SUCCESS;
    }

  switch( mScene->GetObjectList()->front()->NDims() )
    {
    case 2:
      {
      return DoIt<2>( argc, argv );
      }

    case 3:
      {
      return DoIt<3>( argc, argv );
      }

    default:
      {
      tubeErrorMacro(
        << "Error: Only 2D and 3D data is currently supported." );
      return EXIT_FAILURE;
      }
    }
}
