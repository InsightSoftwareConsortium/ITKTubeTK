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

#ifndef __itktubeCropTubesFilter_hxx
#define __itktubeCropTubesFilter_hxx

#include "itktubeCropTubesFilter.h"
#include "tubeMacro.h"
#include "tubeTubeMath.h"

namespace itk
{
namespace tube
{

//----------------------------------------------------------------------------
template< unsigned int VDimension >
CropTubesFilter< VDimension >
::CropTubesFilter( void )
{
  m_CropTubes = false;
  m_UseMaskImage = false;
  m_MaskImage = NULL;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
CropTubesFilter< VDimension >
::~CropTubesFilter( void )
{
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
CropTubesFilter< VDimension >
::GenerateData( void )
{
  const TubeGroupType* pSourceTubeGroup = this->GetInput();
  typename TubeGroupType::ChildrenListPointer pSourceTubeList =
    pSourceTubeGroup->GetChildren();

  TubeGroupType* pTargetTubeGroup = this->GetOutput();

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

  int targetTubeId=0;
  typename TubePointType::PointType worldBoxposition;
  typename TubePointType::VectorType worldBoxSize;

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
      return;
      }
    //Compute Tangent and Normals
    pCurSourceTube->ComputeTangentAndNormals();
    pCurSourceTube->ComputeObjectToWorldTransform();
    //Point List for TargetTube
    typename TubeType::PointListType targetPointList;
    //Get points in current source tube
    typename TubeType::PointListType pointList =
      pCurSourceTube->GetPoints();

    //Get Index to World Transformation
    typename TubeType::TransformType * pTubeObjectPhysTransform =
      pCurSourceTube->GetObjectToWorldTransform();
    typename TubeType::TransformType * pTubeIndexPhysTransform =
      pCurSourceTube->GetIndexToWorldTransform();

    worldBoxposition =
        pTubeObjectPhysTransform->TransformPoint( m_BoxPosition );
    worldBoxSize =
        pTubeObjectPhysTransform->TransformVector( m_BoxSize );

    for( typename TubeType::PointListType::const_iterator
      pointList_it = pointList.begin();
      pointList_it != pointList.end(); ++pointList_it )
      {
      TubePointType curSourcePoint = *pointList_it;
      //Transform parameters in physical space
      typename TubePointType::PointType curSourcePos =
        pTubeObjectPhysTransform->TransformPoint(
          curSourcePoint.GetPosition() );
      typename TubePointType::CovariantVectorType curTubeNormal1 =
        pTubeObjectPhysTransform->TransformCovariantVector(
          curSourcePoint.GetNormal1() );
      typename TubePointType::CovariantVectorType curTubeNormal2 =
        pTubeObjectPhysTransform->TransformCovariantVector(
          curSourcePoint.GetNormal2() );
      VectorType curRadius;
      for ( unsigned int i = 0; i < VDimension; i++ )
        {
        curRadius[i] = curSourcePoint.GetRadius();
        }
      typename TubePointType::VectorType curRadiusVector =
        pTubeObjectPhysTransform->TransformVector( curRadius );
      //Save Normals in a vector to pass it as an argument for IsIside()
      std::vector<typename TubePointType::CovariantVectorType> normalList;
      normalList.push_back( curTubeNormal1 );
      if( VDimension == 3 )
        {
        normalList.push_back( curTubeNormal2 );
        }
      bool volumeMaskFlag = false;
      typename TubePointType::PointType curSourcePosIndexSpace =
        pTubeIndexPhysTransform->TransformPoint(
        curSourcePoint.GetPosition() );
      if( m_UseMaskImage )
        {
        typename ImageType::IndexType imageIndex;
        if( m_MaskImage->TransformPhysicalPointToIndex( curSourcePosIndexSpace,
          imageIndex ) )
          {
          double val = 0;
          val = m_MaskImage->GetPixel( imageIndex );
          if ( val != 0 )
            {
            volumeMaskFlag = true;
            }
          }
        }
      //Save point in target tube if it belongs to the box
      if( volumeMaskFlag || IsInside( curSourcePos, curRadiusVector[0],
        worldBoxposition, worldBoxSize, normalList ) )
        {
        if( m_CropTubes )
          {
          targetPointList.push_back( curSourcePoint );
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
        if( targetPointList.size() > 0 )
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
          //Save cropped tube
          pTargetTube->SetPoints( targetPointList );
          pTargetTubeGroup->AddSpatialObject( pTargetTube );

          targetPointList.clear();
          }
        }
      }
    if( targetPointList.size() > 0 )
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
      //Save cropped tube
      pTargetTube->SetPoints( targetPointList );
      pTargetTubeGroup->AddSpatialObject( pTargetTube );

      targetPointList.clear();
      }
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
CropTubesFilter< VDimension >
::IsInside( itk::Point< double, VDimension > pointPos,
  double tubeRadius,
  itk::Point< double, VDimension > boxPos,
  itk::Vector< double, VDimension > boxSize,
  std::vector<  typename itk::VesselTubeSpatialObjectPoint
    < VDimension >::CovariantVectorType > normalList )
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
    switch( VDimension )
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

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
CropTubesFilter< VDimension >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

}

} // End namespace tube
} // End namespace itk

#endif // End !defined(__itktubeCropTubesFilter_hxx)
