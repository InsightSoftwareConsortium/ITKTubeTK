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

#ifndef __itktubeCropTubesFilter_hxx
#define __itktubeCropTubesFilter_hxx

#include "itktubeCropTubesFilter.h"
#include "tubeMacro.h"
#include "tubeTubeMathFilters.h"

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
  m_MaskImage = nullptr;
  m_BoxPositionInWorldSpace.Fill( 0 );
  m_BoxSizeInWorldSpace.Fill( 0 );
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
    pSourceTubeGroup->GetChildren( TubeGroupType::MaximumDepth, "Tube" );

  TubeGroupType* pTargetTubeGroup = this->GetOutput();

  pTargetTubeGroup->CopyInformation( pSourceTubeGroup );
  pTargetTubeGroup->Update();

  int targetTubeId=0;
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
    pCurSourceTube->Update();
    //Point List for TargetTube
    typename TubeType::TubePointListType targetPointList;
    //Get points in current source tube
    typename TubeType::TubePointListType pointList =
      pCurSourceTube->GetPoints();

    for( typename TubeType::TubePointListType::const_iterator
      pointList_it = pointList.begin();
      pointList_it != pointList.end(); ++pointList_it )
      {
      TubePointType curSourcePoint = *pointList_it;
      typename TubePointType::PointType curSourcePos =
          curSourcePoint.GetPositionInWorldSpace();
      typename TubePointType::CovariantVectorType curTubeNormal1 =
          curSourcePoint.GetNormal1InWorldSpace();
      typename TubePointType::CovariantVectorType curTubeNormal2 =
          curSourcePoint.GetNormal2InWorldSpace();
      double curRadius = curSourcePoint.GetRadiusInWorldSpace();

      std::vector<typename TubePointType::CovariantVectorType> normalList;
      normalList.push_back( curTubeNormal1 );
      if( VDimension == 3 )
        {
        normalList.push_back( curTubeNormal2 );
        }

      bool volumeMaskFlag = false;
      if( m_UseMaskImage )
        {
        typename ImageType::IndexType imageIndex;
        if( m_MaskImage->TransformPhysicalPointToIndex(
          curSourcePoint.GetPositionInWorldSpace(),
          imageIndex ) )
          {
          double val = 0;
          val = m_MaskImage->GetPixel( imageIndex );
          if( val != 0 )
            {
            volumeMaskFlag = true;
            }
          }
        }

      //Save point in target tube if it belongs to the box
      if( volumeMaskFlag || IsInsideInWorldSpace( curSourcePos, curRadius,
        m_BoxPositionInWorldSpace, m_BoxSizeInWorldSpace, normalList ) )
        {
        if( m_CropTubes )
          {
          targetPointList.push_back( curSourcePoint );
          }
        else
          {
          pCurSourceTube->SetId( targetTubeId );
          ++targetTubeId;
          pTargetTubeGroup->AddChild( pCurSourceTube );
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
          pTargetTube->Update();

          pTargetTube->ComputeTangentAndNormals();

          pTargetTube->SetId( targetTubeId );
          ++targetTubeId;
          //Save cropped tube
          pTargetTube->SetPoints( targetPointList );
          pTargetTubeGroup->AddChild( pTargetTube );

          targetPointList.clear();
          }
        }
      }
    if( targetPointList.size() > 0 )
      {
      //**** Target Tube **** :
      typename TubeType::Pointer pTargetTube = TubeType::New();

      pTargetTube->CopyInformation( pCurSourceTube );
      pTargetTube->Update();

      pTargetTube->ComputeTangentAndNormals();

      pTargetTube->SetId( targetTubeId );
      ++targetTubeId;
      //Save cropped tube
      pTargetTube->SetPoints( targetPointList );
      pTargetTubeGroup->AddChild( pTargetTube );

      targetPointList.clear();
      }
    }
  delete pSourceTubeList; 
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
CropTubesFilter< VDimension >
::IsInsideInWorldSpace( itk::Point< double, VDimension > pointPos,
  double tubeRadius,
  itk::Point< double, VDimension > boxPos,
  itk::Vector< double, VDimension > boxSize,
  std::vector<  typename itk::TubeSpatialObjectPoint
    < VDimension >::CovariantVectorType > normalList )
{
  // Return true if a slice of a tube is within the box.
  //   A slice is defined as a center point and its radius in the normal
  //     directions.
  for( unsigned int i = 0; i < normalList.size(); i++ )
    {
    for( unsigned int d = 0; d < VDimension; ++d )
      {
      if( ( pointPos[d] + tubeRadius * normalList[i][d] < boxPos[d]
          && pointPos[d] - tubeRadius * normalList[i][d] < boxPos[d] )
        || ( pointPos[d] + tubeRadius * normalList[i][d] > boxPos[d] + boxSize[d]
          && pointPos[d] - tubeRadius * normalList[i][d] > boxPos[d] + boxSize[d] ) )
        {
        return false;
        }
      }
    }
  return true;
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

#endif // End !defined( __itktubeCropTubesFilter_hxx )
