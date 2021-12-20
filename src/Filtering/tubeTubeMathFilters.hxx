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
#ifndef __tubeTubeMathFilters_hxx
#define __tubeTubeMathFilters_hxx

#include "itkNumericTraits.h"
#include "itkImageRegionIterator.h"
#include "itkDanielssonDistanceMapImageFilter.h"

#include "tubeImageMathFilters.h"
#include "itktubeTubeSpatialObjectToImageFilter.h"

namespace tube
{

template< unsigned int DimensionT, class ImagePixelT >
TubeMathFilters< DimensionT, ImagePixelT >::
TubeMathFilters()
{
  m_InputTubeGroup = nullptr;
  m_InputTube = nullptr;

  m_CurrentTubeId = -1;

  m_TubePointIdImage = nullptr;
  m_TubeRadiusImage = nullptr;
  m_TubeDistanceImage = nullptr;
  m_TubeDirectionImage = nullptr;
}

template< unsigned int DimensionT, class ImagePixelT >
TubeMathFilters< DimensionT, ImagePixelT >::
~TubeMathFilters()
{
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetInputTubeGroup( TubeGroupType * inputTubeGroup )
{
  m_InputTubeGroup = inputTubeGroup;
  m_InputTube = nullptr;

  m_CurrentTubeId = -1;

  m_TubePointIdImage = nullptr;
  m_TubeRadiusImage = nullptr;
  m_TubeDistanceImage = nullptr;
  m_TubeDirectionImage = nullptr;
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetInputTube( TubeType * inputTube )
{
  m_InputTube = inputTube;
  m_InputTubeGroup = TubeGroupType::New();
  m_InputTubeGroup->AddChild( m_InputTube );

  m_CurrentTubeId = -1;

  m_TubePointIdImage = nullptr;
  m_TubeRadiusImage = nullptr;
  m_TubeDistanceImage = nullptr;
  m_TubeDirectionImage = nullptr;
}

template< unsigned int DimensionT, class ImagePixelT >
typename itk::GroupSpatialObject<DimensionT>::Pointer &
TubeMathFilters< DimensionT, ImagePixelT >::
GetOutputTubeGroup( void )
{
  return m_InputTubeGroup;
}

template< unsigned int DimensionT, class ImagePixelT >
typename itk::TubeSpatialObject<DimensionT>::Pointer &
TubeMathFilters< DimensionT, ImagePixelT >::
GetOutputTube( void )
{
  return m_InputTube;
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetCurrentTubeId( int currentTubeId )
{
  m_CurrentTubeId = currentTubeId;
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetUseAllTubes( )
{
  m_CurrentTubeId = -1;
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetPointValues( std::string propertyId, double val, double blend )
{
  typename TubeType::ChildrenListType::iterator tubeIterator;
  typename TubeType::ChildrenListPointer inputTubeList =
    m_InputTubeGroup->GetChildren( m_InputTubeGroup->GetMaximumDepth(),
    "Tube" );
  for( tubeIterator = inputTubeList->begin(); tubeIterator !=
    inputTubeList->end(); tubeIterator++ )
    {
    typename TubeType::Pointer inputTube = ( ( TubeType * )( tubeIterator->
      GetPointer() ) );

    if( m_CurrentTubeId == -1 || inputTube->GetId() == m_CurrentTubeId )
      {
      inputTube->Update();

      unsigned int pointListSize = inputTube->GetNumberOfPoints();
      for( unsigned int pointNum = 0; pointNum < pointListSize; ++pointNum )
        {
        TubePointType * currentPoint = static_cast< TubePointType * >(
          inputTube->GetPoint( pointNum ) );

        if( propertyId == "Ridgeness" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetRidgeness();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetRidgeness( val );
          }
        else if( propertyId == "Medialness" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetMedialness();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetMedialness( val );
          }
        else if( propertyId == "Branchness" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetBranchness();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetBranchness( val );
          }
        else if( propertyId == "Radius" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetRadiusInObjectSpace();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetRadiusInObjectSpace( val );
          }
        else
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetTagScalarValue( propertyId );
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetTagScalarValue( propertyId, val );
          }
        }
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetPointValuesFromImage(
  const itk::Image< ImagePixelT, DimensionT> * inputImage,
  std::string propertyId, double blend )
{
  typename TubeType::ChildrenListType::iterator tubeIterator;
  typename TubeType::ChildrenListPointer inputTubeList =
    m_InputTubeGroup->GetChildren( m_InputTubeGroup->GetMaximumDepth(),
    "Tube" );
  for( tubeIterator = inputTubeList->begin(); tubeIterator !=
    inputTubeList->end(); tubeIterator++ )
    {
    typename TubeType::Pointer inputTube = ( ( TubeType * )( tubeIterator->
      GetPointer() ) );

    if( m_CurrentTubeId == -1 || inputTube->GetId() == m_CurrentTubeId )
      {
      inputTube->Update();

      unsigned int pointListSize = inputTube->GetNumberOfPoints();
      for( unsigned int pointNum = 0; pointNum < pointListSize; ++pointNum )
        {
        TubePointType * currentPoint = static_cast< TubePointType * >(
          inputTube->GetPoint( pointNum ) );
        typename TubeType::PointType pointIndex;
        pointIndex = currentPoint->GetPositionInObjectSpace();

        typename TubeType::PointType pointWorld;
        pointWorld = inputTube->GetObjectToWorldTransform()->TransformPoint(
          pointIndex );

        typename ImageType::IndexType imageIndex;
        double val = 0;
        if( inputImage->TransformPhysicalPointToIndex( pointWorld,
          imageIndex ) )
          {
          val = inputImage->GetPixel( imageIndex );
          }
        if( propertyId == "Ridgeness" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetRidgeness();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetRidgeness( val );
          }
        else if( propertyId == "Medialness" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetMedialness();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetMedialness( val );
          }
        else if( propertyId == "Branchness" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetBranchness();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetBranchness( val );
          }
        else if( propertyId == "Radius" )
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetRadiusInObjectSpace();
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetRadiusInObjectSpace( val );
          }
        else
          {
          if( blend != 1 )
            {
            double val2 = currentPoint->GetTagScalarValue( propertyId );
            val = val * blend + (1 - blend) * val2;
            }
          currentPoint->SetTagScalarValue( propertyId, val );
          }
        }
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}

template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetPointValuesFromImageMean( 
  const itk::Image< ImagePixelT, DimensionT> * inputImage,
  std::string propertyId )
{
  typename TubeType::ChildrenListType::iterator tubeIterator;

  typename TubeType::ChildrenListPointer inputTubeList =
    m_InputTubeGroup->GetChildren( m_InputTubeGroup->GetMaximumDepth(),
    "Tube" );
  for( tubeIterator = inputTubeList->begin(); tubeIterator !=
    inputTubeList->end(); tubeIterator++ )
    {
    typename TubeType::Pointer inputTube = ( ( TubeType * )( tubeIterator->
      GetPointer() ) );

    if( m_CurrentTubeId == -1 || inputTube->GetId() == m_CurrentTubeId )
      {
      inputTube->Update();

      double valAvg = 0;
      unsigned int valCount = 0;
      unsigned int pointListSize = inputTube->GetNumberOfPoints();
      for( unsigned int pointNum = 0; pointNum < pointListSize; ++pointNum )
        {
        TubePointType * currentPoint = static_cast< TubePointType * >(
          inputTube->GetPoint( pointNum ) );
        typename TubeType::PointType pointIndex;
        pointIndex = currentPoint->GetPositionInObjectSpace();

        typename TubeType::PointType pointWorld;
        pointWorld = inputTube->GetObjectToWorldTransform()->TransformPoint(
          pointIndex );

        typename ImageType::IndexType imageIndex;
        if( inputImage->TransformPhysicalPointToIndex( pointWorld,
          imageIndex ) )
          {
          valAvg += inputImage->GetPixel( imageIndex );
          ++valCount;
          }
        }
      valAvg /= valCount;
      for( unsigned int pointNum = 0; pointNum < pointListSize; ++pointNum )
        {
        TubePointType * currentPoint = static_cast< TubePointType * >(
          inputTube->GetPoint( pointNum ) );
        typename TubeType::PointType pointIndex;
        pointIndex = currentPoint->GetPositionInObjectSpace();

        typename TubeType::PointType pointWorld;
        pointWorld = inputTube->GetObjectToWorldTransform()->TransformPoint(
          pointIndex );

        if( propertyId == "Ridgeness" )
          {
          currentPoint->SetRidgeness( valAvg );
          }
        else if( propertyId == "Medialness" )
          {
          currentPoint->SetMedialness( valAvg );
          }
        else if( propertyId == "Branchness" )
          {
          currentPoint->SetBranchness( valAvg );
          }
        else if( propertyId == "Radius" )
          {
          currentPoint->SetRadiusInObjectSpace( valAvg );
          }
        else
          {
          currentPoint->SetTagScalarValue( propertyId, valAvg );
          }
        }
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}

//------------------------------------------------------------------------
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
ComputeTubeRegions( const ImageType * referenceImage )
{
  this->RenumberTubes();
  this->RenumberPoints();

  std::cout << "Computing Tube Regions" << std::endl;

  std::cout << "   Rendering Tube Centerlines" << std::endl;
  typedef itk::tube::TubeSpatialObjectToImageFilter< DimensionT,
    FloatImageType > TubeToImageFilterType;
  typename TubeToImageFilterType::Pointer tubeToImageFilter =
    TubeToImageFilterType::New();
  tubeToImageFilter->SetInput( m_InputTubeGroup );
  tubeToImageFilter->SetColorByPointId(true);
  tubeToImageFilter->SetUseRadius(false);
  tubeToImageFilter->SetBuildRadiusImage(true);
  tubeToImageFilter->SetBuildTangentImage(false);
  tubeToImageFilter->SetReferenceImage<const ImageType>( referenceImage );
  tubeToImageFilter->Update();

  m_TubePointIdImage = tubeToImageFilter->GetOutput();
  m_TubeRadiusImage  = tubeToImageFilter->GetRadiusImage();

  std::cout << "   Computing Tube Distance Maps" << std::endl;
  typedef itk::DanielssonDistanceMapImageFilter< FloatImageType, FloatImageType>
    DanielssonFilterType;
  typename DanielssonFilterType::Pointer distF = DanielssonFilterType::New();
  distF->SetInput( m_TubePointIdImage );
  distF->SetUseImageSpacing( true );
  distF->SetInputIsBinary( true );
  distF->Update();

  m_TubeDistanceImage  = distF->GetDistanceMap();
  m_TubeDirectionImage = distF->GetVectorDistanceMap();
  std::cout << "   Done." << std::endl;
}

//------------------------------------------------------------------------
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetPointValuesFromTubeRegions(
  const itk::Image< ImagePixelT, DimensionT> * inputImage,
  const std::string & propertyId,
  double minRFactor, double maxRFactor )
{
  std::cout << "Setting Tube Point Values from Regions." << std::endl;
  this->SetPointValues( propertyId, 0 );

  typedef itk::ImageRegionIterator<VectorImageType>   VectorImageIteratorType;
  typedef itk::ImageRegionIterator<FloatImageType>    FloatImageIteratorType;

  FloatImageIteratorType itDistImage( m_TubeDistanceImage,
      m_TubeDistanceImage->GetLargestPossibleRegion() );
  itDistImage.GoToBegin();

  VectorImageIteratorType itDirImage( m_TubeDirectionImage,
      m_TubeDirectionImage->GetLargestPossibleRegion() );
  itDirImage.GoToBegin();

  std::cout << "   Allocating accumulator images." << std::endl;
  typename FloatImageType::Pointer imVal = FloatImageType::New();
  imVal->CopyInformation( m_TubeDistanceImage );
  imVal->SetRegions( m_TubeDistanceImage->GetLargestPossibleRegion() );
  imVal->Allocate();
  imVal->FillBuffer(0);

  typename FloatImageType::Pointer imCount = FloatImageType::New();
  imCount->CopyInformation( m_TubeDistanceImage );
  imCount->SetRegions( m_TubeDistanceImage->GetLargestPossibleRegion() );
  imCount->Allocate();
  imCount->FillBuffer(0);

  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  unsigned int lastIndex = 1;
  unsigned int maxIndex =
    imVal->GetLargestPossibleRegion().GetSize()[DimensionT-1];
  std::cout << "   Parsing regions." << std::endl;
  while( !itDirImage.IsAtEnd() )
    {
    double dist = itDistImage.Value();

    typename VectorImageType::IndexType index = itDirImage.GetIndex();
    VectorPixelType v = itDirImage.Value();

    typename VectorImageType::IndexType tubeIndex;
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
      tubeIndex[i] = index[i] + v[i];
      }

    double val = inputImage->GetPixel( index );

    if( val != 0 )
      {
      std::cout << tubeIndex << " = " << val << std::endl;
      }
 
    if( imVal->GetLargestPossibleRegion().IsInside(tubeIndex) )
      {
      imVal->SetPixel(tubeIndex, imVal->GetPixel(tubeIndex)+val);
      imCount->SetPixel(tubeIndex, imCount->GetPixel(tubeIndex)+1);
      }

    ++itDistImage;
    ++itDirImage;
    }

  FloatImageIteratorType itId( m_TubePointIdImage,
    m_TubePointIdImage->GetLargestPossibleRegion() );
  FloatImageIteratorType itVal( imVal, imVal->GetLargestPossibleRegion() );
  FloatImageIteratorType itCount( imCount, imCount->GetLargestPossibleRegion() );

  itId.GoToBegin();
  itVal.GoToBegin();
  itCount.GoToBegin();

  lastIndex = 1;
  std::cout << "   Storing results." << std::endl;
  while( !itId.IsAtEnd() )
    {
    if( itCount.Value() != 0 )
      {
      double pointId = itId.Value();

      double val = itVal.Value() / itCount.Value();

      typename VectorImageType::IndexType index = itVal.GetIndex();
      if( index[DimensionT-1] != lastIndex )
        {
        std::cout << "  Slice = " << index[DimensionT-1] << " of " << maxIndex
          << std::endl;
        lastIndex = index[DimensionT-1];
        }

      int tube = static_cast<int>(pointId);
      int point = 0;
      if( pointId - tube != 0 )
        {
        point = static_cast<int>(1 / (pointId-tube))-1;
        }

      if( val != 0 )
        {
        static int lcount = 0;
        if( lcount < 100 )
          {
          std::cout << tube << ", " << point << " = " << val << std::endl;
          lcount++;
          }
        }

      typename TubeType::ChildrenListType::iterator tubeIterator;
      tubeIterator = inputTubeList->begin();
      for(int tubeI=0; tubeI<tube; ++tubeI)
        {
        ++tubeIterator;
        }
      typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
        tubeIterator->GetPointer() );

      TubePointType * currentPoint = static_cast< TubePointType * >(
        curTube->GetPoint( point ) );
      if( propertyId == "Ridgeness" )
        {
        currentPoint->SetRidgeness(val);
        }
      else if( propertyId == "Medialness" )
        {
        currentPoint->SetMedialness(val);
        }
      else if( propertyId == "Branchness" )
        {
        currentPoint->SetBranchness(val);
        }
      else if( propertyId == "Radius" )
        {
        currentPoint->SetRadiusInObjectSpace(val);
        }
      else
        {
        currentPoint->SetTagScalarValue( propertyId, val );
        }
      }
    ++itId;
    ++itVal;
    ++itCount;
    }
  inputTubeList->clear();
  delete inputTubeList;
}


//------------------------------------------------------------------------
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SetPointValuesFromTubeRadius(
  const itk::Image< ImagePixelT, DimensionT> * inputImage,
  const std::string & propertyId,
  double minRFactor, double maxRFactor )
{
  this->SetPointValues( propertyId, 0 );

  typedef itk::ContinuousIndex<double, DimensionT> ContinuousIndexType;

  // Get the list of tubes
  char tubeName[] = "Tube";
  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  typename TubeGroupType::ChildrenListType::iterator tubeIterator =
    inputTubeList->begin();

  ContinuousIndexType pointI;
  typename FloatImageType::IndexType index;
  typename FloatImageType::IndexType index2;
  while( tubeIterator != inputTubeList->end() )
    {
    TubeType * tube = ( TubeType * )tubeIterator->GetPointer();

    tube->Update();

    for( unsigned int k=0; k < tube->GetNumberOfPoints(); k++ )
      {
      typedef typename TubeType::TubePointType TubePointType;
      TubePointType * tubePoint = dynamic_cast<TubePointType *>(
        tube->GetPoint( k ) );
      inputImage->TransformPhysicalPointToContinuousIndex(
        tubePoint->GetPositionInWorldSpace(),
        pointI );
      for( unsigned int i=0; i<DimensionT; i++ )
        {
        index[i] = ( long int )( pointI[i]+0.5 );
        }
      double val = 0;
      unsigned int count = 0;
      if( inputImage->GetLargestPossibleRegion().IsInside( index ) )
        {
        double radius = tubePoint->GetRadiusInWorldSpace();
        ContinuousIndexType radiusI;
        ContinuousIndexType maxRadiusI;
        for( unsigned int i = 0; i < DimensionT; i++ )
          {
          maxRadiusI[i] = (maxRFactor * radius) / inputImage->GetSpacing()[i];
          radiusI[i] = radius / inputImage->GetSpacing()[i];
          if(radiusI[i] < 0.25)
            {
            radiusI[i] = 0.25;
            }
          }

        if( DimensionT == 2 )
          {
          for( double x=-maxRadiusI[0]; x<=maxRadiusI[0]; x += 0.5 )
            {
            for( double y=-maxRadiusI[1]; y<=maxRadiusI[1]; y += 0.5 )
              {
              index2[0]=( long )( pointI[0]+x+0.5 );
              index2[1]=( long )( pointI[1]+y+0.5 );
              if( inputImage->GetLargestPossibleRegion().IsInside( index2 ) )
                {
                double dist = m_TubeDistanceImage->GetPixel(index2);
                if( dist >= minRFactor && dist <= maxRFactor )
                  {
                  val += inputImage->GetPixel(index2);
                  ++count;
                  }
                }
              }
            }
          }
        else if( DimensionT == 3 )
          {
          for( double x=-maxRadiusI[0]; x<=maxRadiusI[0]; x+=0.5 )
            {
            for( double y=-maxRadiusI[1]; y<=maxRadiusI[1]; y+=0.5 )
              {
              for( double z=-maxRadiusI[2]; z<=maxRadiusI[2]; z+=0.5 )
                {
                index2[0]=( long )( pointI[0]+x+0.5 );
                index2[1]=( long )( pointI[1]+y+0.5 );
                index2[2]=( long )( pointI[2]+z+0.5 );
                if( inputImage->GetLargestPossibleRegion().IsInside( index2 ) )
                  {
                  double dist = m_TubeDistanceImage->GetPixel(index2);
                  if( dist >= minRFactor && dist <= maxRFactor )
                    {
                    val += inputImage->GetPixel(index2);
                    ++count;
                    }
                  }
                }
              }
            }
          }
        }
      if( count > 0 )
        {
        val /= count;
        if( propertyId == "Ridgeness" )
          {
          tubePoint->SetRidgeness(val);
          }
        else if( propertyId == "Medialness" )
          {
          tubePoint->SetMedialness(val);
          }
        else if( propertyId == "Branchness" )
          {
          tubePoint->SetBranchness(val);
          }
        else if( propertyId == "Radius" )
          {
          tubePoint->SetRadiusInObjectSpace(val);
          }
        else
          {
          tubePoint->SetTagScalarValue( propertyId, val );
          }
        }
      }
    ++tubeIterator;
    }
  inputTubeList->clear();
  delete inputTubeList;
} 


//------------------------------------------------------------------------
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
FillGapToParent( double stepSize )
{
  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
      itCurTube->GetPointer() );
    if( m_CurrentTubeId == -1 || curTube->GetId() == m_CurrentTubeId )
      {
      TubeIdType curParentTubeId = curTube->GetParentId();
      TubePointType* parentNearestPoint = NULL;
  
      if( curTube->GetRoot() == false && curParentTubeId != curTube->GetId() )
        {
        //find parent target tube
        for( typename TubeGroupType::ChildrenListType::iterator itTubes =
          inputTubeList->begin(); itTubes != inputTubeList->end(); ++itTubes )
          {
          typename TubeType::Pointer tube = dynamic_cast< TubeType * >(
            itTubes->GetPointer() );
          if( tube->GetId() == curParentTubeId )
            {
            double minDistance = itk::NumericTraits<double>::max();
            int flag =-1;
            for( unsigned int index = 0; index < tube->GetNumberOfPoints();
              ++index )
              {
              TubePointType* tubePoint = dynamic_cast< TubePointType* >(
                tube->GetPoint( index ) );
              PositionType tubePointPosition =
                tubePoint->GetPositionInObjectSpace();
              double distance = tubePointPosition.EuclideanDistanceTo(
                curTube->GetPoint( 0 )->GetPositionInObjectSpace() );
              if( minDistance > distance )
                {
                minDistance = distance;
                parentNearestPoint = tubePoint;
                flag = 1;
                }
              distance = tubePointPosition.EuclideanDistanceTo(
                curTube->GetPoint( curTube->GetNumberOfPoints() - 1 )
                ->GetPositionInObjectSpace() );
              if( minDistance > distance )
                {
                minDistance = distance;
                parentNearestPoint = tubePoint;
                flag = 2;
                }
              }
  
            TubePointListType newTubePoints;
            if( flag == 1 )
              {
              TubePointType* childTubeStartPoint = dynamic_cast<
                TubePointType* >( curTube->GetPoint( 0 ) );
              InterpolatePath( parentNearestPoint,
                childTubeStartPoint, stepSize, newTubePoints );
              TubePointListType targetTubePoints = curTube->GetPoints();
              curTube->Clear();
              for( unsigned int index = 0; index < newTubePoints.size();
                ++index )
                {
                curTube->GetPoints().push_back( newTubePoints[ index ] );
                }
              for( unsigned int i = 0; i < targetTubePoints.size(); ++i )
                {
                curTube->GetPoints().push_back( targetTubePoints[ i ] );
                }
              }
            if( flag == 2 )
              {
              TubePointType* childTubeEndPoint =
                dynamic_cast< TubePointType* >
                ( curTube->GetPoint( curTube->GetNumberOfPoints() - 1 ) );
              InterpolatePath( parentNearestPoint,
                childTubeEndPoint, stepSize, newTubePoints );
              for( int index = newTubePoints.size() - 1; index >= 0; index-- )
                {
                curTube->GetPoints().push_back( newTubePoints[ index ] );
                }
              }
            break;
            }
          }
        }
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}

//------------------------------------------------------------------------
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
InterpolatePath(
  typename TubeType::TubePointType * parentNearestTubePoint,
  typename TubeType::TubePointType * childEndTubePoint,
  float stepSize,
  typename TubeType::TubePointListType & newTubePoints )
{
  double projLen = 0;
  typename TubePointType::PointType midPoint;
  for( unsigned int i=0; i<DimensionT; ++i )
    {
    double tf = (parentNearestTubePoint->GetPositionInObjectSpace()[i]
      - childEndTubePoint->GetPositionInObjectSpace()[i]);
    midPoint[i] = tf/2;
    tf = tf * childEndTubePoint->GetTangentInObjectSpace()[i];
    projLen += tf * tf;
    }
  double midProjLen = std::sqrt( projLen ) / 2;
  for( unsigned int i=0; i<DimensionT; ++i )
    {
    midPoint[i] = midPoint[i] + childEndTubePoint->GetTangentInObjectSpace()[i] * midProjLen;
    midPoint[i] /= 2;
    }
  TubePointType midTubePoint;
  midTubePoint.SetPositionInObjectSpace( midPoint );
  if( stepSize > 0 && midProjLen > stepSize )
    {
    InterpolatePath( &midTubePoint, childEndTubePoint, stepSize,
                     newTubePoints );
    }
  newTubePoints.push_back( midTubePoint );
  if( stepSize > 0 && midProjLen > stepSize )
    {
    InterpolatePath( parentNearestTubePoint, &midTubePoint, stepSize,
                     newTubePoints );
    }
  newTubePoints.push_back( *parentNearestTubePoint );
  return;
}


/**
 * Smooth a tube */
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SmoothTube( double h, SmoothTubeFunctionEnum smoothFunction )
{
  if( h <= 0 )
    {
    return;
    }

  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
      itCurTube->GetPointer() );
    if( m_CurrentTubeId == -1 || curTube->GetId() == m_CurrentTubeId )
      {
      typename TubeType::PointType avgPos;
      std::vector< double > avgTagScalar;

      typename TubeType::TubePointListType::iterator pointItr;
      typename TubeType::TubePointListType::iterator tmpPointItr;

      typename TubeType::TubePointListType::iterator beginItr = curTube->GetPoints().begin();
      typename TubeType::TubePointListType::iterator endItr = curTube->GetPoints().end();

      typename TubeType::TubePointListType newPointList;

      std::vector< double > w;
      int wSize = 0;

      // Calculate the weighing window w
      if( smoothFunction == SMOOTH_TUBE_USING_INDEX_AVERAGE )
        {
        int maxIndex = static_cast< int >( h );
        wSize = 2 * maxIndex + 1;
        w.resize( wSize, 1.0 );
        }
      else
        {
        // Standard Deviation
        double sigma = h;
        // Consider the points until 3*sigma
        int maxIndex = static_cast< int >( 3*sigma );
        wSize = 2 * maxIndex + 1;
        w.resize( wSize, 0.0 );
        for( int i = 0; i <= maxIndex; i++ )
          {
          // The multiplication term 1/sigma*sqrt( 2*pi ) isn't necessary
          // since we normalize at the end by the sum of w
          w[maxIndex+i] = exp( -i*i/( 2.0*sigma*sigma ) );
          w[maxIndex-i] = w[maxIndex+i];
          }
        }
  
      // Apply the weighing window
      int count = 0;
      unsigned int pointDimension = TubeType::ObjectDimension;
      for( pointItr = beginItr; pointItr != endItr; ++pointItr )
        {
        typename TubeType::TubePointType newPoint = *pointItr;
        double wTotal = 0;
        avgPos.Fill( 0 );
        double avgIntensity = 0;
        double avgRadius = 0;
        double avgMedialness = 0;
        double avgRidgeness = 0;
        double avgRoundness = 0;
        double avgCurvature = 0;
        double avgLevelness = 0;
        unsigned int dictSize = pointItr->GetTagScalarDictionary().size();
        avgTagScalar.resize( dictSize, 0 );

        tmpPointItr = pointItr;
        int wCenter = ( wSize-1 )/2;
  
        // Place the tmpPointItr at the beginning of the window
        tmpPointItr -= std::min( count, wCenter );
  
        // Place the window iterator so that the window center
        // is aligned to the current point.
        int pos = std::max( wCenter - count, 0 );
  
        // Compute the average over the window, weighing with w
        while( pos < wSize && tmpPointItr != endItr )
          {
          for( unsigned int j=0; j<pointDimension; ++j )
            {
            avgPos[j] += w[pos] * tmpPointItr->GetPositionInObjectSpace()[j];
            }
          avgIntensity += w[pos] * tmpPointItr->GetIntensity();
          avgRadius += w[pos] * tmpPointItr->GetRadiusInObjectSpace();
          avgMedialness += w[pos] * tmpPointItr->GetMedialness();
          avgRidgeness += w[pos] * tmpPointItr->GetRidgeness();
          avgRoundness += w[pos] * tmpPointItr->GetRoundness();
          avgCurvature += w[pos] * tmpPointItr->GetCurvature();
          avgLevelness += w[pos] * tmpPointItr->GetLevelness();
          unsigned int d = 0;
          std::map<std::string,double>::iterator iter =
            tmpPointItr->GetTagScalarDictionary().begin();
          while( d < dictSize )
            {
            avgTagScalar[d] += w[pos] * iter->second;
            ++iter;
            ++d;
            }

          wTotal += w[pos];
          ++pos;
          ++tmpPointItr;
          }
  
        // Divide by sum of weights -> finish average
        if( wTotal > 0 )
          {
          for( unsigned int i=0; i<pointDimension; ++i )
            {
            avgPos[i] /= wTotal;
            }
          avgIntensity /= wTotal;
          avgRadius /= wTotal;
          avgMedialness /= wTotal;
          avgRidgeness /= wTotal;
          avgRoundness /= wTotal;
          avgCurvature /= wTotal;
          avgLevelness /= wTotal;
          for( unsigned int d=0; d<dictSize; ++d )
            {
            avgTagScalar[d] /= wTotal;
            }
          // Update the new point coordinates
          newPoint.SetPositionInObjectSpace( avgPos );
          newPoint.SetIntensity( avgIntensity );
          newPoint.SetRadiusInObjectSpace( avgRadius );
          newPoint.SetMedialness( avgMedialness );
          newPoint.SetRidgeness( avgRidgeness );
          newPoint.SetRoundness( avgRoundness );
          newPoint.SetCurvature( avgCurvature );
          newPoint.SetLevelness( avgLevelness );
          unsigned int d = 0;
          std::map<std::string,double>::iterator iter =
            pointItr->GetTagScalarDictionary().begin();
          while( d < dictSize )
            {
            newPoint.SetTagScalarValue( iter->first, avgTagScalar[d] );
            ++iter;
            ++d;
            }
          }
        newPointList.push_back( newPoint );
        ++count;
        }
    
      curTube->SetPoints( newPointList );
      curTube->ComputeTangentsAndNormals();
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}


/**
 * Smooth a tube */
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SmoothTubeProperty( const std::string & propertyId,
  double h, SmoothTubeFunctionEnum smoothFunction )
{
  if( h <= 0 )
    {
    return;
    }

  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
      itCurTube->GetPointer() );
    if( m_CurrentTubeId == -1 || curTube->GetId() == m_CurrentTubeId )
      {
      typename TubeType::TubePointListType::iterator pointItr;
      typename TubeType::TubePointListType::iterator tmpPointItr;

      typename TubeType::TubePointListType::iterator beginItr = curTube->GetPoints().begin();
      typename TubeType::TubePointListType::iterator endItr = curTube->GetPoints().end();

      typename TubeType::TubePointListType newPointList;

      std::vector< double > w;
      int wSize = 0;

      // Calculate the weighing window w
      if( smoothFunction == SMOOTH_TUBE_USING_INDEX_AVERAGE )
        {
        int maxIndex = static_cast< int >( h );
        wSize = 2 * maxIndex + 1;
        w.resize( wSize, 1.0 );
        }
      else
        {
        // Standard Deviation
        double sigma = h;
        // Consider the points until 3*sigma
        int maxIndex = static_cast< int >( 3*sigma );
        wSize = 2 * maxIndex + 1;
        w.resize( wSize, 0.0 );
        for( int i = 0; i <= maxIndex; i++ )
          {
          // The multiplication term 1/sigma*sqrt( 2*pi ) isn't necessary
          // since we normalize at the end by the sum of w
          w[maxIndex+i] = exp( -i*i/( 2.0*sigma*sigma ) );
          w[maxIndex-i] = w[maxIndex+i];
          }
        }
  
      // Apply the weighing window
      int count = 0;
      unsigned int pointDimension = TubeType::ObjectDimension;
      for( pointItr = beginItr; pointItr != endItr; ++pointItr )
        {
        typename TubeType::TubePointType newPoint = *pointItr;

        double wTotal = 0;

        typename TubeType::PointType avgPropertyPoint;
        typename TubeType::VectorType avgPropertyVector;
        typename TubeType::CovariantVectorType avgPropertyCoVector;
        avgPropertyPoint.Fill( 0 );
        avgPropertyVector.Fill( 0 );
        avgPropertyCoVector.Fill( 0 );

        double avgPropertyValue = 0;

        tmpPointItr = pointItr;
        int wCenter = ( wSize-1 )/2;
  
        // Place the tmpPointItr at the beginning of the window
        tmpPointItr -= std::min( count, wCenter );
  
        // Place the window iterator so that the window center
        // is aligned to the current point.
        int pos = std::max( wCenter - count, 0 );
  
        // Compute the average over the window, weighing with w
        while( pos < wSize && tmpPointItr != endItr )
          {
          if( propertyId == "Position" )
            {
            for( unsigned int j=0; j<pointDimension; ++j )
              {
              avgPropertyPoint[j] += w[pos] * tmpPointItr->GetPositionInObjectSpace()[j];
              }
            }
          else if( propertyId == "Tangent" )
            {
            for( unsigned int j=0; j<pointDimension; ++j )
              {
              avgPropertyVector[j] += w[pos] * tmpPointItr->GetTangentInObjectSpace()[j];
              }
            }
          else if( propertyId == "Normal1" )
            {
            for( unsigned int j=0; j<pointDimension; ++j )
              {
              avgPropertyCoVector[j] += w[pos] * tmpPointItr->GetNormal1InObjectSpace()[j];
              }
            }
          else if( propertyId == "Normal2" )
            {
            for( unsigned int j=0; j<pointDimension; ++j )
              {
              avgPropertyCoVector[j] += w[pos] * tmpPointItr->GetNormal2InObjectSpace()[j];
              }
            }
          else if( propertyId == "Intensity" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetIntensity();
            }
          else if( propertyId == "Radius" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetRadiusInObjectSpace();
            }
          else if( propertyId == "Medialness" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetMedialness();
            }
          else if( propertyId == "Ridgeness" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetRidgeness();
            }
          else if( propertyId == "Roundness" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetRoundness();
            }
          else if( propertyId == "Curvature" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetCurvature();
            }
          else if( propertyId == "Levelness" )
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetLevelness();
            }
          else
            {
            avgPropertyValue += w[pos] * tmpPointItr->GetTagScalarValue(
              propertyId );
            }

          wTotal += w[pos];
          ++pos;
          ++tmpPointItr;
          }
  
        // Divide by sum of weights -> finish average
        if( wTotal > 0 )
          {
          for( unsigned int i=0; i<pointDimension; ++i )
            {
            avgPropertyPoint[i] /= wTotal;
            avgPropertyVector[i] /= wTotal;
            avgPropertyCoVector[i] /= wTotal;
            }
          avgPropertyValue /= wTotal;

          if( propertyId == "Position" )
            {
            newPoint.SetPositionInObjectSpace( avgPropertyPoint );
            }
          else if( propertyId == "Tangent" )
            {
            newPoint.SetTangentInObjectSpace( avgPropertyVector );
            }
          else if( propertyId == "Normal1" )
            {
            newPoint.SetNormal1InObjectSpace( avgPropertyCoVector );
            }
          else if( propertyId == "Normal2" )
            {
            newPoint.SetNormal2InObjectSpace( avgPropertyCoVector );
            }
          else if( propertyId == "Intensity" )
            {
            newPoint.SetIntensity( avgPropertyValue );
            }
          else if( propertyId == "Radius" )
            {
            newPoint.SetRadiusInObjectSpace( avgPropertyValue );
            }
          else if( propertyId == "Medialness" )
            {
            newPoint.SetMedialness( avgPropertyValue );
            }
          else if( propertyId == "Ridgeness" )
            {
            newPoint.SetRidgeness( avgPropertyValue );
            }
          else if( propertyId == "Roundness" )
            {
            newPoint.SetRoundness( avgPropertyValue );
            }
          else if( propertyId == "Curvature" )
            {
            newPoint.SetCurvature( avgPropertyValue );
            }
          else if( propertyId == "Levelness" )
            {
            newPoint.SetLevelness( avgPropertyValue );
            }
          else
            {
            newPoint.SetTagScalarValue( propertyId, avgPropertyValue );
            }
          }
        newPointList.push_back( newPoint );
        ++count;
        }
    
      curTube->SetPoints( newPointList );
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}

/**
 * Subsample a tube */
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
SubsampleTube( int N )
{
  // Cannot subsample by 0 or 1 or less
  if( N <= 1 )
    {
    return;
    }

  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
      itCurTube->GetPointer() );
    if( m_CurrentTubeId == -1 || curTube->GetId() == m_CurrentTubeId )
      {
      typename TubeType::TubePointListType::iterator pointItr;

      typename TubeType::TubePointListType::iterator beginItr = curTube->GetPoints().begin();
      typename TubeType::TubePointListType::iterator endItr = curTube->GetPoints().end();

      typename TubeType::TubePointListType newPointList;

      int count = 0;
      for( pointItr = beginItr; pointItr != endItr; ++pointItr )
        {
        typename TubeType::TubePointType newPoint = *pointItr;
    
        // An offset of N/2 on each side of the tube is chosen to
        // delete the end and the beginning of it,
        // which usually aren't very good.
        if( ( count - N/2 ) % N == 0 )
          {
          newPointList.push_back( newPoint );
          }
    
        ++count;
        }
    
      curTube->SetPoints( newPointList );
      curTube->ComputeTangentsAndNormals();
      }
    }

  inputTubeList->clear();
  delete inputTubeList;
}
    
/**
 * Compute tube length */
template< unsigned int DimensionT, class ImagePixelT >
double
TubeMathFilters< DimensionT, ImagePixelT >::
ComputeTubeLength( void )
{
  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  double tubeLength = 0;
  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
      itCurTube->GetPointer() );
    if( m_CurrentTubeId == -1 || curTube->GetId() == m_CurrentTubeId )
      {
      if( curTube->GetNumberOfPoints() <= 1 )
        {
        continue;
        }
      typedef typename PositionType::VectorType PositionVectorType;
    
      typename TubePointListType::const_iterator itTubePoints =
        curTube->GetPoints().begin();

      PositionVectorType ptPrevPosVec =
        itTubePoints->GetPositionInObjectSpace().GetVectorFromOrigin();

      ++itTubePoints;

      while( itTubePoints != curTube->GetPoints().end() )
        {
        PositionVectorType ptCurPosVec =
          itTubePoints->GetPositionInObjectSpace().GetVectorFromOrigin();
    
        tubeLength += ( ptCurPosVec - ptPrevPosVec ).GetNorm();
    
        ptPrevPosVec = ptCurPosVec;
        ++itTubePoints;
        }
      }
    }
    
  inputTubeList->clear();
  delete inputTubeList;

  return tubeLength;
}


/**
 * Assign each point in a tube a unique (sequential) Id */
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
RenumberTubes( void )
{
  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  unsigned int tubeCount = 0;
  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    (*itCurTube)->SetId( tubeCount++ );
    }
  inputTubeList->clear();
  delete inputTubeList;
}


/**
 * Assign each point in a tube a unique (sequential) Id */
template< unsigned int DimensionT, class ImagePixelT >
void
TubeMathFilters< DimensionT, ImagePixelT >::
RenumberPoints( void )
{
  TubeListPointerType inputTubeList = m_InputTubeGroup->GetChildren(
    m_InputTubeGroup->GetMaximumDepth(), "Tube" );

  for( typename TubeGroupType::ChildrenListType::iterator itCurTube =
    inputTubeList->begin(); itCurTube != inputTubeList->end(); ++itCurTube )
    {
    typename TubeType::Pointer curTube = dynamic_cast< TubeType * >(
      itCurTube->GetPointer() );
    if( m_CurrentTubeId == -1 || curTube->GetId() == m_CurrentTubeId )
      {
      if( curTube->GetNumberOfPoints() <= 1 )
        {
        continue;
        }
      typename TubePointListType::iterator itTubePoints =
        curTube->GetPoints().begin();

      unsigned int pointCount = 0;
      while( itTubePoints != curTube->GetPoints().end() )
        {
        itTubePoints->SetId( pointCount++ );
        ++itTubePoints;
        }
      //curTube->SetPoints( tubePointsList );
      }
    }
  inputTubeList->clear();
  delete inputTubeList;
}

} // End namespace tube

#endif // End !defined( __tubeTubeMathFilters_hxx )
