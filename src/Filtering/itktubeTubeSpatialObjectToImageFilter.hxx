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

#ifndef __itktubeTubeSpatialObjectToImageFilter_hxx
#define __itktubeTubeSpatialObjectToImageFilter_hxx

#include "itktubeTubeSpatialObjectToImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>

#include <vnl/vnl_vector.h>

/** Constructor */
template< unsigned int ObjectDimension, class TOutputImage,
  class TRadiusImage, class TTangentImage >
TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage, TRadiusImage,
  TTangentImage>
::TubeSpatialObjectToImageFilter( void )
{
  m_UseRadius = false;
  m_Cumulative = false;
  m_BuildRadiusImage = false;
  m_BuildTangentImage = false;

  m_ColorByTubeId = true;
  m_ColorByPointId = false;
  m_ColorByRadius = false;
  m_ColorByRidgeness = false;
  m_ColorByMedialness = false;
  m_ColorByBranchness = false;
  m_ColorByCurvature = false;
  m_ColorByLevelness = false;
  m_ColorByRoundness = false;
  m_ColorByIntensity = false;

  // This is a little bit tricky since the 2nd an 3rd outputs are
  //   not always computed
  this->SetNumberOfRequiredOutputs( 3 );
  m_RadiusImage = TRadiusImage::New();
  this->SetNthOutput( 1, m_RadiusImage );

  m_TangentImage = TTangentImage::New();
  this->SetNthOutput( 2, m_TangentImage );
}

/** Destructor */
template< unsigned int ObjectDimension, class TOutputImage,
  class TRadiusImage, class TTangentImage >
TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage, TRadiusImage,
  TTangentImage >
::~TubeSpatialObjectToImageFilter( void )
{
}

/** Return the Radius Image */
template< unsigned int ObjectDimension, class TOutputImage,
  class TRadiusImage, class TTangentImage >
typename TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage,
  TRadiusImage, TTangentImage >::RadiusImagePointer
TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage, TRadiusImage,
  TTangentImage>
::GetRadiusImage( void )
{
  return dynamic_cast< TRadiusImage * >( this->ProcessObject::GetOutput( 1 ) );
}

/** Return the tangent Image */
template< unsigned int ObjectDimension, class TOutputImage,
  class TRadiusImage, class TTangentImage >
typename TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage,
  TRadiusImage, TTangentImage >::TangentImagePointer
TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage, TRadiusImage,
  TTangentImage >
::GetTangentImage( void )
{
  return dynamic_cast< TTangentImage * >(
    this->ProcessObject::GetOutput( 2 ) );
}

/** Update */
template< unsigned int ObjectDimension, class TOutputImage,
  class TRadiusImage, class TTangentImage >
void
TubeSpatialObjectToImageFilter< ObjectDimension, TOutputImage, TRadiusImage,
  TTangentImage >
::GenerateData( void )
{
  itkDebugMacro( << "TubeSpatialObjectToImageFilter::Update() called." );

  //Get the input and output pointers
  const typename SuperClass::InputSpatialObjectType   * InputTube =
    this->GetInput();
  typename SuperClass::OutputImagePointer             OutputImage =
    this->GetOutput();

  // Generate the image
  typename OutputImageType::RegionType region;
  if( this->m_Size[0] == 0 )
    {
    std::cout << "WARNING: itktubeTubeSpacialObjectToImageFilter: Size not set."
      << std::endl;
    std::cout << "   Reverting to an incorrect method to compute region."
      << std::endl;
    typename OutputImageType::SizeType size;

    typename SuperClass::InputSpatialObjectType::BoundingBoxType::PointType
      maxPoint;
    InputTube->ComputeFamilyBoundingBox( 9999, "Tube" );
    maxPoint = InputTube->GetFamilyBoundingBoxInWorldSpace()->GetMaximum();

    typename OutputImageType::PointType   physicalSize;

    unsigned int buffer = 4;

    for( unsigned int i=0; i<ObjectDimension; i++ )
      {
      maxPoint[i] = maxPoint[i];
      physicalSize[i] = maxPoint[i] - this->m_Origin[i];

      /** Get the origin point within the image so that the object
       * remains in the image **/
      size[i] = ( long unsigned int )(
        physicalSize[i] / this->m_Spacing[i] ) + buffer;
      }
    region.SetSize( size );
    }
  else
    {
    region.SetSize( this->m_Size );
    }

#if ITK_VERSION_MAJOR >= 5 
#if ITK_VERSION_MINOR > 3
  region.SetIndex( this->m_Index );
#endif
#endif

  OutputImage->SetRegions( region );
  OutputImage->SetSpacing( this->m_Spacing );
  OutputImage->SetOrigin( this->m_Origin );
  OutputImage->SetDirection( this->m_Direction );
  OutputImage->Allocate();
  OutputImage->FillBuffer( 0 );
  typedef typename OutputImageType::PixelType PixelType;

  typedef itk::ContinuousIndex<double, ObjectDimension> ContinuousIndexType;
  ContinuousIndexType pointI;

  m_RadiusImage = this->GetRadiusImage();
  //Build radius image for processing
  if( m_BuildRadiusImage )
    {
    m_RadiusImage->SetRegions( region );
    m_RadiusImage->SetSpacing( this->m_Spacing );
    m_RadiusImage->SetOrigin( this->m_Origin );
    m_RadiusImage->SetDirection( this->m_Direction );
    m_RadiusImage->Allocate();
    m_RadiusImage->FillBuffer( 0 );
    }

  m_TangentImage = this->GetTangentImage();
  //Build radius image for processing
  if( m_BuildTangentImage )
    {
    m_TangentImage->SetRegions( region );
    m_TangentImage->SetSpacing( this->m_Spacing );
    m_TangentImage->SetOrigin( this->m_Origin );
    m_TangentImage->SetDirection( this->m_Direction );
    m_TangentImage->Allocate();
    TangentPixelType v;
    v.Fill( 0 );
    m_TangentImage->FillBuffer( v );
    }

  // Get the list of tubes
  char tubeName[] = "Tube";
  ChildrenListType* tubeList = InputTube->GetChildren(
    this->m_ChildrenDepth, tubeName );

  //int size = tubeList->size();

  typedef typename ChildrenListType::iterator ChildrenIteratorType;
  ChildrenIteratorType TubeIterator = tubeList->begin();

  typename OutputImageType::IndexType index;
  typename OutputImageType::IndexType index2;

  while( TubeIterator != tubeList->end() )
    {
    TubeType * tube = ( TubeType * )TubeIterator->GetPointer();

    tube->Update();

    // Force the computation of the tangents
    if( m_BuildTangentImage )
      {
      tube->RemoveDuplicatePointsInObjectSpace();
      tube->ComputeTangentsAndNormals();
      }

    for( unsigned int k=0; k < tube->GetNumberOfPoints(); k++ )
      {
      typedef typename TubeType::TubePointType TubePointType;
      const TubePointType* tubePoint = static_cast<const TubePointType*>(
        tube->GetPoint( k ) );
      OutputImage->TransformPhysicalPointToContinuousIndex(
        tubePoint->GetPositionInWorldSpace(),
        pointI );
      for( unsigned int i=0; i<ObjectDimension; i++ )
        {
        index[i] = ( long int )( pointI[i]+0.5 );
        }
      bool IsInside = OutputImage->GetLargestPossibleRegion().IsInside( index );

      if( IsInside )
        {
        double val = 1;
        if( m_ColorByTubeId )
          {
          val = tube->GetId();
          }
        else if( m_ColorByPointId )
          {
          val = tube->GetId();
          if(tubePoint->GetId() > 0)
            {
            val += 1.0/(tubePoint->GetId()+1);
            }
          }
        else if( m_ColorByRadius )
          {
          val = tubePoint->GetRadiusInWorldSpace();
          }
        else if( m_ColorByRidgeness )
          {
          val = tubePoint->GetRidgeness();
          }
        else if( m_ColorByMedialness )
          {
          val = tubePoint->GetMedialness();
          }
        else if( m_ColorByBranchness )
          {
          val = tubePoint->GetBranchness();
          }
        else if( m_ColorByCurvature )
          {
          val = tubePoint->GetCurvature();
          }
        else if( m_ColorByLevelness )
          {
          val = tubePoint->GetLevelness();
          }
        else if( m_ColorByRoundness )
          {
          val = tubePoint->GetRoundness();
          }
        else if( m_ColorByIntensity )
          {
          val = tubePoint->GetIntensity();
          }
        if( m_Cumulative )
          {
          OutputImage->SetPixel( index, OutputImage->GetPixel( index )+val );
          }
        else
          {
          OutputImage->SetPixel( index, val );
          }

        // Tangent Image
        if( m_BuildTangentImage )
          {
          // Convert the tangent type to the actual tangent image pixel type
          typename TubeType::VectorType t = tubePoint->GetTangentInWorldSpace();
          TangentPixelType tp;
          for( unsigned int tpind = 0;tpind<ObjectDimension;tpind++ )
            {
            tp[tpind] = t[tpind];
            }
          m_TangentImage->SetPixel( index, tp );
          }

        // Radius Image and Density image with radius
        if( m_UseRadius )
          {
          RadiusPixelType radius = tubePoint->GetRadiusInWorldSpace();

          if( m_BuildRadiusImage )
            {
            m_RadiusImage->SetPixel( index, radius );
            }

          ContinuousIndexType radiusI;
          for( unsigned int i = 0; i < ObjectDimension; i++ )
            {
            radiusI[i] = radius / this->m_Spacing[i];
            if(radiusI[i] < 0.25)
              {
              radiusI[i] = 0.25;
              }
            }

          if( ObjectDimension == 2 )
            {
            for( double x=-radiusI[0]; x<=radiusI[0]; x += 0.5 )
              {
              double xr = x / radiusI[0];
              for( double y=-radiusI[1]; y<=radiusI[1]; y += 0.5 )
                {
                double yr = y / radiusI[1];
                if( ( (xr*xr)+(yr*yr) ) <= 1 )
                  // test  inside the sphere
                  {
                  index2[0]=( long )( pointI[0]+x+0.5 );
                  index2[1]=( long )( pointI[1]+y+0.5 );
                  if( OutputImage->GetLargestPossibleRegion().IsInside(
                    index2 ) )
                    {
                    if( m_Cumulative )
                      {
                      OutputImage->SetPixel( index2, ( PixelType )(
                          OutputImage->GetPixel( index2 ) + val ) );
                      }
                    else
                      {
                      OutputImage->SetPixel( index2, val );
                      }
                    if( m_BuildRadiusImage )
                      {
                      m_RadiusImage->SetPixel( index2, radius );
                      }
                    }
                  }
                }
              }
            }
          else if( ObjectDimension == 3 )
            {
            for( double x=-radiusI[0]; x<=radiusI[0]; x+=0.5 )
              {
              double xr = x / radiusI[0];
              for( double y=-radiusI[1]; y<=radiusI[1]; y+=0.5 )
                {
                double yr = y / radiusI[1];
                for( double z=-radiusI[2]; z<=radiusI[2]; z+=0.5 )
                  {
                  double zr = z / radiusI[2];
                  if( ( (xr*xr) + (yr*yr) +(zr*zr) ) <= 1 )
                    {
                    index2[0]=( long )( pointI[0]+x+0.5 );
                    index2[1]=( long )( pointI[1]+y+0.5 );
                    index2[2]=( long )( pointI[2]+z+0.5 );

                    // Test that point is within the output image boundries
                    if( OutputImage->GetLargestPossibleRegion().IsInside(
                      index2 ) )
                      {
                      if( m_Cumulative )
                        {
                        OutputImage->SetPixel( index2, ( PixelType )(
                            OutputImage->GetPixel( index2 ) + val ) );
                        }
                      else
                        {
                        OutputImage->SetPixel( index2, val );
                        }
                      if( m_BuildRadiusImage )
                        {
                        m_RadiusImage->SetPixel( index2, radius );
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    ++TubeIterator;
    }

  delete tubeList;

  itkDebugMacro( << "TubeSpatialObjectToImageFilter::Update() finished." );

} // End update function

#endif // End !defined( __itktubeTubeSpatialObjectToImageFilter_hxx )
