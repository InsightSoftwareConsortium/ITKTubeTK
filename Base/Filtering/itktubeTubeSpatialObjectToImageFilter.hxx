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
  m_FallOff = 0.0;

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
    std::cout << "WARNING: Size not set." << std::endl;
    std::cout << "   Reverting to an incorrect method to compute region."
      << std::endl;
    SizeType size;

    typename SuperClass::InputSpatialObjectType::BoundingBoxType::PointType
      maxPoint;
    maxPoint = InputTube->GetBoundingBox()->GetMaximum();

    typename OutputImageType::PointType   physicalSize;

    unsigned int buffer = 4;

    for( unsigned int i=0; i<ObjectDimension; i++ )
      {
      maxPoint[i] = maxPoint[i] *
                    ( InputTube->GetIndexToObjectTransform() )
                      ->GetScaleComponent()[i];
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

  typename OutputImageType::IndexType index;
  index.Fill( 0 );
  region.SetIndex( index );

  OutputImage->SetLargestPossibleRegion( region );  //
  OutputImage->SetBufferedRegion( region );         // set the region
  OutputImage->SetRequestedRegion( region );        //
  OutputImage->SetSpacing( this->m_Spacing );
  OutputImage->SetOrigin( this->m_Origin );
  OutputImage->Allocate();
  OutputImage->FillBuffer( 0 );

  itk::ContinuousIndex<double, ObjectDimension> point;

  m_RadiusImage = this->GetRadiusImage();
  //Build radius image for processing
  if( m_BuildRadiusImage )
    {
    m_RadiusImage->SetRegions( region );
    m_RadiusImage->SetSpacing( this->m_Spacing );
    m_RadiusImage->SetOrigin( this->m_Origin );
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

  typename OutputImageType::IndexType index2;

  while( TubeIterator != tubeList->end() )
    {
    TubeType * tube = ( TubeType * )TubeIterator->GetPointer();

    typename TubeType::TransformType * tubeIndexPhysTransform =
      tube->GetIndexToWorldTransform();

    // Force the computation of the tangents
    if( m_BuildTangentImage )
      {
      tube->RemoveDuplicatePoints();
      tube->ComputeTangentAndNormals();
      }

    for( unsigned int k=0; k < tube->GetNumberOfPoints(); k++ )
      {
      typedef typename TubeType::TubePointType TubePointType;
      const TubePointType* tubePoint = static_cast<const TubePointType*>(
        tube->GetPoint( k ) );
      OutputImage->TransformPhysicalPointToContinuousIndex(
	tubeIndexPhysTransform->TransformPoint( tubePoint->GetPosition() ),
	point );
      for( unsigned int i=0; i<ObjectDimension; i++ )
        {
        index[i] = ( long int )( point[i]+0.5 );
        }
      bool IsInside = OutputImage->GetLargestPossibleRegion().IsInside(index);

      if( IsInside )
        {
        // Density Image
        if( m_Cumulative )
          {
          OutputImage->SetPixel( index, OutputImage->GetPixel( index )+1 );
          }
        else
          {
          OutputImage->SetPixel( index, 1 );
          }

        // Tangent Image
        if( m_BuildTangentImage )
          {
          // Convert the tangent type to the actual tangent image pixel type
          typename TubeType::VectorType t = tubePoint->GetTangent();
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
          double phys_pt_radius = tubePoint->GetRadius() *
                                      tube
                                      ->GetIndexToObjectTransform()
                                      ->GetScaleComponent()[0];
          if( m_BuildRadiusImage )
            {
            m_RadiusImage->SetPixel( index,
                                  static_cast<RadiusPixelType>(
                                      phys_pt_radius ) );
            }

          long radius = ( long int )( phys_pt_radius / this->m_Spacing[0] );


          double step = radius/2;
          while( step > 1 )
            {
            step /= 2;
            }
          if( step < 0.5 )
            {
            step = 0.5;
            }
          if( ObjectDimension == 2 )
            {
            for( double x=-radius; x<=radius+step/2; x+=step )
              {
              for( double y=-radius; y<=radius+step/2; y+=step )
                {
                if( ( ( x*x ) +( y*y ) ) <= ( radius*radius ) )
                  // test  inside the sphere
                  {
                  index2[0]=( long )( point[0]+x+0.5 );
                  index2[1]=( long )( point[1]+y+0.5 );
                  if( OutputImage->GetLargestPossibleRegion().IsInside( index2 ) )
                    {
                    typedef typename OutputImageType::PixelType PixelType;
                    if( m_Cumulative )
                      {
                      OutputImage->SetPixel( index2,
                                            ( PixelType )( OutputImage
                                                        ->GetPixel( index2 )
                                            + 0.5 ) );
                      }
                    else
                      {
                      OutputImage->SetPixel( index2, 1 );
                      }
                    if( m_BuildRadiusImage )
                      {
                      m_RadiusImage->SetPixel( index2, phys_pt_radius );
                      }
                    }
                  }
                }
              }
            }
          else if( ObjectDimension == 3 )
            {
            for( double x=-radius; x<=radius+step/2; x+=step )
              {
              for( double y=-radius; y<=radius+step/2; y+=step )
                {
                for( double z=-radius; z<=radius+step/2; z+=step )
                  {
                  if( ( ( x*x ) +( y*y ) +( z*z ) ) <= ( radius*radius ) )
                    // test  inside the sphere
                    {
                    index2[0]=( long )( point[0]+x+0.5 );
                    index2[1]=( long )( point[1]+y+0.5 );
                    index2[2]=( long )( point[2]+z+0.5 );

                    // Test that point is within the output image boundries
                    if( OutputImage->GetLargestPossibleRegion().IsInside( index2 ) )
                      {
                      OutputImage->SetPixel( index2, 1 );
                      if( m_BuildRadiusImage )
                        {
                        m_RadiusImage->SetPixel( index2, phys_pt_radius );
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
