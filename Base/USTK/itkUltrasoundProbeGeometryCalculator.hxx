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

#ifndef __itkUltrasoundProbeGeometryCalculator_hxx
#define __itkUltrasoundProbeGeometryCalculator_hxx

#include "itkUltrasoundProbeGeometryCalculator.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionReverseConstIterator.h>
#include <itkListSample.h>
#include <itkStatisticsAlgorithm.h>

#include <vnl/algo/vnl_qr.h>

namespace itk
{

namespace tube
{

template< class TInputImage >
UltrasoundProbeGeometryCalculator< TInputImage >
::UltrasoundProbeGeometryCalculator( void )
  : m_GeneralBeamDirection( 1 ),
    m_BackgroundValue( NumericTraits< InputPixelType >::Zero )
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 2 );
  this->SetPrimaryOutputName( "UltrasoundProbeOrigin" );
  this->ProcessObject::SetOutput( "UltrasoundProbeOrigin",
    this->MakeOutput( 0 ) );
  this->ProcessObject::SetOutput( "StartOfAcquisitionRadius",
    this->MakeOutput( 1 ) );
}


template< class TInputImage >
UltrasoundProbeGeometryCalculator< TInputImage >
::~UltrasoundProbeGeometryCalculator( void )
{
}


template< class TInputImage >
DataObject::Pointer
UltrasoundProbeGeometryCalculator< TInputImage >
::MakeOutput( DataObjectPointerArraySizeType index )
{
  if( index == 0 )
    {
    return DecoratedOriginType::New().GetPointer();
    }
  return DecoratedRadiusType::New().GetPointer();
}


template< class TInputImage >
void
UltrasoundProbeGeometryCalculator< TInputImage >
::SetInput( const InputImageType * image )
{
  this->ProcessObject::SetPrimaryInput( const_cast< InputImageType * >( image ) );
}


template< class TInputImage >
const typename UltrasoundProbeGeometryCalculator< TInputImage >::InputImageType *
UltrasoundProbeGeometryCalculator< TInputImage >
::GetInput( void ) const
{
  return itkDynamicCastInDebugMode< const TInputImage * >( this->GetPrimaryInput() );
}


template< class TInputImage >
const typename UltrasoundProbeGeometryCalculator< TInputImage >::OriginType &
UltrasoundProbeGeometryCalculator< TInputImage >
::GetUltrasoundProbeOrigin( void ) const
{
  typename DecoratedOriginType::ConstPointer decoratedOrigin =
    static_cast< const DecoratedOriginType * >(
      this->ProcessObject::GetOutput( "UltrasoundProbeOrigin" ) );
  return decoratedOrigin->Get();
}


template< class TInputImage >
const typename UltrasoundProbeGeometryCalculator< TInputImage >::RadiusType &
UltrasoundProbeGeometryCalculator< TInputImage >
::GetStartOfAcquisitionRadius( void ) const
{
  typename DecoratedRadiusType::ConstPointer decoratedRadius =
    static_cast< const DecoratedRadiusType * >(
      this->ProcessObject::GetOutput( "StartOfAcquisitionRadius" ) );
  return decoratedRadius->Get();
}


template< class TInputImage >
void
UltrasoundProbeGeometryCalculator< TInputImage >
::GenerateData( void )
{
  OriginType probeOrigin;
  probeOrigin.Fill( 0.0 );
  RadiusType startOfAcquisitionRadius;

  typename InputImageType::ConstPointer inputImage = this->GetInput();
  const typename InputImageType::RegionType inputRegion =
    inputImage->GetLargestPossibleRegion();
  const typename InputImageType::SizeType   inputSize =
    inputRegion.GetSize();
  const typename InputImageType::IndexType  inputIndex =
    inputRegion.GetIndex();

  // Find the directions orthogonal to the general beam direction.
  unsigned int orthogonalDirections[InputImageType::ImageDimension - 1];
  for( unsigned int ii = 0, jj = 0; ii < InputImageType::ImageDimension; ++ii )
    {
    if( ii != m_GeneralBeamDirection )
      {
      orthogonalDirections[jj] = ii;
      ++jj;
      }
    }

  // Number of points to examine on each side of the imaged sector.  This could
  // be a class parameter.
  const unsigned int pointsToExamine = 8;
  typedef std::vector< OriginType > PointsContainerType;
  PointsContainerType pointsOnSide1( pointsToExamine );
  PointsContainerType pointsOnSide2( pointsToExamine );

  // The indices along the beam direction to walk in from along the orthogonal
  // direction.
  IndexValueType beamDirectionIndicesToExamine[pointsToExamine];
  for( unsigned int ii = 0; ii < pointsToExamine; ++ii )
    {
    // Fraction along each side of the beam direction to ignore in the image
    // because of padding, the sector geometry, etc.  This is also a potential
    // class parameter.
    const float ignoreFraction = 0.3f;
    const float keepFraction = 1.0f - 2*ignoreFraction;
    const float keepStep = keepFraction / ( pointsToExamine - 1 );
    const SizeValueType size = inputSize[m_GeneralBeamDirection];
    const SizeValueType offset =
      static_cast< SizeValueType >( size * ( ignoreFraction + keepStep * ii ) );
    beamDirectionIndicesToExamine[ii] = inputIndex[m_GeneralBeamDirection] + offset;
    }

  // Come in from the sides, and find the sides of the real data by locating
  // points that do not have the BackgroundValue
  // loop will work for 2D, 3D, probably not other D's... or E's
  for( unsigned int ii = 0; ii < InputImageType::ImageDimension - 1; ++ii )
    {
    // The orthogonal direction we are working with.  We work in the
    // GeneralBeamDirection x orthogonalDirection plane.
    const unsigned int orthogonalDirection = orthogonalDirections[ii];

    // For every point, walk in from the edge of the image, and find the side of
    // the sector image by locating the first non-background value.
    for( unsigned int pointIndex = 0; pointIndex < pointsToExamine; ++pointIndex )
      {
      // Define an ImageRegion that is a line along one of the
      // beamDirectionIndicesToExamine.
      typename InputImageType::RegionType lineRegion = inputRegion;
      typename InputImageType::IndexType  regionIndex = lineRegion.GetIndex();
      regionIndex[m_GeneralBeamDirection] = beamDirectionIndicesToExamine[pointIndex];
      typename InputImageType::SizeType   regionSize = lineRegion.GetSize();
      regionSize[m_GeneralBeamDirection] = 1;
      // use the middle plane in the 3D case
      for( unsigned int jj = 0; jj < InputImageType::ImageDimension; ++jj )
        {
        if( jj != m_GeneralBeamDirection && jj != orthogonalDirection )
          {
          regionIndex[jj] = regionIndex[jj] + regionSize[jj] / 2;
          regionSize[jj] = 1;
          }
        }
      lineRegion.SetIndex( regionIndex );
      lineRegion.SetSize( regionSize );

      // Walk in from one edge.
      typedef ImageRegionConstIterator< InputImageType > ImageIteratorType;
      ImageIteratorType imageIt( inputImage, lineRegion );
      for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
        {
        if( imageIt.Value() != m_BackgroundValue )
          {
          OriginType point;
          inputImage->TransformIndexToPhysicalPoint( imageIt.GetIndex(), point );
          pointsOnSide1[pointIndex] = point;
          break;
          }
        }

      // Walk in from the opposite edge.
      imageIt.GoToEnd();
      for( --imageIt; !imageIt.IsAtBegin(); --imageIt )
        {
        if( imageIt.Value() != m_BackgroundValue )
          {
          OriginType point;
          inputImage->TransformIndexToPhysicalPoint( imageIt.GetIndex(), point );
          pointsOnSide2[pointIndex] = point;
          break;
          }
        }
      } // end for each point to examine

    // Compute the line slope and intercept for each pair of points that were
    // previously detected.
    typedef typename OriginType::RealType                   MeasurementType;
    typedef Vector< MeasurementType, 2 >                    MeasurementVectorType;
    typedef Statistics::ListSample< MeasurementVectorType > SampleType;
    typename SampleType::Pointer side1LineParameters = SampleType::New();
    typename SampleType::Pointer side2LineParameters = SampleType::New();
    side1LineParameters->SetMeasurementVectorSize( 2 );
    side2LineParameters->SetMeasurementVectorSize( 2 );
    for( unsigned int pointIndex = 0; pointIndex < pointsToExamine - 1; ++pointIndex )
      {
      for( unsigned int pointIndexInc = pointIndex+1;
           pointIndexInc < pointsToExamine;
           ++pointIndexInc )
        {
        MeasurementVectorType side1LineParameter;
        MeasurementVectorType side2LineParameter;
        const MeasurementType slope1 =
          ( pointsOnSide1[pointIndexInc][orthogonalDirection] -
            pointsOnSide1[pointIndex][orthogonalDirection] ) /
          ( pointsOnSide1[pointIndexInc][m_GeneralBeamDirection] -
            pointsOnSide1[pointIndex][m_GeneralBeamDirection] );
        const MeasurementType slope2 =
          ( pointsOnSide2[pointIndexInc][orthogonalDirection] -
            pointsOnSide2[pointIndex][orthogonalDirection] ) /
          ( pointsOnSide2[pointIndexInc][m_GeneralBeamDirection] -
            pointsOnSide2[pointIndex][m_GeneralBeamDirection] );

        side1LineParameter[0] = slope1;
        side2LineParameter[0] = slope2;

        const MeasurementType intercept1 =
          pointsOnSide1[pointIndex][orthogonalDirection] -
            slope1 * pointsOnSide1[pointIndex][m_GeneralBeamDirection];
        const MeasurementType intercept2 =
          pointsOnSide2[pointIndex][orthogonalDirection] -
            slope2 * pointsOnSide2[pointIndex][m_GeneralBeamDirection];
        side1LineParameter[1] = intercept1;
        side2LineParameter[1] = intercept2;
        side1LineParameters->PushBack( side1LineParameter );
        side2LineParameters->PushBack( side2LineParameter );
        }
      }

    // Get the median line parameters.  Each slope and intercept are precise,
    // but since there may be some pixel within the imaged sector that happen to
    // have BackgroundValue values, they will not be detected and their line
    // parameters will be wrong -- it will have the wrong slope.  These outliers
    // are avoided by taking the median value.
    typedef Statistics::Subsample< SampleType > SubsampleType;
    typename SubsampleType::Pointer side1LineParametersSS = SubsampleType::New();
    typename SubsampleType::Pointer side2LineParametersSS = SubsampleType::New();
    side1LineParametersSS->SetSample( side1LineParameters );
    side2LineParametersSS->SetSample( side2LineParameters );
    side1LineParametersSS->InitializeWithAllInstances();
    side2LineParametersSS->InitializeWithAllInstances();
    const unsigned int activeDimension = 0;
    Statistics::Algorithm::InsertSort< SubsampleType >( side1LineParametersSS,
      activeDimension, 0, side1LineParametersSS->Size() );
    Statistics::Algorithm::InsertSort< SubsampleType >( side2LineParametersSS,
      activeDimension, 0, side2LineParametersSS->Size() );
    const MeasurementVectorType side1MedianLineParameter =
      side1LineParametersSS->GetMeasurementVectorByIndex(
        side1LineParametersSS->Size() / 2 );
    const MeasurementVectorType side2MedianLineParameter =
      side2LineParametersSS->GetMeasurementVectorByIndex(
        side2LineParametersSS->Size() / 2 );

    // Get the intersection of the lines -- defines the ProbeOrigin
    typedef vnl_vector< MeasurementType > VnlVectorType;
    VnlVectorType bVector( 2 );
    bVector[0] = side1MedianLineParameter[1];
    bVector[1] = side2MedianLineParameter[1];
    typedef vnl_matrix< MeasurementType > VnlMatrixType;
    VnlMatrixType mMatrix( 2, 2 );
    mMatrix[0][0] = 1.0;
    mMatrix[0][1] = -side1MedianLineParameter[0];
    mMatrix[1][0] = 1.0;
    mMatrix[1][1] = -side2MedianLineParameter[0];
    typedef vnl_qr< MeasurementType > VnlQRType;
    VnlQRType qr( mMatrix );
    const VnlVectorType planeOrigin = qr.solve( bVector );

    probeOrigin[orthogonalDirection] = planeOrigin[0];
    // Take the average in the GeneralBeamDirection.
    probeOrigin[m_GeneralBeamDirection] +=
      planeOrigin[1] / ( InputImageType::ImageDimension - 1 );
    } // end for each plane

  // Now, find the radius.
  const unsigned int radiusPointsToExamine = 5;

  // Just look at one orthogonal direction because the radius is assumed to be
  // the same in all planes.
  unsigned int orthogonalDirection = 0;
  for( unsigned int ii = 0; ii < InputImageType::ImageDimension; ++ii )
    {
    if( ii != m_GeneralBeamDirection )
      {
      orthogonalDirection = ii;
      }
    }

  // Look at points right along the center of the sector
  IndexValueType orthogonalIndicesToExamine[radiusPointsToExamine];
  for( unsigned int ii = 0; ii < radiusPointsToExamine; ++ii )
    {
    const float ignoreFraction = 0.9f;
    const float keepFraction = 1.0f - 2*ignoreFraction;
    const float keepStep = keepFraction / ( radiusPointsToExamine - 1 );
    const SizeValueType size = inputSize[orthogonalDirection];
    const SizeValueType offset = static_cast< SizeValueType >(
      size * ( ignoreFraction + keepStep * ii ) );
    orthogonalIndicesToExamine[ii] = inputIndex[m_GeneralBeamDirection] + offset;
    }

  // Compute the radius for each point found by casting along the beam
  // direction.
  typedef Vector< RadiusType, 1 > RadiusMeasurementVectorType;
  typedef Statistics::ListSample< RadiusMeasurementVectorType > RadiusSampleType;
  typename RadiusSampleType::Pointer radiiSamples = RadiusSampleType::New();
  for( unsigned int pointIndex = 0; pointIndex < radiusPointsToExamine; ++pointIndex )
    {
    // Create an ImageRegion along a line defined by the
    // orthogonalIndicesToExamine
    typename InputImageType::RegionType lineRegion = inputRegion;
    typename InputImageType::IndexType  regionIndex = lineRegion.GetIndex();
    regionIndex[orthogonalDirection] = orthogonalIndicesToExamine[pointIndex];
    typename InputImageType::SizeType   regionSize = lineRegion.GetSize();
    regionSize[orthogonalDirection] = 1;
    // use the middle plane in the 3D case
    for( unsigned int jj = 0; jj < InputImageType::ImageDimension; ++jj )
      {
      if( jj != m_GeneralBeamDirection && jj != orthogonalDirection )
        {
        regionIndex[jj] = regionIndex[jj] + regionSize[jj] / 2;
        regionSize[jj] = 1;
        }
      }
    lineRegion.SetIndex( regionIndex );
    lineRegion.SetSize( regionSize );

    // Walk along the beam direction and locate the sector as the first
    // non-BackgroundValue pixel.
    typedef ImageRegionConstIterator< InputImageType > ImageIteratorType;
    ImageIteratorType imageIt( inputImage, lineRegion );
    for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
      {
      if( imageIt.Value() != m_BackgroundValue )
        {
        OriginType point;
        inputImage->TransformIndexToPhysicalPoint( imageIt.GetIndex(), point );
        const typename OriginType::VectorType radiusVector = point - probeOrigin;
        radiiSamples->PushBack( radiusVector.GetNorm() );
        break;
        }
      }
    }

  // Use the median radii to avoid outliers due to BackgroundValue pixels in the
  // sector.
  typedef Statistics::Subsample< RadiusSampleType > RadiusSubsampleType;
  typename RadiusSubsampleType::Pointer radiiSubsample = RadiusSubsampleType::New();
  radiiSubsample->SetSample( radiiSamples );
  radiiSubsample->InitializeWithAllInstances();
  Statistics::Algorithm::InsertSort< RadiusSubsampleType >( radiiSubsample,
    0, 0, static_cast< int >( radiusPointsToExamine ) );
  startOfAcquisitionRadius =
    radiiSubsample->GetMeasurementVectorByIndex( radiusPointsToExamine / 2 )[0];

  // Set the values to our ProcessObject outputs.
  typename DecoratedOriginType::Pointer decoratedProbeOrigin =
    static_cast< DecoratedOriginType * >(
      this->ProcessObject::GetOutput( "UltrasoundProbeOrigin" ) );
  decoratedProbeOrigin->Set( probeOrigin );
  typename DecoratedRadiusType::Pointer decoratedStartOfAcquisitionRadius =
    static_cast< DecoratedRadiusType * >(
      this->ProcessObject::GetOutput( "StartOfAcquisitionRadius" ) );
  decoratedStartOfAcquisitionRadius->Set( startOfAcquisitionRadius );
}


template< class TInputImage >
void
UltrasoundProbeGeometryCalculator< TInputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "GeneralBeamDirection: "
     << m_GeneralBeamDirection << std::endl;
  os << indent << "BackgroundValue: "
     << m_BackgroundValue << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itkUltrasoundProbeGeometryCalculator_hxx )
