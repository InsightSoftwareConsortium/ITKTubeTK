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

#ifndef __itktubeAnisotropicDiffusiveSparseRegistrationFilter_hxx
#define __itktubeAnisotropicDiffusiveSparseRegistrationFilter_hxx


#include "itktubeAnisotropicDiffusiveSparseRegistrationFilter.h"

#include "itktubeDiffusiveRegistrationFilterUtils.h"
#include "tubeTubeMath.h"

#include <itkImageRegionSplitter.h>
#include <itktubeSmoothingRecursiveGaussianImageFilter.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkVersion.h>

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
AnisotropicDiffusiveSparseRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::AnisotropicDiffusiveSparseRegistrationFilter( void )
{
  // Initialize attributes to NULL
  m_BorderSurface                               = 0;
  m_TubeList                                    = 0;
  m_TubeSurface                                 = 0;
  m_NormalMatrixImage                           = 0;
  m_WeightStructuresImage                       = 0;
  m_WeightRegularizationsImage                  = 0;
  m_HighResolutionNormalMatrixImage             = 0;
  m_HighResolutionWeightStructuresImage         = 0;
  m_HighResolutionWeightRegularizationsImage    = 0;

  // Lambda/gamma used to calculate weight from distance
  m_Lambda  = 0.01;
  m_Gamma   = -1.0;
  //Use the ITKv4 Threading Model
  //  (call ThreadedGenerateData instead of DynamicThreadedGenerateData)
  this->DynamicMultiThreadingOff();
}

/**
 * PrintSelf
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  if( m_BorderSurface )
    {
    os << indent << "Border surface:" << std::endl;
    m_BorderSurface->Print( os );
    }
//  if( m_TubeList )
//    {
//    os << indent << "Tube list:" << std::endl;
//    m_TubeList->Print( os );
//    }
  if( m_TubeSurface )
    {
    os << indent << "Tube surface:" << std::endl;
    m_TubeSurface->Print( os );
    }
  if( m_NormalMatrixImage )
    {
    os << indent << "Normal vector image:" << std::endl;
    m_NormalMatrixImage->Print( os, indent );
    }
  if( m_WeightStructuresImage )
    {
    os << indent << "Weight structures image:" << std::endl;
    m_WeightStructuresImage->Print( os, indent );
    }
  if( m_WeightRegularizationsImage )
    {
    os << indent << "Weight regularizations image:" << std::endl;
    m_WeightRegularizationsImage->Print( os, indent );
    }
  os << indent << "lambda: " << m_Lambda << std::endl;
  os << indent << "gamma: " << m_Gamma << std::endl;
  if( m_HighResolutionNormalMatrixImage )
    {
    os << indent << "High resolution normal vector image:" << std::endl;
    m_HighResolutionNormalMatrixImage->Print( os, indent );
    }
  if( m_HighResolutionWeightStructuresImage )
    {
    os << indent << "High resolution weight structures image:"
      << std::endl;
    m_HighResolutionWeightStructuresImage->Print( os, indent );
    }
  if( m_HighResolutionWeightRegularizationsImage )
    {
    os << indent << "High resolution weight regularizations image:"
      << std::endl;
    m_HighResolutionWeightRegularizationsImage->Print( os, indent );
    }
}

/**
 * Setup the pointers for the deformation component images
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeDeformationComponentAndDerivativeImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // Setup pointers to the deformation component images - we have the
  // TANGENTIAL component, which is the entire deformation field, and the
  // NORMAL component, which is the deformation vectors projected onto their
  // normals
  // There are four terms that share these two deformation component images
  this->SetDeformationComponentImage( SMOOTH_TANGENTIAL,
    this->GetOutput() );
  this->SetDeformationComponentImage( PROP_TANGENTIAL, this->GetOutput() );

  DeformationFieldPointer normalDeformationField =
    DeformationFieldType::New();
  DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
    normalDeformationField, output );
  this->SetDeformationComponentImage( SMOOTH_NORMAL,
    normalDeformationField );
  this->SetDeformationComponentImage( PROP_NORMAL,
    normalDeformationField );

  // Setup the first and second order deformation component image
  // derivatives. The two TANGENTIAL and two NORMAL components share
  // images.
  int termOrder[4] = { SMOOTH_TANGENTIAL, PROP_TANGENTIAL, SMOOTH_NORMAL,
    PROP_NORMAL };
  int t = 0;
  ScalarDerivativeImagePointer firstOrder = nullptr;
  TensorDerivativeImagePointer secondOrder = nullptr;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for( int j = 0; j < this->GetNumberOfTerms(); j++ )
      {
      t = termOrder[j];
      if( t == SMOOTH_TANGENTIAL || t == SMOOTH_NORMAL )
        {
        firstOrder = ScalarDerivativeImageType::New();
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          firstOrder, output );
        secondOrder = TensorDerivativeImageType::New();
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          secondOrder, output );
        }
      this->SetDeformationComponentFirstOrderDerivative( t, i,
        firstOrder );
      this->SetDeformationComponentSecondOrderDerivative( t, i,
        secondOrder );
      }
    }

  // If required, allocate and compute the normal matrix and weight images
  this->SetupNormalMatrixAndWeightImages();
}

/**
 * All other initialization done before the initialize / calculate change /
 * apply update loop
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::SetupNormalMatrixAndWeightImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // Whether or not we must compute the normal vector and/or weight images
  bool computeNormals = !m_NormalMatrixImage;
  bool computeWeightStructures = !m_WeightStructuresImage;
  bool computeWeightRegularizations = !m_WeightRegularizationsImage;

  // If we have a template for image attributes, use it.  The normal and
  // weight
  // images will be stored at their full resolution.  The diffusion tensor,
  // deformation component, derivative and multiplication vector images are
  // recalculated every time Initialize() is called to regenerate them at
  // the correct resolution.
  FixedImagePointer highResolutionTemplate =
    this->GetHighResolutionTemplate();

  // If we don't have a template:
  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // Compute the normal vector and/or weight images if required
  if( computeNormals || computeWeightStructures
      || computeWeightRegularizations )
    {
    // Ensure we have a border surface or a tube list to work with
    if( !this->GetBorderSurface() && !this->GetTubeList() )
      {
      itkExceptionMacro( << "You must provide a border surface, or a tube "
        << "list, or both a normal matrix image and a weight "
        << "image" );
      }

    // Compute the normals for the surface
    if( this->GetBorderSurface() )
      {
      this->ComputeBorderSurfaceNormals();
      }
    // Compute the normals for the tube list
    if( this->GetTubeList() )
      {
      this->ComputeTubeNormals();
      }

    // Allocate the normal vector and/or weight images
    if( computeNormals )
      {
      m_NormalMatrixImage = NormalMatrixImageType::New();
      if( highResolutionTemplate )
        {
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          m_NormalMatrixImage, highResolutionTemplate );
        }
      else
        {
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          m_NormalMatrixImage, output );
        }
      }
    if( computeWeightStructures )
      {
      m_WeightStructuresImage = WeightMatrixImageType::New();
      if( highResolutionTemplate )
        {
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          m_WeightStructuresImage, highResolutionTemplate );
        }
      else
        {
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          m_WeightStructuresImage, output );
        }
      }
    if( computeWeightRegularizations )
      {
      m_WeightRegularizationsImage = WeightComponentImageType::New();
      if( highResolutionTemplate )
        {
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          m_WeightRegularizationsImage, highResolutionTemplate );
        }
      else
        {
        DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
          m_WeightRegularizationsImage, output );
        }
      }

    // Actually compute the normal vectors and/or weights
    this->ComputeNormalMatrixAndWeightImages( computeNormals,
      computeWeightStructures, computeWeightRegularizations );
    }

  // On the first iteration of the first level, the normal and weight
  // images
  // will contain the highest resolution images.  We need to save the high
  // resolution images, and then resample them down to correspond to this
  // level.
  // On subsequent iterations, we just do the resampling.

  // Set the high resolution images only once
  if( !m_HighResolutionNormalMatrixImage )
    {
    m_HighResolutionNormalMatrixImage = m_NormalMatrixImage;
    }
  if( !m_HighResolutionWeightStructuresImage )
    {
    m_HighResolutionWeightStructuresImage = m_WeightStructuresImage;
    }
  if( !m_HighResolutionWeightRegularizationsImage )
    {
    m_HighResolutionWeightRegularizationsImage =
      m_WeightRegularizationsImage;
    }

  // If we are using a template or getting an image from the user, we need
  // to
  // make sure that the attributes of the member images match those of the
  // current output, so that they can be used to calclulate the diffusion
  // tensors, deformation components, etc
  if( !DiffusiveRegistrationFilterUtils::CompareImageAttributes(
        m_NormalMatrixImage.GetPointer(), output.GetPointer() ) )
    {
    DiffusiveRegistrationFilterUtils::ResampleImageNearestNeighbor(
          m_HighResolutionNormalMatrixImage, output, m_NormalMatrixImage );
    assert( DiffusiveRegistrationFilterUtils::CompareImageAttributes(
             m_NormalMatrixImage.GetPointer(), output.GetPointer() ) );
    }
  if( !DiffusiveRegistrationFilterUtils::CompareImageAttributes(
        m_WeightStructuresImage.GetPointer(), output.GetPointer() ) )
    {
    DiffusiveRegistrationFilterUtils::ResampleImageNearestNeighbor(
      m_HighResolutionWeightStructuresImage, output,
      m_WeightStructuresImage );
    assert( DiffusiveRegistrationFilterUtils::CompareImageAttributes(
      m_WeightStructuresImage.GetPointer(), output.GetPointer() ) );
    }
  if( !DiffusiveRegistrationFilterUtils::CompareImageAttributes(
        m_WeightRegularizationsImage.GetPointer(), output.GetPointer() ) )
    {
    DiffusiveRegistrationFilterUtils::ResampleImageLinear(
          m_HighResolutionWeightRegularizationsImage,
          output,
          m_WeightRegularizationsImage );
    assert( DiffusiveRegistrationFilterUtils::CompareImageAttributes(
      m_WeightRegularizationsImage.GetPointer(), output.GetPointer() ) );
    }
}

/**
 * Compute the normals for the border surface
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeBorderSurfaceNormals( void )
{
  assert( m_BorderSurface );
  vtkSmartPointer< vtkPolyDataNormals > normalsFilter =
    vtkSmartPointer< vtkPolyDataNormals >::New();
  normalsFilter->ComputePointNormalsOn();
  normalsFilter->ComputeCellNormalsOff();
  //normalsFilter->SetFeatureAngle( 30 ); // TODO
#if VTK_MAJOR_VERSION > 5
  normalsFilter->SetInputData( m_BorderSurface );
#else
  normalsFilter->SetInput( m_BorderSurface );
#endif
  normalsFilter->Update();
  m_BorderSurface = normalsFilter->GetOutput();

  // Make sure we now have the normals
  if( !m_BorderSurface->GetPointData() )
    {
    itkExceptionMacro( << "Border surface does not contain point data." );
    }
  else if( !m_BorderSurface->GetPointData()->GetNormals() )
    {
    itkExceptionMacro(
      << "Border surface point data does not have normals." );
    }
}

/**
 * Compute the normals for tubes
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeTubeNormals( void )
{
  assert( m_TubeList );
  assert( !m_TubeSurface ); // We only want to compute this once

  if( m_TubeList->empty() )
    {
    itkExceptionMacro( << "Tube list does not contain tubes." );
    }

  // Iterate through the tubes to compute tangents and normals
  int numPoints = 0;
  for( typename TubeListType::iterator tubeIt = m_TubeList->begin();
       tubeIt != m_TubeList->end();
       ++tubeIt )
    {
    typename TubeType::Pointer tube
        = static_cast< TubeType * >( tubeIt->GetPointer() );
    ::tube::ComputeTubeTangentsAndNormals< TubeType >( tube );
    numPoints += tube->GetPoints().size();
    }

  // We store the tube point positions and two normals in a vtkPolyData,
  // so that
  // later we can search through them using a vtkPointLocator.  Otherwise,
  // to
  // determine the normal matrix and weightings later on we will have a
  // nested
  // loop iterating through each tube point for each voxel coordinate
  // - which
  // takes forever.  As far as I know, the ITK point locator is not in the
  // Slicer ITK that tubeTK relies on.

  // Setup the normal float arrays for the tubes
  vtkSmartPointer< vtkFloatArray > positionFloatArray =
    vtkSmartPointer< vtkFloatArray >::New();
  positionFloatArray->SetNumberOfComponents( ImageDimension );
  positionFloatArray->SetNumberOfTuples( numPoints );

  vtkSmartPointer< vtkFloatArray > normal1FloatArray =
    vtkSmartPointer< vtkFloatArray >::New();
  normal1FloatArray->SetNumberOfComponents( ImageDimension );
  normal1FloatArray->SetNumberOfTuples( numPoints );
  normal1FloatArray->SetName( "normal1" );

  vtkSmartPointer< vtkFloatArray > normal2FloatArray =
    vtkSmartPointer< vtkFloatArray >::New();
  normal2FloatArray->SetNumberOfComponents( ImageDimension );
  normal2FloatArray->SetNumberOfTuples( numPoints );
  normal2FloatArray->SetName( "normal2" );

  vtkSmartPointer< vtkFloatArray > radiusFloatArray =
    vtkSmartPointer< vtkFloatArray >::New();
  radiusFloatArray->SetNumberOfValues( numPoints );
  radiusFloatArray->SetName( "radius " );

  // Fill in the normal float array for the tubes
  TubePointType * point;
  typename TubePointType::PointType pointPosition;
  typename TubePointType::CovariantVectorType pointNormal1;
  typename TubePointType::CovariantVectorType pointNormal2;
  float radius;
  int pointCounter = 0;
  for( typename TubeListType::iterator tubeIt = m_TubeList->begin();
       tubeIt != m_TubeList->end();
       ++tubeIt )
    {
    typename TubeType::Pointer tube
        = static_cast< TubeType * >( tubeIt->GetPointer() );
    TubePointListType tubePointList = tube->GetPoints();
    for( typename TubePointListType::iterator pointIt =
      tubePointList.begin();
      pointIt != tubePointList.end();
      ++pointIt )
      {
      point = static_cast< TubePointType * >( &( *pointIt ) );
      pointPosition = point->GetPositionInObjectSpace();
      pointNormal1 = point->GetNormal1InObjectSpace();
      pointNormal2 = point->GetNormal2InObjectSpace();
      radius = point->GetRadiusInObjectSpace();
      positionFloatArray->SetTuple( pointCounter,
                                    pointPosition.GetDataPointer() );
      normal1FloatArray->SetTuple( pointCounter,
                                   pointNormal1.GetDataPointer() );
      normal2FloatArray->SetTuple( pointCounter,
                                   pointNormal2.GetDataPointer() );
      radiusFloatArray->SetValue( pointCounter, radius );
      pointCounter++;
      }
    }

  // Add the position array to a vtkPoints
  vtkSmartPointer< vtkPoints > points =
    vtkSmartPointer< vtkPoints >::New();
  points->SetData( positionFloatArray );

  // Add the normal arrays to a vtkFieldData
  vtkSmartPointer< vtkFieldData > fieldData =
    vtkSmartPointer< vtkFieldData >::New();
  fieldData->AddArray( normal1FloatArray );
  fieldData->AddArray( normal2FloatArray );
  fieldData->AddArray( radiusFloatArray );

  // Finally, we store everything within the vtkPolyData
  m_TubeSurface = vtkSmartPointer< BorderSurfaceType >::New();
  m_TubeSurface->SetPoints( points );
  m_TubeSurface->SetFieldData( fieldData );
}

/**
 * Computes the normal vectors and distances to the closest point
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetNormalsAndDistancesFromClosestSurfacePoint(
    bool computeNormals,
    bool computeWeightStructures,
    bool computeWeightRegularizations )
{
  // Setup the point locator and get the normals from the surface polydata
  vtkSmartPointer< vtkPointLocator > surfacePointLocator = nullptr;
  vtkFloatArray * surfaceNormalData = nullptr;
  if( this->GetBorderSurface() )
    {
    surfacePointLocator = vtkSmartPointer< vtkPointLocator >::New();
    surfacePointLocator->SetDataSet( m_BorderSurface );
    surfacePointLocator->Initialize();
    surfacePointLocator->BuildLocator();
    surfaceNormalData = static_cast< vtkFloatArray * >(
        m_BorderSurface->GetPointData()->GetNormals() );
    assert( surfaceNormalData );
    }

  // Create a vtk polydata representing the tube points and associated normals
  vtkSmartPointer< vtkPointLocator > tubePointLocator = nullptr;
  vtkFloatArray * tubeNormal1Data = nullptr;
  vtkFloatArray * tubeNormal2Data = nullptr;
  vtkFloatArray * tubeRadiusData = nullptr;
  if( this->GetTubeSurface() )
    {
    tubePointLocator = vtkSmartPointer< vtkPointLocator >::New();
    tubePointLocator->SetDataSet( m_TubeSurface );
    tubePointLocator->Initialize();
    tubePointLocator->BuildLocator();
    tubeNormal1Data = static_cast< vtkFloatArray * >(
        m_TubeSurface->GetFieldData()->GetArray( "normal1" ) );
    tubeNormal2Data = static_cast< vtkFloatArray * >(
        m_TubeSurface->GetFieldData()->GetArray( "normal2" ) );
    tubeRadiusData = static_cast< vtkFloatArray * >(
        m_TubeSurface->GetFieldData()->GetArray( "radius " ) );
    assert( tubeNormal1Data );
    assert( tubeNormal2Data );
    assert( tubeRadiusData );
    }

  // Set up struct for multithreaded processing.
  AnisotropicDiffusiveSparseRegistrationFilterThreadStruct str;
  str.Filter = this;
  str.SurfacePointLocator = surfacePointLocator;
  str.SurfaceNormalData = surfaceNormalData;
  str.TubePointLocator = tubePointLocator;
  str.TubeNormal1Data = tubeNormal1Data;
  str.TubeNormal2Data = tubeNormal2Data;
  str.TubeRadiusData = tubeRadiusData;
  str.NormalMatrixImageLargestPossibleRegion
      = m_NormalMatrixImage->GetLargestPossibleRegion();
  str.WeightStructuresImageLargestPossibleRegion
      = m_WeightStructuresImage->GetLargestPossibleRegion();
  str.WeightRegularizationsImageLargestPossibleRegion
      = m_WeightRegularizationsImage->GetLargestPossibleRegion();
  str.ComputeNormals = computeNormals;
  str.ComputeWeightStructures = computeWeightStructures;
  str.ComputeWeightRegularizations = computeWeightRegularizations;

  // Multithread the execution
  this->GetMultiThreader()->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  this->GetMultiThreader()->SetSingleMethod(
      this->GetNormalsAndDistancesFromClosestSurfacePointThreaderCallback,
      & str );
  this->GetMultiThreader()->SingleMethodExecute();

  // Explicitly call Modified on the normal and weight images here, since
  // ThreadedGetNormalsAndDistancesFromClosestSurfacePoint changes these
  // buffers
  // through iterators which do not increment the update buffer timestamp
  this->m_NormalMatrixImage->Modified();
  this->m_WeightStructuresImage->Modified();
  this->m_WeightRegularizationsImage->Modified();
}

/**
 * Calls ThreadedGetNormalsAndDistancesFromClosestSurfacePoint for
 * processing
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetNormalsAndDistancesFromClosestSurfacePointThreaderCallback(
  void * arg )
{
  int threadId =
    ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->WorkUnitID;
  int threadCount =
    ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->NumberOfWorkUnits;

  AnisotropicDiffusiveSparseRegistrationFilterThreadStruct * str
      = ( AnisotropicDiffusiveSparseRegistrationFilterThreadStruct * )
            ( ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // First find out how many pieces extent can be split into.
  // We don't want to use the SplitRequestedRegion method from
  // itk::ImageSource
  // because we might be calculating the normals and weights of a high res
  // template, where the image extent will not match that of the output
  typedef itk::ImageRegionSplitter< ImageDimension > SplitterType;
  typename SplitterType::Pointer splitter = SplitterType::New();

  int normalTotal = splitter->GetNumberOfSplits(
      str->NormalMatrixImageLargestPossibleRegion, threadCount );
  ThreadNormalMatrixImageRegionType splitNormalRegion = splitter->GetSplit(
      threadId, normalTotal, str->NormalMatrixImageLargestPossibleRegion );

  int weightMatrixTotal = splitter->GetNumberOfSplits(
      str->WeightStructuresImageLargestPossibleRegion, threadCount );
  ThreadWeightMatrixImageRegionType splitWeightMatrixRegion
      = splitter->GetSplit( threadId,
        weightMatrixTotal,
        str->WeightStructuresImageLargestPossibleRegion );

  int weightComponentTotal = splitter->GetNumberOfSplits(
      str->WeightRegularizationsImageLargestPossibleRegion, threadCount );
  ThreadWeightComponentImageRegionType splitWeightComponentMatrixRegion
      = splitter->GetSplit(
          threadId,
          weightComponentTotal,
          str->WeightRegularizationsImageLargestPossibleRegion );

  // Assert we could split all of the images equally
  assert( normalTotal == weightMatrixTotal );
  assert( normalTotal == weightComponentTotal );

  if( threadId < normalTotal )
    {
    str->Filter->ThreadedGetNormalsAndDistancesFromClosestSurfacePoint(
        str->SurfacePointLocator,
        str->SurfaceNormalData,
        str->TubePointLocator,
        str->TubeNormal1Data,
        str->TubeNormal2Data,
        str->TubeRadiusData,
        splitNormalRegion,
        splitWeightMatrixRegion,
        splitWeightComponentMatrixRegion,
        str->ComputeNormals,
        str->ComputeWeightStructures,
        str->ComputeWeightRegularizations,
        threadId );
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

/**
 * Does the actual work of computing the normal vectors and distances to
 * the
 * closest point given an initialized vtkPointLocator and the surface
 * border normals
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedGetNormalsAndDistancesFromClosestSurfacePoint(
    vtkPointLocator * surfacePointLocator,
    vtkFloatArray * surfaceNormalData,
    vtkPointLocator * tubePointLocator,
    vtkFloatArray * tubeNormal1Data,
    vtkFloatArray * tubeNormal2Data,
    vtkFloatArray * tubeRadiusData,
    ThreadNormalMatrixImageRegionType & normalRegionToProcess,
    ThreadWeightMatrixImageRegionType & weightMatrixRegionToProcess,
    ThreadWeightComponentImageRegionType & weightComponentRegionToProcess,
    bool computeNormals,
    bool computeWeightStructures,
    bool computeWeightRegularizations,
    int )
{
  assert( ( surfacePointLocator && surfaceNormalData )
          || ( tubePointLocator && tubeNormal1Data && tubeNormal2Data
               && tubeRadiusData ) );

  // Setup iterators over the normal vector and weight images
  NormalMatrixImageRegionType normalIt(
      m_NormalMatrixImage, normalRegionToProcess );
  WeightComponentImageRegionType weightRegularizationsIt(
      m_WeightRegularizationsImage, weightComponentRegionToProcess );
  WeightMatrixImageRegionType weightStructuresIt(
      m_WeightStructuresImage, weightMatrixRegionToProcess );

  // The normal vector image will hold the normal of the closest point of
  // the
  // surface polydata, and the weight image will be a function of the
  // distance
  // between the voxel and this closest point, and the weight structure
  // image
  // will be a function of the structure type

  itk::Point< double, ImageDimension >  imageCoord;
  imageCoord.Fill( 0 );

  vtkIdType                             surfaceId = -1;
  WeightComponentType                   surfaceDistance = 0;
  vtkIdType                             tubeId = -1;
  WeightComponentType                   tubeDistance = 0;

  NormalMatrixType                      normalMatrix;
  normalMatrix.Fill( 0 );

  WeightMatrixType                      surfaceWeightMatrix;
  surfaceWeightMatrix.Fill( 0 );
  surfaceWeightMatrix( 0, 0 ) = 1.0;
  WeightMatrixType                      tubeWeightMatrix;
  tubeWeightMatrix.Fill( 0 );
  tubeWeightMatrix( 0, 0 ) = 1.0;
  tubeWeightMatrix( 1, 1 ) = 1.0;

  // Determine the normals of and the distances to the nearest border point
  for( normalIt.GoToBegin(),
       weightRegularizationsIt.GoToBegin(), weightStructuresIt.GoToBegin();
       !normalIt.IsAtEnd();
       ++normalIt, ++weightRegularizationsIt, ++weightStructuresIt )
    {
    surfaceDistance = 100000000.0;
    tubeDistance = 100000000.0;

    // Find the id of the closest surface point to the current voxel
    m_NormalMatrixImage->TransformIndexToPhysicalPoint( normalIt.GetIndex(),
                                                        imageCoord );
    if( surfacePointLocator )
      {
      surfaceId = surfacePointLocator->FindClosestPoint(
          imageCoord.GetDataPointer() );
      double borderCoord[ImageDimension];
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        borderCoord[i] = 0.0;
        }
      m_BorderSurface->GetPoint( surfaceId, borderCoord );
      surfaceDistance = 0.0;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        surfaceDistance += std::pow( imageCoord[i] - borderCoord[i], 2 );
        }
      surfaceDistance = std::sqrt( surfaceDistance );
      }
    if( tubePointLocator )
      {
      tubeId = tubePointLocator->FindClosestPoint(
          imageCoord.GetDataPointer() );
      double centerlineCoord[ImageDimension];
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        centerlineCoord[i] = 0.0;
        }
      m_TubeSurface->GetPoint( tubeId, centerlineCoord );

      // We want the distance to the tube surface, not the centerline point
      // Project the current index coordinate onto the plane defined by the
      // tube centerline point and its two normals, and consider the
      // distance
      // to the surface via the radius.
      float normal1[ImageDimension];
      float normal2[ImageDimension];
#if VTK_MAJOR_VERSION < 7
      tubeNormal1Data->GetTupleValue( tubeId, normal1 );
      tubeNormal2Data->GetTupleValue( tubeId, normal2 );
#else
      tubeNormal1Data->GetTypedTuple( tubeId, normal1 );
      tubeNormal2Data->GetTypedTuple( tubeId, normal2 );
#endif

      double distanceToCenterCoord = ComputeDistanceToPointOnPlane(
            centerlineCoord, normal1, normal2, imageCoord );
      tubeDistance =
        std::abs( distanceToCenterCoord
          - tubeRadiusData->GetValue( tubeId ) );
      }

    // Find the normal of the surface point that is closest to the
    // current voxel
    if( computeNormals )
      {
      normalMatrix.Fill( 0 );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        if( surfaceDistance <= tubeDistance )
          {
          normalMatrix( i, 0 ) =
            surfaceNormalData->GetValue( surfaceId * ImageDimension + i );
          }
        else
          {
          normalMatrix( i, 0 ) =
            tubeNormal1Data->GetValue( tubeId * ImageDimension + i );
          normalMatrix( i, 1 ) =
            tubeNormal2Data->GetValue( tubeId * ImageDimension + i );
          }
        }
      normalIt.Set( normalMatrix );
      }

    // Calculate distance between the current coordinate and the border
    // surface
    // coordinate
    if( computeWeightRegularizations )
      {
      if( surfaceDistance <= tubeDistance )
        {
        weightRegularizationsIt.Set( surfaceDistance );
        }
      else
        {
        weightRegularizationsIt.Set( tubeDistance );
        }
      }

    // Determine the weight structures based on the structure type
    if( computeWeightStructures )
      {
      if( surfaceDistance <= tubeDistance )
        {
        weightStructuresIt.Set( surfaceWeightMatrix );
        }
      else
        {
        weightStructuresIt.Set( tubeWeightMatrix );
        }
      }
    }
}

/**
 * Computes the distance from the point to the planePoint in the plane
 * defined
 * by two tangent vectors.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
double
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDistanceToPointOnPlane( double * planePoint,
   float * tangentVector1,
   float * tangentVector2,
   itk::Point< double, ImageDimension> otherPoint ) const
{
  double relativePoint[ImageDimension];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    relativePoint[i] = otherPoint[i] - planePoint[i];
    }

  double projection1 = 0.0;
  double projection2 = 0.0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    projection1 += std::pow( relativePoint[i] * tangentVector1[i], 2 );
    projection2 += std::pow( relativePoint[i] * tangentVector2[i], 2 );
    }
  projection1 = std::sqrt( projection1 );
  projection2 = std::sqrt( projection2 );

  double pointOnPlane[ImageDimension];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    pointOnPlane[i] = projection1 * tangentVector1[i]
        + projection2 * tangentVector2[i];
    }

  double distance = 0.0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    distance += std::pow( pointOnPlane[i], 2 );
    }
  distance = std::sqrt( distance );
  return distance;
}

/**
 * Updates the border normals and the weighting factor w
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeNormalMatrixAndWeightImages( bool computeNormals,
                                      bool computeWeightStructures,
                                      bool computeWeightRegularizations )
{
  assert( this->GetComputeRegularizationTerm() );
  //assert( m_BorderSurface->GetPointData()->GetNormals()
  //|| m_TubeSurface );
  // TODO put back
  assert( m_NormalMatrixImage );
  assert( m_WeightStructuresImage );
  assert( m_WeightRegularizationsImage );

  std::cout << "Computing normals and weights... " << std::endl;

  // The normal vector image will hold the normal of the closest point of
  // the
  // surface polydata, and the weight image will be a function of the
  // distance
  // between the voxel and this closest point
  this->GetNormalsAndDistancesFromClosestSurfacePoint( computeNormals,
    computeWeightStructures, computeWeightRegularizations );

  // Smooth the normals to handle corners ( because we are choosing the
  // closest point in the polydata
  if( computeNormals )
    {
    //  typedef itk::RecursiveGaussianImageFilter
    //      < NormalVectorImageType, NormalVectorImageType >
    //      NormalSmoothingFilterType;
    //  typename NormalSmoothingFilterType::Pointer normalSmooth
    //      = NormalSmoothingFilterType::New();
    //  normalSmooth->SetInput( m_NormalVectorImage );
    //  double normalSigma = 3.0;
    //  normalSmooth->SetSigma( normalSigma );
    //  normalSmooth->Update();
    //  m_NormalVectorImage = normalSmooth->GetOutput();
    }

  // Smooth the distance image to avoid "streaks" from faces of the
  // polydata
  // ( because we are choosing the closest point in the polydata )
  if( computeWeightRegularizations )
    {
    double weightSmoothingSigma = 1.0;
    typedef itk::tube::SmoothingRecursiveGaussianImageFilter
        < WeightComponentImageType, WeightComponentImageType >
        WeightSmoothingFilterType;
    typename WeightSmoothingFilterType::Pointer weightSmooth
        = WeightSmoothingFilterType::New();
    weightSmooth->SetInput( m_WeightRegularizationsImage );
    weightSmooth->SetSigma( weightSmoothingSigma );
    weightSmooth->Update();
    m_WeightRegularizationsImage = weightSmooth->GetOutput();

    // Iterate through the weight image and compute the weight from the
    // distance
    WeightComponentType weight = 0;
    WeightComponentImageRegionType weightRegularizationsIt(
        m_WeightRegularizationsImage,
        m_WeightRegularizationsImage->GetLargestPossibleRegion() );
    bool useExponential = ( m_Gamma == -1.0 );
    if( useExponential )
      {
      for( weightRegularizationsIt.GoToBegin();
      !weightRegularizationsIt.IsAtEnd();
      ++weightRegularizationsIt )
        {
        weight = this->ComputeWeightFromDistanceExponential(
            weightRegularizationsIt.Get() );
        weightRegularizationsIt.Set( weight );
        }
      }
    else
      {
      for( weightRegularizationsIt.GoToBegin();
      !weightRegularizationsIt.IsAtEnd();
      ++weightRegularizationsIt )
        {
        weight = this->ComputeWeightFromDistanceDirac(
            weightRegularizationsIt.Get() );
        weightRegularizationsIt.Set( weight );
        }
      }
    }

  std::cout << "Finished computing normals and weights." << std::endl;
}

/**
 * Calculates the weighting between the anisotropic diffusive and diffusive
 * regularizations, based on a given distance from a voxel to the border,
 * using
 * exponential decay.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::WeightComponentType
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeWeightFromDistanceExponential( const WeightComponentType
  distance ) const
{
  return std::exp( -1.0 * m_Lambda * distance );
}

/**
 * Calculates the weighting between the anisotropic diffusive and diffusive
 * regularizations, based on a given distance from a voxel to the border,
 * using
 * a Dirac-shaped function
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::WeightComponentType
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeWeightFromDistanceDirac( const WeightComponentType
  distance ) const
{
  return 1.0 - ( 1.0 / ( 1.0 + m_Lambda * m_Gamma *
    std::exp( -1.0 * m_Lambda * distance * distance ) ) );
}

/**
 * Updates the diffusion tensor image before each run of the registration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( m_NormalMatrixImage );
  assert( m_WeightStructuresImage );
  assert( m_WeightRegularizationsImage );

  // For the anisotropic diffusive regularization, we need to setup the
  // tangential diffusion tensors and the normal diffusion tensors

  // Used to compute the tangential and normal diffusion tensor images:
  // P = NAN^T

  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
      DiffusionTensorImageRegionType;
  typedef itk::Matrix
      < DeformationVectorComponentType, ImageDimension, ImageDimension >
      MatrixType;

  NormalMatrixType        N;
  NormalMatrixType        N_transpose;
  WeightMatrixType        A;
  WeightComponentType     w;
  MatrixType              P;
  MatrixType              wP;   // wP
  MatrixType              wP_transpose;
  MatrixType              I_wP; // I-wP
  MatrixType              I_wP_transpose;

  MatrixType              smoothTangentialMatrix;  // ( I-wP )^T * ( I-wP )
  MatrixType              smoothNormalMatrix;      // ( I-wP )^T * wP
  MatrixType              propTangentialMatrix;    // ( wP )^T * ( I-wP )
  MatrixType              propNormalMatrix;        // ( wP )^T * wP

  DiffusionTensorType     smoothTangentialDiffusionTensor;
  DiffusionTensorType     smoothNormalDiffusionTensor;
  DiffusionTensorType     propTangentialDiffusionTensor;
  DiffusionTensorType     propNormalDiffusionTensor;

  // Setup iterators
  NormalMatrixImageRegionType normalIt = NormalMatrixImageRegionType(
    m_NormalMatrixImage, m_NormalMatrixImage->GetLargestPossibleRegion() );
  WeightMatrixImageRegionType weightStructuresIt
      = WeightMatrixImageRegionType(
          m_WeightStructuresImage,
          m_WeightStructuresImage->GetLargestPossibleRegion() );
  WeightComponentImageRegionType weightRegularizationsIt
      = WeightComponentImageRegionType(
      m_WeightRegularizationsImage,
      m_WeightRegularizationsImage->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType smoothTangentialTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( SMOOTH_TANGENTIAL ),
          this->GetDiffusionTensorImage( SMOOTH_TANGENTIAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType smoothNormalTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( SMOOTH_NORMAL ),
          this->GetDiffusionTensorImage( SMOOTH_NORMAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType propTangentialTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( PROP_TANGENTIAL ),
          this->GetDiffusionTensorImage( PROP_TANGENTIAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType propNormalTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( PROP_NORMAL ),
          this->GetDiffusionTensorImage( PROP_NORMAL )
          ->GetLargestPossibleRegion() );

  for( normalIt.GoToBegin(),
    weightStructuresIt.GoToBegin(), weightRegularizationsIt.GoToBegin(),
    smoothTangentialTensorIt.GoToBegin(), smoothNormalTensorIt.GoToBegin(),
    propTangentialTensorIt.GoToBegin(), propNormalTensorIt.GoToBegin();
    !smoothTangentialTensorIt.IsAtEnd();
    ++normalIt,
    ++weightStructuresIt, ++weightRegularizationsIt,
    ++smoothTangentialTensorIt, ++smoothNormalTensorIt,
    ++propTangentialTensorIt, ++propNormalTensorIt )
    {
    N = normalIt.Get();
    A = weightStructuresIt.Get();
    w = weightRegularizationsIt.Get();

    // The matrices are used for calculations, and will be copied to the
    // diffusion tensors afterwards.  The matrices are guaranteed to be
    // symmetric.
    P = N * A * N.GetTranspose(); // NAN^T
    wP = P * w; // wP
    I_wP.SetIdentity();
    I_wP = I_wP - wP; // I-wP
    I_wP_transpose = I_wP.GetTranspose(); // ( I-wP )^T
    smoothTangentialMatrix = I_wP_transpose * I_wP;
    // ( I-wP )^T * ( I-wP )
    smoothNormalMatrix = I_wP_transpose * wP; // ( I-wP )^T * wP
    wP_transpose = wP.GetTranspose(); // ( wP )^T
    propTangentialMatrix = wP_transpose * I_wP; // ( wP )^T * ( I-wP )
    propNormalMatrix = wP_transpose * wP; // ( wP )^T * wP

    // Copy the matrices to the diffusion tensor
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        smoothTangentialDiffusionTensor( i, j ) =
          smoothTangentialMatrix( i, j );
        smoothNormalDiffusionTensor( i, j ) = smoothNormalMatrix( i, j );
        propTangentialDiffusionTensor( i, j ) =
          propTangentialMatrix( i, j );
        propNormalDiffusionTensor( i, j ) = propNormalMatrix( i, j );
        }
      }

    // Copy the diffusion tensors to their images
    smoothTangentialTensorIt.Set( smoothTangentialDiffusionTensor );
    smoothNormalTensorIt.Set( smoothNormalDiffusionTensor );
    propTangentialTensorIt.Set( propTangentialDiffusionTensor );
    propNormalTensorIt.Set( propNormalDiffusionTensor );
    }
}

/** Computes the multiplication vectors that the div( Tensor /grad u )
 * values
 *  are multiplied by.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeMultiplicationVectorImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );
  assert( this->GetNormalMatrixImage() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // Allocate the images needed when using the anisotropic diffusive
  // regularization
  // There are no multiplication vectors for the smooth terms, and the
  // superclass will init them to zeros for us.
  // The prop multiplication vector is NAN_l

  // Iterate over the normal matrix and weight images
  NormalMatrixImageRegionType normalIt = NormalMatrixImageRegionType(
    m_NormalMatrixImage, m_NormalMatrixImage->GetLargestPossibleRegion() );
  WeightMatrixImageRegionType weightStructuresIt =
    WeightMatrixImageRegionType( m_WeightStructuresImage,
      m_WeightStructuresImage->GetLargestPossibleRegion() );

  DeformationVectorType multVector;
  multVector.Fill( 0.0 );
  NormalMatrixType N;
  N.Fill( 0.0 );
  WeightMatrixType A;
  A.Fill( 0.0 );
  DeformationVectorType N_l;
  N_l.Fill( 0.0 );

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    // Create the multiplication vector image that is shared between the
    // PROPs
    DeformationFieldPointer normalMultsImage = DeformationFieldType::New();
    DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
      normalMultsImage, output );

    // Calculate NAN_l
    DeformationVectorImageRegionType multIt =
      DeformationVectorImageRegionType( normalMultsImage,
        normalMultsImage->GetLargestPossibleRegion() );
    for( normalIt.GoToBegin(), weightStructuresIt.GoToBegin(),
      multIt.GoToBegin();
      !multIt.IsAtEnd();
      ++normalIt, ++weightStructuresIt, ++multIt )
      {
      multVector.Fill( 0.0 );
      N = normalIt.Get();
      A = weightStructuresIt.Get();
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        N_l[j] = N[i][j];
        }
      multVector = N * A * N_l;
      multIt.Set( multVector );
      }

    // Set the multiplication vector image
    this->SetMultiplicationVectorImage( PROP_TANGENTIAL, i,
      normalMultsImage );
    this->SetMultiplicationVectorImage( PROP_NORMAL, i, normalMultsImage );
    }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    assert( this->GetMultiplicationVectorImage( SMOOTH_TANGENTIAL, i )
            == this->GetMultiplicationVectorImage( SMOOTH_NORMAL, i ) );
    assert( this->GetMultiplicationVectorImage( PROP_TANGENTIAL, i )
            == this->GetMultiplicationVectorImage( PROP_NORMAL, i ) );
    }
}

/**
 * Updates the deformation vector component images before each iteration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationComponentImages( OutputImageType * output )
{
  assert( this->GetComputeRegularizationTerm() );

  // The normal component of u_l is ( NAN_l )^T * u
  // Conveniently, ( NAN_l ) is the PROP multiplication vector ( for both
  // SMOOTH
  // and PROP ), so we don't need to compute it again here
  DeformationVectorImageRegionArrayType NAN_lRegionArray;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    DeformationFieldPointer propMultImage
        = this->GetMultiplicationVectorImage( PROP_TANGENTIAL, i );
    assert( propMultImage );
    NAN_lRegionArray[i] = DeformationVectorImageRegionType(
        propMultImage, propMultImage->GetLargestPossibleRegion() );
    }

  // Setup the iterators for the deformation field
  OutputImageRegionType outputRegion( output,
                                     output->GetLargestPossibleRegion() );

  // We want to update the deformation component images for both
  // SMOOTH_NORMAL
  // and PROP_NORMAL, but they point to the same image, so we will grab the
  // pointer from SMOOTH_NORMAL to update both
  DeformationFieldPointer normalDeformationField
      = this->GetDeformationComponentImage( SMOOTH_NORMAL );
  OutputImageRegionType normalDeformationRegion(
      normalDeformationField,
      normalDeformationField->GetLargestPossibleRegion() );

  // Iterate over the deformation field and calculate the new normal
  // components
  DeformationVectorType  u; // deformation vector
  u.Fill( 0.0 );
  DeformationVectorType  normalDeformationVector;
  normalDeformationVector.Fill( 0.0 );
  DeformationVectorType  tangentialDeformationVector;
  tangentialDeformationVector.Fill( 0.0 );

  outputRegion.GoToBegin();
  normalDeformationRegion.GoToBegin();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    NAN_lRegionArray[i].GoToBegin();
    }

  while( !outputRegion.IsAtEnd() )
    {
    u = outputRegion.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      normalDeformationVector[i] = NAN_lRegionArray[i].Get() * u;
      }
    normalDeformationRegion.Set( normalDeformationVector );

    // Test that the normal and tangential components were computed corectly
    // ( they should be orthogonal )

    // tangential component = u - normal component
    tangentialDeformationVector = u - normalDeformationVector;

    if( normalDeformationVector * tangentialDeformationVector > 0.005 )
      {
      itkWarningMacro( << "Normal and tangential deformation field "
                       << "components are not orthogonal" );
      this->StopRegistration();
      }

    ++outputRegion;
    ++normalDeformationRegion;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      ++NAN_lRegionArray[i];
      }
    }
  normalDeformationField->Modified();

  assert( this->GetDeformationComponentImage( SMOOTH_TANGENTIAL )
          == this->GetDeformationComponentImage( PROP_TANGENTIAL ) );
  assert( this->GetDeformationComponentImage( SMOOTH_NORMAL )
          == this->GetDeformationComponentImage( PROP_NORMAL ) );
}

/**
 * Calculates the derivatives of the deformation vector derivatives after
 * each iteration.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDeformationComponentDerivativeImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // Get the spacing and the radius
  SpacingType spacing = this->GetOutput()->GetSpacing();
  const RegistrationFunctionType * df =
    this->GetRegistrationFunctionPointer();
  typename OutputImageType::SizeType radius = df->GetRadius();

  // Calculate first and second order deformation component image derivatives.
  // The two TANGENTIAL and two NORMAL components share images.  We will
  // calculate derivatives only for the SMOOTH terms, and the PROP terms
  // will be automatically updated since they point to the same image data.
  DeformationComponentImageArrayType deformationComponentImageArray;
  deformationComponentImageArray.Fill( nullptr );
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    if( i == SMOOTH_TANGENTIAL || i == SMOOTH_NORMAL )
      {
      DiffusiveRegistrationFilterUtils
        ::ExtractXYZComponentsFromDeformationField(
            this->GetDeformationComponentImage( i ),
            deformationComponentImageArray );

      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        this->ComputeDeformationComponentDerivativeImageHelper(
            deformationComponentImageArray[j], i, j, spacing, radius );
        }
      }
    }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    assert( this->GetDeformationComponentFirstOrderDerivative(
        SMOOTH_TANGENTIAL, i )
          == this->GetDeformationComponentFirstOrderDerivative(
              PROP_TANGENTIAL, i ) );
    assert( this->GetDeformationComponentFirstOrderDerivative(
        SMOOTH_NORMAL, i )
    == this->GetDeformationComponentFirstOrderDerivative(
              PROP_NORMAL, i ) );

    assert( this->GetDeformationComponentSecondOrderDerivative(
        SMOOTH_TANGENTIAL, i )
          == this->GetDeformationComponentSecondOrderDerivative(
              PROP_TANGENTIAL, i ) );
    assert( this->GetDeformationComponentSecondOrderDerivative(
        SMOOTH_NORMAL, i )
          == this->GetDeformationComponentSecondOrderDerivative(
              PROP_NORMAL, i ) );
    }
}

/**
 * Get the normal matrix image as a vector image.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetHighResolutionNormalVectorImage(
    NormalVectorImagePointer & normalImage,
    int dim,
    bool getHighResolutionNormalVectorImage ) const
{
  if( getHighResolutionNormalVectorImage
    && !m_HighResolutionNormalMatrixImage )
    {
    return;
    }
  if( !getHighResolutionNormalVectorImage && !m_NormalMatrixImage )
    {
    return;
    }

  // Allocate the vector image and iterate over the normal matrix image
  NormalMatrixImageRegionType normalMatrixIt;
  if( getHighResolutionNormalVectorImage )
    {
    DiffusiveRegistrationFilterUtils::AllocateSpaceForImage( normalImage,
                                 m_HighResolutionNormalMatrixImage );
    normalMatrixIt = NormalMatrixImageRegionType(
        m_HighResolutionNormalMatrixImage,
        m_HighResolutionNormalMatrixImage->GetLargestPossibleRegion() );
    }
  else
    {
    DiffusiveRegistrationFilterUtils::AllocateSpaceForImage( normalImage,
                                 m_NormalMatrixImage );
    normalMatrixIt = NormalMatrixImageRegionType(
        m_NormalMatrixImage,
        m_NormalMatrixImage->GetLargestPossibleRegion() );
    }
  typedef itk::ImageRegionIterator< NormalVectorImageType >
      NormalVectorImageRegionType;
  NormalVectorImageRegionType normalVectorIt = NormalVectorImageRegionType(
      normalImage, normalImage->GetLargestPossibleRegion() );

  NormalMatrixType matrix;
  matrix.Fill( 0.0 );
  NormalVectorType vector;
  vector.Fill( 0.0 );

  // Extract the normal vectors from the matrix
  for( normalMatrixIt.GoToBegin(), normalVectorIt.GoToBegin();
       !normalMatrixIt.IsAtEnd();
       ++normalMatrixIt, ++normalVectorIt )
    {
    matrix = normalMatrixIt.Get();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      vector[i] = matrix( i, dim );
      normalVectorIt.Set( vector );
      }
    }
}

} // End namespace tube

} // End namespace itk

// End !defined( __itktubeAnisotropicDiffusiveSparseRegistrationFilter_hxx )
#endif
