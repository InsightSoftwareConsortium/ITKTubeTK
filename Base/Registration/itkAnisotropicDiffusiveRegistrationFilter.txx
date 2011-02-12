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
#ifndef __itkAnisotropicDiffusiveRegistrationFilter_txx
#define __itkAnisotropicDiffusiveRegistrationFilter_txx

#include "itkAnisotropicDiffusiveRegistrationFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPolyDataNormals.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
AnisotropicDiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::AnisotropicDiffusiveRegistrationFilter()
{
  // Initialize attributes to NULL
  m_BorderSurface                               = 0;
  m_NormalVectorImage                           = 0;
  m_WeightImage                                 = 0;
  m_NormalDeformationField                      = 0;
  m_NormalDiffusionTensorImage                  = 0;
  m_NormalDiffusionTensorDerivativeImage        = 0;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_NormalDeformationComponentImages[i]       = 0;
    }

  // Lambda for exponential decay used to calculate weight from distance
  m_lambda = -0.01;
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  if( m_BorderSurface )
    {
    os << indent << "Border surface:" << std::endl;
    m_BorderSurface->Print( os );
    }
  if( m_NormalVectorImage )
    {
    os << indent << "Normal vector image:" << std::endl;
    m_NormalVectorImage->Print( os, indent );
    }
  if( m_WeightImage )
    {
    os << indent << "Weight image:" << std::endl;
    m_WeightImage->Print( os, indent );
    }
  if( m_NormalDiffusionTensorImage )
    {
    os << indent << "Normal diffusion tensor image:" << std::endl;
    m_NormalDiffusionTensorImage->Print( os, indent );
    }
  if( m_NormalDiffusionTensorDerivativeImage )
    {
    os << indent << "Normal diffusion tensor derivative image:" << std::endl;
    m_NormalDiffusionTensorDerivativeImage->Print( os, indent  );
    }
  if( m_NormalDeformationField )
    {
    os << indent << "Normal deformation field:" << std::endl;
    m_NormalDeformationField->Print( os, indent );
    }
  if( m_NormalDeformationComponentImages.Length != 0 )
    {
    os << indent << "Normal deformation component images:" << std::endl;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( m_NormalDeformationComponentImages[i] )
        {
        m_NormalDeformationComponentImages[i]->Print( os, indent );
        }
      }
    }
  os << indent << "lambda: " << m_lambda << std::endl;
}

/**
 * Create the registration function
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CreateRegistrationFunction()
{
  typename RegistrationFunctionType::Pointer registrationFunction =
      RegistrationFunctionType::New();
  registrationFunction->SetComputeRegularizationTerm( true );
  registrationFunction->SetComputeIntensityDistanceTerm( true );
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
      registrationFunction.GetPointer() ) );
}

/**
 * Get the registration function pointer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::RegistrationFunctionType *
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetRegistrationFunctionPointer() const
{
  RegistrationFunctionType * df = dynamic_cast< RegistrationFunctionType * >
       ( this->GetDifferenceFunction().GetPointer() );
  return df;
}

/**
 * All other initialization done before the initialize / calculate change /
 * apply update loop
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateImages()
{
  Superclass::AllocateImages();

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the images needed when using the anisotropic diffusive
  // regularization

  m_NormalDeformationField = OutputImageType::New();
  this->AllocateSpaceForImage( m_NormalDeformationField,
                               output );

  m_NormalDiffusionTensorImage = DiffusionTensorImageType::New();
  this->AllocateSpaceForImage( m_NormalDiffusionTensorImage,
                               output );

  m_NormalDiffusionTensorDerivativeImage = TensorDerivativeImageType::New();
  this->AllocateSpaceForImage( m_NormalDiffusionTensorDerivativeImage,
                               output );

  // If a normal vector image or weight image was supplied by the user, check
  // that it matches the output
  if( ( m_NormalVectorImage
        && !this->CompareImageAttributes( m_NormalVectorImage, output ) )
    || ( m_WeightImage
         && !this->CompareImageAttributes( m_WeightImage, output ) ) )
      {
    itkExceptionMacro( << "Normal vector image and/or weight image must have "
                       << "the same attributes as the output deformation "
                       << "field" );
    }

  // Whether or not we must compute the normal vector and/or weight images
  bool computeNormals = !m_NormalVectorImage;
  bool computeWeights = !m_WeightImage;

  // Compute the normal vector and/or weight images if required
  if( computeNormals || computeWeights )
    {
    // Ensure we have a border surface to work with
    if( !this->GetBorderSurface() )
      {
      itkExceptionMacro( << "You must provide a border surface, or both a "
                         << "normal vector image and a weight image" );
      }

    // Compute the normals for the surface
    this->ComputeBorderSurfaceNormals();

    // Allocate the normal vector and/or weight images
    if( computeNormals )
      {
      m_NormalVectorImage = NormalVectorImageType::New();
      this->AllocateSpaceForImage( m_NormalVectorImage, output );
      }
    if( computeWeights )
      {
      m_WeightImage = WeightImageType::New();
      this->AllocateSpaceForImage( m_WeightImage, output );
      }

    // Actually compute the normal vectors and/or weights
    this->ComputeNormalVectorAndWeightImages( computeNormals,
                                              computeWeights );
    }
}

/**
 * Compute the normals for the border surface
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeBorderSurfaceNormals()
{
  assert( m_BorderSurface );
  vtkPolyDataNormals * normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->ComputePointNormalsOn();
  normalsFilter->ComputeCellNormalsOff();
  //normalsFilter->SetFeatureAngle(30); // TODO
  normalsFilter->SetInput( m_BorderSurface );
  normalsFilter->Update();
  m_BorderSurface = normalsFilter->GetOutput();
  normalsFilter->Delete();

  // Make sure we now have the normals
  if ( !m_BorderSurface->GetPointData() )
    {
    itkExceptionMacro( << "Border surface does not contain point data" );
    }
  else if ( !m_BorderSurface->GetPointData()->GetNormals() )
    {
    itkExceptionMacro( << "Border surface point data does not have normals" );
    }
}
/**
 * Updates the border normals and the weighting factor w
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeNormalVectorAndWeightImages( bool computeNormals, bool computeWeights )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( m_BorderSurface->GetPointData()->GetNormals() );
  assert( m_NormalVectorImage );
  assert( m_WeightImage );

  std::cout << "Computing normals and weights... " << std::endl;

  // Setup iterators over the normal vector and weight images
  NormalVectorImageRegionType normalIt(
      m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );
  WeightImageRegionType weightIt(m_WeightImage,
                                 m_WeightImage->GetLargestPossibleRegion() );

  // Get the normals from the polydata
  vtkPointLocator * pointLocator = vtkPointLocator::New();
  pointLocator->SetDataSet( m_BorderSurface );
  vtkSmartPointer< vtkDataArray > normalData
      = m_BorderSurface->GetPointData()->GetNormals();

  // The normal vector image will hold the normal of the closest point of the
  // surface polydata, and the weight image will be a function of the distance
  // between the voxel and this closest point

  itk::Point< double, ImageDimension >  imageCoordAsPoint;
  imageCoordAsPoint.Fill( 0 );
  double                                imageCoord[ImageDimension];
  double                                borderCoord[ImageDimension];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    imageCoord[i] = 0;
    borderCoord[i] = 0;
    }
  vtkIdType                             id = 0;
  WeightType                            distance = 0;
  NormalVectorType                      normal;
  normal.Fill(0);
  WeightType                            weight = 0;

  // Determine the normals of and the distances to the nearest border
  for( normalIt.GoToBegin(), weightIt.GoToBegin();
       !normalIt.IsAtEnd();
       ++normalIt, ++weightIt )
    {

    // Find the normal of the surface point that is closest to the current voxel
    m_NormalVectorImage->TransformIndexToPhysicalPoint( normalIt.GetIndex(),
                                                        imageCoordAsPoint );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      imageCoord[i] = imageCoordAsPoint[i];
      }
    id = pointLocator->FindClosestPoint( imageCoord );
    normal = normalData->GetTuple( id );
    if( computeNormals )
      {
      normalIt.Set( normal );
      }

    // Calculate distance between the current coordinate and the border surface
    // coordinate
    m_BorderSurface->GetPoint( id, borderCoord );
    distance = 0.0;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      distance += pow( imageCoord[i] - borderCoord[i], 2 );
      }
    distance = sqrt( distance );

    // The weight image will temporarily store distances
    if( computeWeights )
      {
      weightIt.Set( distance );
      }
    }

  // Clean up memory
  pointLocator->Delete();

  // Smooth the normals to handle corners (because we are choosing the
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

  // Smooth the distance image to avoid "streaks" from faces of the polydata
  // (because we are choosing the closest point in the polydata)
  if( computeWeights )
    {
    double weightSmoothingSigma = 1.0;
    typedef itk::SmoothingRecursiveGaussianImageFilter
        < WeightImageType, WeightImageType > WeightSmoothingFilterType;
    typename WeightSmoothingFilterType::Pointer weightSmooth
        = WeightSmoothingFilterType::New();
    weightSmooth->SetInput( m_WeightImage );
    weightSmooth->SetSigma( weightSmoothingSigma );
    weightSmooth->Update();
    m_WeightImage = weightSmooth->GetOutput();

    // Iterate through the weight image and compute the weight from the
    WeightImageRegionType weightIt(
        m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );
    for( weightIt.GoToBegin(); !weightIt.IsAtEnd(); ++weightIt )
      {
      weight = this->ComputeWeightFromDistance( weightIt.Get() );
      weightIt.Set( weight );
      }
    }

  std::cout << "Finished computing normals and weights." << std::endl;
}

/**
 * Calculates the weighting between the anisotropic diffusive and diffusive
 * regularizations, based on a given distance from a voxel to the border
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::WeightType
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeWeightFromDistance( WeightType distance )
{
  WeightType weight = exp( m_lambda * distance );
  return weight;
}

/**
 * Updates the diffusion tensor image before each run of the registration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetTangentialDiffusionTensorImage() );
  assert( m_NormalVectorImage );
  assert( m_WeightImage );
  assert( m_NormalDiffusionTensorImage );

  // For the anisotropic diffusive regularization, we need to setup the
  // tangential diffusion tensors and the normal diffusion tensors

  // Used to compute the tangential and normal diffusion tensor images:
  // P = I - wnn^T
  // tangentialMatrix = tangentialDiffusionTensor = P^TP
  // normalMatrix = normalDiffusionTensor = w^2nn^T

  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
      DiffusionTensorImageRegionType;
  typedef itk::Matrix
      < DeformationVectorComponentType, ImageDimension, ImageDimension >
      MatrixType;

  NormalVectorType                              n;
  WeightType                                    w;
  MatrixType                                    P;
  MatrixType                                    normalMatrix;
  MatrixType                                    tangentialMatrix;
  typename DiffusionTensorImageType::PixelType  tangentialDiffusionTensor;
  typename DiffusionTensorImageType::PixelType  normalDiffusionTensor;

  // Setup iterators
  NormalVectorImageRegionType normalIt = NormalVectorImageRegionType(
      m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );
  WeightImageRegionType weightIt = WeightImageRegionType(
      m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType tangentialTensorIt
      = DiffusionTensorImageRegionType(
          this->GetTangentialDiffusionTensorImage(),
          this->GetTangentialDiffusionTensorImage()
            ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType normalTensorIt
      = DiffusionTensorImageRegionType(
          m_NormalDiffusionTensorImage,
          m_NormalDiffusionTensorImage->GetLargestPossibleRegion() );

  for( normalIt.GoToBegin(), weightIt.GoToBegin(),
       tangentialTensorIt.GoToBegin(), normalTensorIt.GoToBegin();
       !tangentialTensorIt.IsAtEnd();
       ++normalIt, ++weightIt, ++tangentialTensorIt, ++normalTensorIt )
    {
    n = normalIt.Get();
    w = weightIt.Get();

    // The matrices are used for calculations, and will be copied to the
    // diffusion tensors afterwards.  The matrices are guaranteed to be
    // symmetric.

    // Create the normalMatrix used to calculate nn^T
    // (The first column is filled with the values of n, the rest are 0s)
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      normalMatrix( i, 0 ) = n[i];
      for ( unsigned int j = 1; j < ImageDimension; j++ )
        {
        normalMatrix( i, j ) = 0;
        }
      }

    // Calculate the normal and tangential diffusion tensors
    normalMatrix = normalMatrix * normalMatrix.GetTranspose(); // nn^T
    normalMatrix = normalMatrix * w; // wnn^T
    P.SetIdentity();
    P = P - normalMatrix; // I - wnn^T
    tangentialMatrix = P.GetTranspose();
    tangentialMatrix = tangentialMatrix * P; // P^TP
    normalMatrix = normalMatrix * w; // w^2nn^T

    // Copy the matrices to the diffusion tensor
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        tangentialDiffusionTensor( i, j ) = tangentialMatrix( i, j );
        normalDiffusionTensor( i, j ) = normalMatrix( i, j );
        }
      }

    // Copy the diffusion tensors to their images
    tangentialTensorIt.Set( tangentialDiffusionTensor );
    normalTensorIt.Set( normalDiffusionTensor );
    }
}

/**
 * Updates the diffusion tensor image derivatives before each run of the
 * registration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorDerivativeImages()
{
  Superclass::ComputeDiffusionTensorDerivativeImages();

  // Compute the diffusion tensor derivative image for the normal plane
  this->ComputeDiffusionTensorDerivativeImageHelper(
      m_NormalDiffusionTensorImage, m_NormalDiffusionTensorDerivativeImage );
}

/**
 * Updates the deformation vector component images before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationVectorComponentImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetNormalVectorImage() );
  assert( this->GetWeightImage() );

  // Setup iterators
  NormalVectorImageRegionType normalVectorNeighborhood(
      m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );

  typename OutputImageType::Pointer output = this->GetOutput();
  OutputImageRegionType outputRegion(output,
                                     output->GetLargestPossibleRegion() );

  OutputImageRegionType normalDeformationRegion(
      m_NormalDeformationField,
      m_NormalDeformationField->GetLargestPossibleRegion() );

  // Calculate the tangential and normal components of the deformation field
  NormalVectorType       n;
  DeformationVectorType  u; // deformation vector
  DeformationVectorType  normalDeformationVector;
  DeformationVectorType  tangentialDeformationVector;

  for( normalVectorNeighborhood.GoToBegin(), outputRegion.GoToBegin(),
       normalDeformationRegion.GoToBegin();
  !outputRegion.IsAtEnd();
  ++normalVectorNeighborhood, ++outputRegion, ++normalDeformationRegion )
    {
    n = normalVectorNeighborhood.Get();
    u = outputRegion.Get();

    // normal component = (u^Tn)n
    normalDeformationVector = ( u * n ) * n;
    normalDeformationRegion.Set( normalDeformationVector );

    // Test that the normal and tangential components were computed corectly
    // (they should be orthogonal)

    // tangential component = u - normal component
    tangentialDeformationVector = u - normalDeformationVector;

    if( normalDeformationVector * tangentialDeformationVector > 0.005 )
      {
      itkExceptionMacro( << "Normal and tangential deformation field "
                         << "components are not orthogonal" );
      }
    }
  m_NormalDeformationField->Modified();

  // Update the tangential extracted components
  Superclass::UpdateDeformationVectorComponentImages();

  // Update the normal extracted components
  this->ExtractXYZComponentsFromDeformationField(
      m_NormalDeformationField, m_NormalDeformationComponentImages );
}

/**
 * Populates the update buffer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateChange()
{
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  // Not used during the calculate change step
  str.TimeStep = NumericTraits< TimeStepType >::Zero;
  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(
      this->CalculateChangeThreaderCallback, & str );

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  int threadCount = this->GetMultiThreader()->GetNumberOfThreads();

  str.TimeStepList = new TimeStepType[threadCount];
  str.ValidTimeStepList = new bool[threadCount];
  for ( int i = 0; i < threadCount; ++i )
    {
    str.ValidTimeStepList[i] = false;
    }

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Resolve the single value time step to return
  TimeStepType dt = this->ResolveTimeStep( str.TimeStepList,
                                           str.ValidTimeStepList,
                                           threadCount);
  delete [] str.TimeStepList;
  delete [] str.ValidTimeStepList;

  // Explicitly call Modified on the UpdateBuffer here, since
  // ThreadedCalculateChange changes this buffer through iterators which do not
  // increment the update buffer timestamp
  this->GetUpdateBuffer()->Modified();

  return dt;
}

/**
 * Calls ThreadedCalculateChange for processing
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_TYPE
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateChangeThreaderCallback( void * arg )
{
  int threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  int threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  DenseFDThreadStruct * str = (DenseFDThreadStruct *)
            (((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  int total;
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
                                             splitRegion );

  ThreadNormalVectorImageRegionType splitNormalVectorRegion;
  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
                                             splitNormalVectorRegion );

  ThreadDiffusionTensorImageRegionType splitTensorRegion;
  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
                                             splitTensorRegion );

  ThreadTensorDerivativeImageRegionType splitTensorDerivativeRegion;
  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
                                             splitTensorDerivativeRegion );

  ThreadDeformationVectorComponentImageRegionType
      splitDeformationComponentRegion;
  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
                                             splitDeformationComponentRegion );

  if (threadId < total)
    {
    str->TimeStepList[threadId] = str->Filter->ThreadedCalculateChange(
      splitRegion,
      splitNormalVectorRegion,
      splitTensorRegion,
      splitTensorDerivativeRegion,
      splitDeformationComponentRegion,
      threadId);
    str->ValidTimeStepList[threadId] = true;
    }

  return ITK_THREAD_RETURN_VALUE;
}

/**
 * Inherited from superclass - do not call
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChange(
    const ThreadRegionType &                      regionToProcess,
    const ThreadDiffusionTensorImageRegionType &  tensorRegionToProcess,
    const ThreadTensorDerivativeImageRegionType & derivativeRegionToProcess,
    const ThreadDeformationVectorComponentImageRegionType &
      deformationComponentRegionToProcess,
    int)
{
  // This function should never be called!
  itkExceptionMacro( << "ThreadedCalculateChange(regionToProcess, "
                     << "tensorRegionToProcess, derivativeRegionToProcess, "
                     << "deformationComponentRegionToProcess, threadId "
                     << "should never be called.  Use the other "
                     << "ThreadedCalculateChange function instead" );
}

/**
 * Does the actual work of calculating change over a region supplied by the
 * multithreading mechanism
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChange(
    const ThreadRegionType &                      regionToProcess,
    const ThreadNormalVectorImageRegionType &     normalVectorRegionToProcess,
    const ThreadDiffusionTensorImageRegionType &  tensorRegionToProcess,
    const ThreadTensorDerivativeImageRegionType & derivativeRegionToProcess,
    const ThreadDeformationVectorComponentImageRegionType &
      deformationComponentRegionToProcess,
    int)
{
  // Get the FiniteDifferenceFunction to use in calculations.
  RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );
  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  void * globalData = df->GetGlobalDataPointer();

  // Get the output and update buffer values
  const typename OutputImageType::SizeType radius = df->GetRadius();
  OutputImagePointer output = this->GetOutput();
  SpacingType spacing = output->GetSpacing();
  typename UpdateBufferType::Pointer updateBuffer = this->GetUpdateBuffer();

  // Get the tangential diffusion tensor image and its derivative image from the
  // superclass
  DiffusionTensorImagePointer tangentialDiffusionTensorImage
      = this->GetTangentialDiffusionTensorImage();
  TensorDerivativeImagePointer tangentialDiffusionTensorDerivativeImage
      = this->GetTangentialDiffusionTensorDerivativeImage();
  DeformationComponentImageArrayType tangentialDeformationComponentImages
      = this->GetTangentialDeformationComponentImages();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because the input has been copied to the output.

  // Setup the types of structs for the face calculations
  // (Struct handles the case where the image pointer doesn't exist)
  FaceStruct< OutputImagePointer > outputStruct(
      output, regionToProcess, radius );
  NeighborhoodType outputNeighborhood;
  UpdateBufferRegionType updateRegion;

  FaceStruct< NormalVectorImagePointer > normalVectorStruct(
      m_NormalVectorImage, normalVectorRegionToProcess, radius );
  NormalVectorNeighborhoodType normalVectorNeighborhood;

  FaceStruct< DiffusionTensorImagePointer > tangentialTensorStruct(
      tangentialDiffusionTensorImage, tensorRegionToProcess, radius );
  DiffusionTensorNeighborhoodType tangentialTensorNeighborhood;
  FaceStruct< DiffusionTensorImagePointer > normalTensorStruct(
      m_NormalDiffusionTensorImage, tensorRegionToProcess, radius );
  DiffusionTensorNeighborhoodType normalTensorNeighborhood;

  FaceStruct< TensorDerivativeImagePointer > tangentialTensorDerivativeStruct(
      tangentialDiffusionTensorDerivativeImage,
      derivativeRegionToProcess,
      radius );
  TensorDerivativeImageRegionType tangentialTensorDerivativeRegion;
  FaceStruct< TensorDerivativeImagePointer > normalTensorDerivativeStruct(
      m_NormalDiffusionTensorDerivativeImage,
      derivativeRegionToProcess,
      radius );
  TensorDerivativeImageRegionType normalTensorDerivativeRegion;

  typedef FaceStruct< DeformationVectorComponentImagePointer >
      DeformationComponentStructType;
  typedef itk::FixedArray< DeformationComponentStructType, ImageDimension >
      DeformationComponentStructArrayType;
  DeformationComponentStructArrayType tangentialDeformationComponentStructs;
  DeformationComponentStructArrayType normalDeformationComponentStructs;
  DeformationVectorComponentNeighborhoodArrayType
      tangentialDeformationComponentNeighborhoods;
  DeformationVectorComponentNeighborhoodArrayType
      normalDeformationComponentNeighborhoods;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    tangentialDeformationComponentStructs[i] = DeformationComponentStructType(
        tangentialDeformationComponentImages[i],
        deformationComponentRegionToProcess,
        radius );
    normalDeformationComponentStructs[i] = DeformationComponentStructType(
        m_NormalDeformationComponentImages[i],
        deformationComponentRegionToProcess,
        radius );
    }

  // Get the type of registration
  bool computeRegularization = this->GetComputeRegularizationTerm();

  // Go to the first face
  outputStruct.begin();
  if( computeRegularization )
    {
    tangentialTensorStruct.begin();
    tangentialTensorDerivativeStruct.begin();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      tangentialDeformationComponentStructs[i].begin();
      }
    normalVectorStruct.begin();
    normalTensorStruct.begin();
    normalTensorDerivativeStruct.begin();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      normalDeformationComponentStructs[i].begin();
      }
    } // end going to first face

  // Iterate over each face
  while( !outputStruct.IsAtEnd() )
    {

    // Set the neighborhood iterators to the current face
    outputNeighborhood = NeighborhoodType( radius,
                                           output,
                                           *outputStruct.faceListIt );
    updateRegion = UpdateBufferRegionType( updateBuffer,
                                           *outputStruct.faceListIt );
    if( computeRegularization )
      {
      tangentialTensorNeighborhood = DiffusionTensorNeighborhoodType(
          radius,
          tangentialDiffusionTensorImage,
          *tangentialTensorStruct.faceListIt );
      tangentialTensorDerivativeRegion = TensorDerivativeImageRegionType(
          tangentialDiffusionTensorDerivativeImage,
          *tangentialTensorDerivativeStruct.faceListIt );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        tangentialDeformationComponentNeighborhoods[i]
            = DeformationVectorComponentNeighborhoodType(
                radius,
                tangentialDeformationComponentImages[i],
                *tangentialDeformationComponentStructs[i].faceListIt );
        }
      normalVectorNeighborhood = NormalVectorNeighborhoodType(
          radius, m_NormalVectorImage, *normalVectorStruct.faceListIt );
      normalTensorNeighborhood = DiffusionTensorNeighborhoodType(
          radius,
          m_NormalDiffusionTensorImage,
          *normalTensorStruct.faceListIt );
      normalTensorDerivativeRegion = TensorDerivativeImageRegionType(
          m_NormalDiffusionTensorDerivativeImage,
          *normalTensorDerivativeStruct.faceListIt );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        normalDeformationComponentNeighborhoods[i]
            = DeformationVectorComponentNeighborhoodType(
                radius,
                m_NormalDeformationComponentImages[i],
                *normalDeformationComponentStructs[i].faceListIt );
        }
      } // end setting neighborhood iterators to the current face

    // Go to the beginning of the neighborhood for this face
    outputNeighborhood.GoToBegin();
    updateRegion.GoToBegin();
    if( computeRegularization )
      {
      tangentialTensorNeighborhood.GoToBegin();
      tangentialTensorDerivativeRegion.GoToBegin();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        tangentialDeformationComponentNeighborhoods[i].GoToBegin();
        }
      normalVectorNeighborhood.GoToBegin();
      normalTensorNeighborhood.GoToBegin();
      normalTensorDerivativeRegion.GoToBegin();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        normalDeformationComponentNeighborhoods[i].GoToBegin();
        }
      } // end going to the beginning of the neighborhood for this face

    // Iterate through the neighborhood for this face and compute updates
    while( !outputNeighborhood.IsAtEnd() )
      {
      updateRegion.Value() = df->ComputeUpdate(
          outputNeighborhood,
          tangentialTensorNeighborhood,
          tangentialTensorDerivativeRegion,
          tangentialDeformationComponentNeighborhoods,
          normalVectorNeighborhood,
          normalTensorNeighborhood,
          normalTensorDerivativeRegion,
          normalDeformationComponentNeighborhoods,
          spacing,
          globalData );

      // Go to the next neighborhood
      ++outputNeighborhood;
      ++updateRegion;
      if( computeRegularization )
        {
        ++tangentialTensorNeighborhood;
        ++tangentialTensorDerivativeRegion;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          ++tangentialDeformationComponentNeighborhoods[i];
          }
        ++normalVectorNeighborhood;
        ++normalTensorNeighborhood;
        ++normalTensorDerivativeRegion;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          ++normalDeformationComponentNeighborhoods[i];
          }
        } // end going to the next neighborhood
      } // end iterating through the neighborhood for this face

    // Go to the next face
    ++outputStruct.faceListIt;
    if( computeRegularization )
      {
      ++tangentialTensorStruct.faceListIt;
      ++tangentialTensorDerivativeStruct.faceListIt;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        ++tangentialDeformationComponentStructs[i].faceListIt;
        }
      ++normalVectorStruct.faceListIt;
      ++normalTensorStruct.faceListIt;
      ++normalTensorDerivativeStruct.faceListIt;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        ++normalDeformationComponentStructs[i].faceListIt;
        }
      }
    } // end iterating over each face

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  TimeStepType timeStep = df->ComputeGlobalTimeStep(globalData);
  df->ReleaseGlobalDataPointer(globalData);

  return timeStep;
}

} // end namespace itk

#endif
