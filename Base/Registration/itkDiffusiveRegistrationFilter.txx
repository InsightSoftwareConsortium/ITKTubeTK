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
#ifndef __itkDiffusiveRegistrationFilter_txx
#define __itkDiffusiveRegistrationFilter_txx

#include "itkDiffusiveRegistrationFilter.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
DiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::DiffusiveRegistrationFilter()
{
  m_UpdateBuffer = UpdateBufferType::New();

  // Initialize attributes to NULL
  m_DiffusionTensorImage              = 0;
  m_DiffusionTensorDerivativeImage    = 0;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_DeformationComponentImages[i]   = 0;
    }

  // We are using our own regularization, so don't use the implementation
  // provided by the PDERegistration framework
  this->SmoothDeformationFieldOff();
  this->SmoothUpdateFieldOff();
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_DiffusionTensorImage )
    {
    os << indent << "Diffusion tensor image:" << std::endl;
    m_DiffusionTensorImage->Print( os, indent );
    }
  if( m_DiffusionTensorDerivativeImage )
    {
    os << indent << "Diffusion tensor derivative image:"
        << std::endl;
    m_DiffusionTensorDerivativeImage->Print( os, indent );
    }
  if( m_DeformationComponentImages.Length != 0 )
    {
    os << indent << "Deformation component images:" << std::endl;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( m_DeformationComponentImages[i] )
        {
        m_DeformationComponentImages[i]->Print( os, indent );
        }
      }
    }
}

/**
 * Create the registration function
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
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
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::RegistrationFunctionType *
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetRegistrationFunctionPointer() const
{
  RegistrationFunctionType * df = dynamic_cast< RegistrationFunctionType * >
       ( this->GetDifferenceFunction().GetPointer() );
  return df;
}

/**
 * Helper function to allocate space for an image given a template image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
template < class UnallocatedImageType, class TemplateImageType >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateSpaceForImage( UnallocatedImageType& image,
                         const TemplateImageType& templateImage )
{
  assert( image );
  assert( templateImage );
  image->SetOrigin( templateImage->GetOrigin() );
  image->SetSpacing( templateImage->GetSpacing() );
  image->SetDirection( templateImage->GetDirection() );
  image->SetLargestPossibleRegion( templateImage->GetLargestPossibleRegion() );
  image->SetRequestedRegion( templateImage->GetRequestedRegion() );
  image->SetBufferedRegion( templateImage->GetBufferedRegion() );
  image->Allocate();
}

/**
 * Helper function to check whether the attributes of an image matches template
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
template < class CheckedImageType, class TemplateImageType >
bool
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CompareImageAttributes( const CheckedImageType& image,
                          const TemplateImageType& templateImage )
{
  assert( image );
  assert( templateImage );
  return image->GetOrigin() == templateImage->GetOrigin()
      && image->GetSpacing() == templateImage->GetSpacing()
      && image->GetDirection() == templateImage->GetDirection()
      && image->GetLargestPossibleRegion()
          == templateImage->GetLargestPossibleRegion()
      && image->GetRequestedRegion() == templateImage->GetRequestedRegion()
      && image->GetBufferedRegion() == templateImage->GetBufferedRegion();
}

/**
 * Allocate space for the update buffer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateUpdateBuffer()
{
  /* The update buffer looks just like the output and holds the voxel changes */
  typename OutputImageType::Pointer output = this->GetOutput();
  assert( output );
  this->AllocateSpaceForImage( m_UpdateBuffer, output );
}

/**
 * All other initialization done before the initialize / calculate change /
 * apply update loop
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::Initialize()
{
  Superclass::Initialize();

  // Ensure we have a fixed image, moving image and deformation field
  if ( !this->GetFixedImage() || !this->GetMovingImage() )
    {
    itkExceptionMacro( << "Fixed image and/or moving image not set" );
    }

  // Check the timestep for stability if we are using the diffusive or
  // anisotropic diffusive regularization terms
  RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );
  df->CheckTimeStepStability( this->GetInput(), this->GetUseImageSpacing() );

  // Update the registration function's deformation field
  assert( this->GetDeformationField() );
  df->SetDeformationField( this->GetDeformationField() );

  // If we're not computing the regularization term, we're done
  if( !this->GetComputeRegularizationTerm() )
    {
    return;
    }

  // Allocate the images we will use to store data computed during the
  // registration
  this->AllocateImages();

  // Compute the diffusion tensors and their derivatives
  this->ComputeDiffusionTensorImages();
  this->ComputeDiffusionTensorDerivativeImages();
}

/**
 * Allocate the images we will use to store data computed during the
 * registration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateImages()
{
  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the diffusion tensor image and its derivative
  m_DiffusionTensorImage = DiffusionTensorImageType::New();
  this->AllocateSpaceForImage( m_DiffusionTensorImage, output );

  m_DiffusionTensorDerivativeImage = TensorDerivativeImageType::New();
  this->AllocateSpaceForImage( m_DiffusionTensorDerivativeImage, output );
}

/**
 * Update x, y, z components of a deformation field
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ExtractXYZComponentsFromDeformationField(
    OutputImagePointer deformationField,
    DeformationComponentImageArrayType& deformationComponentImages )
{
  typename VectorIndexSelectionFilterType::Pointer indexSelector;
  assert( deformationField );
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    indexSelector = VectorIndexSelectionFilterType::New();
    indexSelector->SetInput( deformationField );
    indexSelector->SetIndex( i );
    deformationComponentImages[i] = indexSelector->GetOutput();
    indexSelector->Update();
    }
}

/**
 * Updates the diffusion tensor image before each run of the registration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( m_DiffusionTensorImage );

  // For the Gaussian regularization, we only need to set the
  // diffusion tensors to the identity
  typename DiffusionTensorImageType::PixelType identityTensor;
  identityTensor.SetIdentity();
  m_DiffusionTensorImage->FillBuffer( identityTensor );
}

/**
 * Updates the diffusion tensor image derivatives before each run of the
 * registration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorDerivativeImages()
{
  assert( this->GetComputeRegularizationTerm() );

  // Compute the diffusion tensor derivative image
  this->ComputeDiffusionTensorDerivativeImageHelper(
      m_DiffusionTensorImage, m_DiffusionTensorDerivativeImage );
}

/**
 * Actually computes the diffusion tensor derivative images
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorDerivativeImageHelper(
    DiffusionTensorImagePointer tensorImage,
    TensorDerivativeImagePointer tensorDerivativeImage )
{
  assert( tensorImage );
  assert( tensorDerivativeImage );

  // Get the FiniteDifferenceFunction to use in calculations.
  const RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );
  const RegularizationFunctionPointer reg
      = df->GetRegularizationFunctionPointer();
  assert( reg );

  // Get the radius
  const typename OutputImageType::SizeType radius = df->GetRadius();

  // Get the spacing
  SpacingType spacing = tensorImage->GetSpacing();

  // Setup the structs for the face calculations, the face iterators, and the
  // iterators over the current face
  FaceStruct< DiffusionTensorImagePointer > tensorStruct( tensorImage, radius );
  DiffusionTensorNeighborhoodType tensorNeighborhood;

  FaceStruct< TensorDerivativeImagePointer > tensorDerivativeStruct(
      tensorDerivativeImage, radius );
  TensorDerivativeImageRegionType tensorDerivativeRegion;

  for(tensorStruct.begin(), tensorDerivativeStruct.begin();
      !tensorDerivativeStruct.IsAtEnd();
      ++tensorStruct.faceListIt, ++tensorDerivativeStruct.faceListIt )
    {
    // Set the neighborhood iterators to the current face
    tensorNeighborhood = DiffusionTensorNeighborhoodType(
        radius, tensorImage, *tensorStruct.faceListIt );
    tensorDerivativeRegion = TensorDerivativeImageRegionType(
        tensorDerivativeImage, *tensorDerivativeStruct.faceListIt );

    // Iterate through the neighborhood for this face and compute derivatives
    for( tensorNeighborhood.GoToBegin(), tensorDerivativeRegion.GoToBegin();
         !tensorNeighborhood.IsAtEnd();
         ++tensorNeighborhood, ++tensorDerivativeRegion )
           {
      reg->ComputeDiffusionTensorFirstOrderPartialDerivatives(
          tensorNeighborhood, tensorDerivativeRegion, spacing );
      }
    }
}

/**
 * Initialize the state of the filter and equation before each iteration.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeIteration()
{
  Superclass::InitializeIteration();

  // Update the deformation field component images
  // Since the components depend on the current deformation field, they must be
  // computed on every registration iteration
  if( this->GetComputeRegularizationTerm() )
    {
    this->UpdateDeformationVectorComponentImages();
    }
}

/**
 * Updates the deformation vector component images before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationVectorComponentImages()
{
  assert( this->GetComputeRegularizationTerm() );

  // Update the extracted components
  this->ExtractXYZComponentsFromDeformationField(
      this->GetOutput(), m_DeformationComponentImages );
}

/**
 * Populates the update buffer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
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

  // Explicitly call Modified on m_UpdateBuffer here, since
  // ThreadedCalculateChange changes this buffer through iterators which do not
  // increment the update buffer timestamp
  this->m_UpdateBuffer->Modified();

  return dt;
}

/**
 * Calls ThreadedCalculateChange for processing
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_TYPE
DiffusiveRegistrationFilter
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
                                             splitRegion );;

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
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChange( const ThreadRegionType &, int)
{
  // This function should never be called!
  itkExceptionMacro( << "ThreadedCalculateChange(regionToProcess, threadId) "
                     << "should never be called.  Use the other "
                     << "ThreadedCalculateChange function instead" );
}

/**
 * Does the actual work of calculating change over a region supplied by the
 * multithreading mechanism
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChange(
    const ThreadRegionType &                      regionToProcess,
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

  // Get the radius and output
  const typename OutputImageType::SizeType radius = df->GetRadius();
  OutputImagePointer output = this->GetOutput();
  SpacingType spacing = output->GetSpacing();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because the input has been copied to the output.

  // Setup the types of structs for the face calculations
  // (Struct handles the case where the image pointer doesn't exist)
  FaceStruct< OutputImagePointer > outputStruct(
      output, regionToProcess, radius );
  NeighborhoodType outputNeighborhood;
  UpdateBufferRegionType updateRegion;

  FaceStruct< DiffusionTensorImagePointer > tensorStruct(
      m_DiffusionTensorImage, tensorRegionToProcess, radius );
  DiffusionTensorNeighborhoodType tensorNeighborhood;

  FaceStruct< TensorDerivativeImagePointer > tensorDerivativeStruct(
      m_DiffusionTensorDerivativeImage,
      derivativeRegionToProcess,
      radius );
  TensorDerivativeImageRegionType tensorDerivativeRegion;

  typedef FaceStruct< DeformationVectorComponentImagePointer >
      DeformationComponentStructType;
  typedef itk::FixedArray< DeformationComponentStructType, ImageDimension >
      DeformationComponentStructArrayType;
  DeformationComponentStructArrayType deformationComponentStructs;
  DeformationVectorComponentNeighborhoodArrayType
      deformationComponentNeighborhoods;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    deformationComponentStructs[i] = DeformationComponentStructType(
        m_DeformationComponentImages[i],
        deformationComponentRegionToProcess,
        radius );
    }

  // Get the type of registration
  bool computeRegularization = this->GetComputeRegularizationTerm();

  // Go to the first face
  outputStruct.begin();
  if( computeRegularization )
    {
    tensorStruct.begin();
    tensorDerivativeStruct.begin();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      deformationComponentStructs[i].begin();
      }
    } // end going to first face

  // Iterate over each face
  while( !outputStruct.IsAtEnd() )
    {

    // Set the neighborhood iterators to the current face
    outputNeighborhood = NeighborhoodType( radius,
                                           output,
                                           *outputStruct.faceListIt );
    updateRegion = UpdateBufferRegionType( m_UpdateBuffer,
                                           *outputStruct.faceListIt );
    if( computeRegularization )
      {
      tensorNeighborhood = DiffusionTensorNeighborhoodType(
          radius, m_DiffusionTensorImage, *tensorStruct.faceListIt );
      tensorDerivativeRegion = TensorDerivativeImageRegionType(
          m_DiffusionTensorDerivativeImage,
          *tensorDerivativeStruct.faceListIt );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationComponentNeighborhoods[i]
            = DeformationVectorComponentNeighborhoodType(
                radius,
                m_DeformationComponentImages[i],
                *deformationComponentStructs[i].faceListIt );
        }
      } // end setting neighborhood iterators to the current face

    // Go to the beginning of the neighborhood for this face
    outputNeighborhood.GoToBegin();
    updateRegion.GoToBegin();
    if( computeRegularization )
      {
      tensorNeighborhood.GoToBegin();
      tensorDerivativeRegion.GoToBegin();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationComponentNeighborhoods[i].GoToBegin();
        }
      } // end going to the beginning of the neighborhood for this face

    // Iterate through the neighborhood for this face and compute updates
    while( !outputNeighborhood.IsAtEnd() )
      {
      updateRegion.Value() = df->ComputeUpdate(
          outputNeighborhood,
          tensorNeighborhood,
          tensorDerivativeRegion,
          deformationComponentNeighborhoods,
          spacing,
          globalData );

      // Go to the next neighborhood
      ++outputNeighborhood;
      ++updateRegion;
      if( computeRegularization )
        {
        ++tensorNeighborhood;
        ++tensorDerivativeRegion;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          ++deformationComponentNeighborhoods[i];
          }
        } // end going to the next neighborhood
      } // end iterating through the neighborhood for this face

    // Go to the next face
    ++outputStruct.faceListIt;
    if( computeRegularization )
      {
      ++tensorStruct.faceListIt;
      ++tensorDerivativeStruct.faceListIt;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        ++deformationComponentStructs[i].faceListIt;
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

/**
 * Applies changes from the update buffer to the output
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdate(TimeStepType dt)
{
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(
    this->ApplyUpdateThreaderCallback, & str);

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Explicitely call Modified on GetOutput here, since ThreadedApplyUpdate
  // changes this buffer through iterators which do not increment the
  // output timestamp
  this->GetOutput()->Modified();
}

/**
 * Calls ThreadedApplyUpdate, need to reimplement here to also split the
 * diffusion tensor image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_TYPE
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdateThreaderCallback( void * arg )
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

  if (threadId < total)
    {
    str->Filter->ThreadedApplyUpdate(str->TimeStep, splitRegion, threadId );
    }

  return ITK_THREAD_RETURN_VALUE;
}

/**
 * Does the actual work of updating the output from the UpdateContainer
 * over an output region supplied by the multithreading mechanism.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::ThreadedApplyUpdate(TimeStepType dt,
  const ThreadRegionType &regionToProcess, int )
{
  UpdateBufferRegionType  u( m_UpdateBuffer, regionToProcess );
  OutputImageRegionType   o( this->GetOutput(), regionToProcess );

  for( u = u.Begin(), o = o.Begin();
       !u.IsAtEnd();
       ++o, ++u )
         {
    o.Value() += static_cast< DeformationVectorType >( u.Value() * dt );
    // no adaptor support here
    }
}

} // end namespace itk

#endif
