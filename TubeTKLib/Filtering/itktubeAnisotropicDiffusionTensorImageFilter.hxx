/*=========================================================================

Library:   TubeTKLib

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

#ifndef __itktubeAnisotropicDiffusionTensorImageFilter_hxx
#define __itktubeAnisotropicDiffusionTensorImageFilter_hxx

#include "itktubeAnisotropicDiffusionTensorImageFilter.h"

#include "itktubeAnisotropicDiffusionTensorFunction.h"

#include <itkFixedArray.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkNumericTraits.h>
#include <itkVector.h>

#include <list>

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TInputImage, class TOutputImage >
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::AnisotropicDiffusionTensorImageFilter( void )
{
  m_UpdateBuffer = UpdateBufferType::New();
  m_DiffusionTensorImage = DiffusionTensorImageType::New();
  this->SetNumberOfIterations( 1 );
  m_TimeStep = 0.11;

  //set the finite difference function object
  typename AnisotropicDiffusionTensorFunction<UpdateBufferType>::Pointer q
      = AnisotropicDiffusionTensorFunction<UpdateBufferType>::New();
  this->SetDifferenceFunction( q );
}

/** Prepare for the iteration process - called at the beginning of each
 * iteration.
 */
 template< class TInputImage, class TOutputImage >
 void
 AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
 ::InitializeIteration( void )
{
  itkDebugMacro( << "InitializeIteration() called ." );

  AnisotropicDiffusionTensorFunction<UpdateBufferType> *f =
     dynamic_cast<AnisotropicDiffusionTensorFunction<UpdateBufferType> *>
     ( this->GetDifferenceFunction().GetPointer() );

  if( !f )
    {
    throw ExceptionObject( __FILE__, __LINE__,
        "Anisotropic diffusion Vessel Enhancement function is not set.",
         ITK_LOCATION );
    }

  f->SetTimeStep( m_TimeStep );

  // Check the timestep for stability
  f->CheckTimeStepStability( this->GetInput(), this->GetUseImageSpacing() );

  f->InitializeIteration();

  if( this->GetNumberOfIterations() != 0 )
    {
    this->UpdateProgress( ( ( float )( this->GetElapsedIterations() ) )
                          /( ( float )( this->GetNumberOfIterations() ) ) );
    }
  else
    {
    this->UpdateProgress( 0 );
    }

  // Update the diffusion tensor image: implemented in subclasses, for example
  // to calculate the structure tensor and its eigenvectors and eigenvalues
  this->UpdateDiffusionTensorImage();
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::CopyInputToOutput( void )
{
  typename TInputImage::ConstPointer  input  = this->GetInput();
  typename TOutputImage::Pointer      output = this->GetOutput();

  if( !input || !output )
    {
    itkExceptionMacro( << "Either input and/or output is NULL." );
    }

  // Check if we are doing in-place filtering
  if( this->GetInPlace() && ( typeid( TInputImage ) ==
    typeid( TOutputImage ) ) )
    {
    typename TInputImage::Pointer tempPtr =
      dynamic_cast<TInputImage *>( output.GetPointer() );
    if( tempPtr && tempPtr->GetPixelContainer() ==
      input->GetPixelContainer() )
      {
      // the input and output container are the same - no need to copy
      return;
      }
    }

  ImageRegionConstIterator<TInputImage> in( input,
    output->GetRequestedRegion() );
  ImageRegionIterator<TOutputImage> out( output,
    output->GetRequestedRegion() );

  while( !out.IsAtEnd() )
    {
    out.Value() = static_cast<PixelType>( in.Get() );
    ++in;
    ++out;
    }
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::AllocateUpdateBuffer( void )
{
  itkDebugMacro( << "AllocateUpdateBuffer() called." );

  /* The update buffer looks just like the output and holds the change in
   the pixel  */

  typename TOutputImage::Pointer output = this->GetOutput();

  m_UpdateBuffer->SetSpacing( output->GetSpacing() );
  m_UpdateBuffer->SetOrigin( output->GetOrigin() );
  m_UpdateBuffer->SetDirection( output->GetDirection() );
  m_UpdateBuffer->SetLargestPossibleRegion(
    output->GetLargestPossibleRegion() );
  m_UpdateBuffer->SetRequestedRegion( output->GetRequestedRegion() );
  m_UpdateBuffer->SetBufferedRegion( output->GetBufferedRegion() );
  m_UpdateBuffer->Allocate();
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::AllocateDiffusionTensorImage( void )
{
  itkDebugMacro( << "AllocateDiffusionTensorImage() called." );

  /* The diffusionTensor image has the same size as the output and holds
     the diffusion tensor matrix for each pixel */

  typename TOutputImage::Pointer output = this->GetOutput();

  m_DiffusionTensorImage->SetSpacing( output->GetSpacing() );
  m_DiffusionTensorImage->SetOrigin( output->GetOrigin() );
  m_DiffusionTensorImage->SetDirection( output->GetDirection() );
  m_DiffusionTensorImage->SetLargestPossibleRegion(
    output->GetLargestPossibleRegion() );
  m_DiffusionTensorImage->SetRequestedRegion(
    output->GetRequestedRegion() );
  m_DiffusionTensorImage->SetBufferedRegion( output->GetBufferedRegion() );
  m_DiffusionTensorImage->Allocate();
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ApplyUpdate( const TimeStepType & dt )
{
  itkDebugMacro( << "ApplyUpdate Invoked with time step size: " << dt );
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfThreads(
    this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(
    this->ApplyUpdateThreaderCallback, &str );
  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

#ifdef INTERMEDIATE_OUTPUTS
  typedef ImageFileWriter< OutputImageType > WriterType;
  typename WriterType::Pointer   writer = WriterType::New();
  writer->SetFileName( "UpdatedOutputImage.mha" );
  writer->SetInput( m_UpdateBuffer );
  writer->Update();
#endif
}

template< class TInputImage, class TOutputImage >
ITK_THREAD_RETURN_TYPE
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ApplyUpdateThreaderCallback( void * arg )
{
  DenseFDThreadStruct * str;
  ThreadIdType total, threadId, threadCount;

  threadId = ( ( MultiThreaderBase::ThreadInfoStruct * )( arg ) )->ThreadID;
  threadCount = ( ( MultiThreaderBase::ThreadInfoStruct * )( arg ) )->
    NumberOfThreads;

  str = ( DenseFDThreadStruct * )( ( ( MultiThreaderBase::ThreadInfoStruct * )
    ( arg ) )->UserData );

  // Execute the actual method with appropriate output region.
  // First find out how many pieces the extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount, splitRegion );

  ThreadDiffusionTensorImageRegionType splitRegionDiffusionTensorImage;
  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitRegionDiffusionTensorImage );
  if( threadId < total )
    {
    str->Filter->ThreadedApplyUpdate( str->TimeStep, splitRegion,
      splitRegionDiffusionTensorImage, threadId );
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template< class TInputImage, class TOutputImage >
typename
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>::TimeStepType
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::CalculateChange( void )
{
  itkDebugMacro( << "CalculateChange called." );

  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = NumericTraits<TimeStepType>::Zero;  // Not used during the
  // calculate change step.
  this->GetMultiThreader()->SetNumberOfThreads(
    this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(
    this->CalculateChangeThreaderCallback, &str );

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  ThreadIdType threadCount = this->GetMultiThreader()->GetNumberOfThreads();
  str.TimeStepList.resize( threadCount );
  str.ValidTimeStepList.resize( threadCount );
  for( ThreadIdType i =0; i < threadCount; ++i )
    {
    str.ValidTimeStepList[i] = false;
    }

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Resolve the single value time step to return
  TimeStepType dt = this->ResolveTimeStep( str.TimeStepList,
    str.ValidTimeStepList );

  return  dt;
}

template< class TInputImage, class TOutputImage >
ITK_THREAD_RETURN_TYPE
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::CalculateChangeThreaderCallback( void * arg )
{
  DenseFDThreadStruct * str;
  ThreadIdType total, threadId, threadCount;

  threadId = ( ( MultiThreaderBase::ThreadInfoStruct * )( arg ) )->ThreadID;
  threadCount = ( ( MultiThreaderBase::ThreadInfoStruct * )( arg ) )->
    NumberOfThreads;

  str = ( DenseFDThreadStruct * )( ( ( MultiThreaderBase::ThreadInfoStruct * )
    ( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;

  str->Filter->SplitRequestedRegion( threadId, threadCount, splitRegion );

  ThreadDiffusionTensorImageRegionType splitDiffusionimageRegion;

  total = str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitDiffusionimageRegion );

  if( threadId < total )
    {
    str->TimeStepList[threadId]
      = str->Filter->ThreadedCalculateChange( splitRegion,
        splitDiffusionimageRegion, threadId );
    str->ValidTimeStepList[threadId] = true;
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ThreadedApplyUpdate( TimeStepType dt,
  const ThreadRegionType &regionToProcess,
  const ThreadDiffusionTensorImageRegionType &, ThreadIdType )
{
  ImageRegionIterator<UpdateBufferType> u( m_UpdateBuffer,
    regionToProcess );
  ImageRegionIterator<OutputImageType>  o( this->GetOutput(),
    regionToProcess );

  u.GoToBegin();
  o.GoToBegin();

  while( !u.IsAtEnd() )
    {

    // no adaptor support
    o.Value() += static_cast<PixelType>( u.Value() * dt );

    ++o;
    ++u;
    }
}

template< class TInputImage, class TOutputImage >
typename AnisotropicDiffusionTensorImageFilter<TInputImage,
  TOutputImage>::TimeStepType
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ThreadedCalculateChange( const ThreadRegionType &regionToProcess,
  const ThreadDiffusionTensorImageRegionType & diffusionRegionToProcess,
  ThreadIdType )
{
  typedef typename OutputImageType::SizeType
    SizeType;
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
    NeighborhoodIteratorType;
  typedef ImageRegionIterator<UpdateBufferType>
    UpdateIteratorType;

  typename OutputImageType::Pointer output = this->GetOutput();
  typename FiniteDifferenceFunctionType::SpacingType spacing =
    output->GetSpacing();

  TimeStepType timeStep;
  void *globalData;

  // Get the FiniteDifferenceFunction to use in calculations.
  const typename FiniteDifferenceFunctionType::Pointer df =
   dynamic_cast<AnisotropicDiffusionTensorFunction<UpdateBufferType> *>(
   this->GetDifferenceFunction().GetPointer() );

  const SizeType  radius = df->GetRadius();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because input has been copied to output.
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
    OutputImageType >   FaceCalculatorType;

  typedef typename FaceCalculatorType::FaceListType FaceListType;

  FaceCalculatorType faceCalculator;

  FaceListType faceList = faceCalculator( output, regionToProcess, radius );
  typename FaceListType::iterator fIt = faceList.begin();

  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
    DiffusionTensorImageType>
    DiffusionTensorFaceCalculatorType;

  typedef typename DiffusionTensorFaceCalculatorType::FaceListType
    DiffusionTensorFaceListType;

  DiffusionTensorFaceCalculatorType diffusionTensorFaceCalculator;

  DiffusionTensorFaceListType diffusionTensorFaceList =
   diffusionTensorFaceCalculator( m_DiffusionTensorImage,
   diffusionRegionToProcess, radius );

  typename DiffusionTensorFaceListType::iterator dfIt =
    diffusionTensorFaceList.begin();

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();

  // Process the non-boundary region.
  NeighborhoodIteratorType nD( radius, output, *fIt );
  UpdateIteratorType       nU( m_UpdateBuffer, *fIt );
  DiffusionTensorNeighborhoodType  dTN( radius, m_DiffusionTensorImage,
    *dfIt );

  nD.GoToBegin();
  nU.GoToBegin();
  dTN.GoToBegin();
  while( !nD.IsAtEnd() )
    {
    nU.Value() = df->ComputeUpdate( nD, dTN, spacing, globalData );
    ++nD;
    ++nU;
    ++dTN;
    }

  // Process each of the boundary faces.
  NeighborhoodIteratorType bD;
  DiffusionTensorNeighborhoodType bDD;
  UpdateIteratorType   bU;
  ++dfIt;
  for( ++fIt; fIt != faceList.end(); ++fIt )
    {
    bD = NeighborhoodIteratorType( radius, output, *fIt );
    bDD = DiffusionTensorNeighborhoodType( radius, m_DiffusionTensorImage,
      *dfIt );
    bU = UpdateIteratorType  ( m_UpdateBuffer, *fIt );

    bD.GoToBegin();
    bU.GoToBegin();
    bDD.GoToBegin();
    while( !bD.IsAtEnd() )
      {
      bU.Value() = df->ComputeUpdate( bD, bDD, spacing, globalData );
      ++bD;
      ++bU;
      ++bDD;
      }
    ++dfIt;
    }

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  timeStep = df->ComputeGlobalTimeStep( globalData );
  df->ReleaseGlobalDataPointer( globalData );

  return timeStep;
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::GenerateData( void )
{
  itkDebugMacro( << "GenerateData is called." );

  if( !this->GetIsInitialized() )
    {
    // Allocate the output image
    this->AllocateOutputs();

    // Copy the input image to the output image.  Algorithms will operate
    // directly on the output image and the update buffer.
    this->CopyInputToOutput();

    // Allocate the internal update buffer.
    this->AllocateUpdateBuffer();

    // Allocate buffer for the diffusion tensor image
    this->AllocateDiffusionTensorImage();

    this->SetStateToInitialized();

    this->SetElapsedIterations( 0 );
    }

  // Iterative algorithm
  TimeStepType dt;
  unsigned int iter = 0;

  while( !this->Halt() )
    {
    this->InitializeIteration(); // An optional method for precalculating
                                 // global values, or otherwise setting up
                                 // for the next iteration
    dt = this->CalculateChange();

    this->ApplyUpdate( dt );

    ++iter;

    this->SetElapsedIterations( iter );

    // Invoke the iteration event.
    this->InvokeEvent( IterationEvent() );
    if( this->GetAbortGenerateData() )
      {
      this->InvokeEvent( IterationEvent() );
      this->ResetPipeline();
      throw ProcessAborted( __FILE__, __LINE__ );
      }
    }
}

template< class TInputImage, class TOutputImage >
typename AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::DiffusionTensorImagePointerType
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::GetDiffusionTensorImage( void )
{
  return m_DiffusionTensorImage;
}

template< class TInputImage, class TOutputImage >
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "TimeStep: " << m_TimeStep  << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeAnisotropicDiffusionTensorImageFilter_hxx )
