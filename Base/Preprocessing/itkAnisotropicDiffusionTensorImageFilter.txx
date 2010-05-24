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
#ifndef __itkAnisotropicDiffusionTensorImageFilter_txx
#define __itkAnisotropicDiffusionTensorImageFilter_txx

#include "itkAnisotropicDiffusionTensorImageFilter.h"
#include "itkAnisotropicDiffusionTensorFunction.h"

#include <list>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkNeighborhoodAlgorithm.h"

#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "itkFixedArray.h"

//#define INTERMEDIATE_OUTPUTS

namespace itk {

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::AnisotropicDiffusionTensorImageFilter()
{
  m_UpdateBuffer = UpdateBufferType::New(); 

  m_DiffusionTensorImage  = DiffusionTensorImageType::New();
 
  this->SetNumberOfIterations(1);

  m_TimeStep = 0.11; 

  //set the function
  typename AnisotropicDiffusionTensorFunction<UpdateBufferType>::Pointer q
      = AnisotropicDiffusionTensorFunction<UpdateBufferType>::New();
  this->SetDifferenceFunction(q);

}

/** Prepare for the iteration process. */
 template <class TInputImage, class TOutputImage>
 void
 AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
 ::InitializeIteration()
{
  //itkDebugMacro( << "InitializeIteration() called " );
  std::cerr << "InitalizeIteration" << std::endl;

  AnisotropicDiffusionTensorFunction<UpdateBufferType> *f = 
     dynamic_cast<AnisotropicDiffusionTensorFunction<UpdateBufferType> *>
     (this->GetDifferenceFunction().GetPointer());

  if (! f)
    {
    throw ExceptionObject(__FILE__, __LINE__, 
        "Anisotropic diffusion Vessel Enhancement function is not set.",
         ITK_LOCATION);
    }
   
  f->SetTimeStep(m_TimeStep);
   
  // Check the timestep for stability
  double minSpacing;
  if (this->GetUseImageSpacing())
    {
    minSpacing = this->GetInput()->GetSpacing()[0];
    for (unsigned int i = 1; i < ImageDimension; i++)
      {
      if (this->GetInput()->GetSpacing()[i] < minSpacing)
        {
        minSpacing = this->GetInput()->GetSpacing()[i];
        }
      }
    }
  else
    {
    minSpacing = 1.0;
    }

  double ratio = 
     minSpacing /vcl_pow(2.0, static_cast<double>(ImageDimension) + 1);

  if ( m_TimeStep > ratio ) 
    {
    itkWarningMacro(<< std::endl << "Anisotropic diffusion unstable time step:" 
                    << m_TimeStep << std::endl << "Minimum stable time step" 
                    << "for this image is " 
                    << ratio ); 
    }
   
  f->InitializeIteration();
   
  if (this->GetNumberOfIterations() != 0)
    {
    this->UpdateProgress(((float)(this->GetElapsedIterations()))
                          /((float)(this->GetNumberOfIterations())));
    }
  else
    {
    this->UpdateProgress(0);
    }

 //Update the Diffusion tensor image
  this->UpdateDiffusionTensorImage();
}

template <class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::CopyInputToOutput()
{
  typename TInputImage::ConstPointer  input  = this->GetInput();
  typename TOutputImage::Pointer      output = this->GetOutput();

  if ( !input || !output )
    {
    itkExceptionMacro(<< "Either input and/or output is NULL.");
    }

  // Check if we are doing in-place filtering
  if ( this->GetInPlace() && (typeid(TInputImage) == typeid(TOutputImage)) )
    {
    typename TInputImage::Pointer tempPtr = 
      dynamic_cast<TInputImage *>( output.GetPointer() );
    if ( tempPtr && tempPtr->GetPixelContainer() == input->GetPixelContainer() )
      {
      // the input and output container are the same - no need to copy
      return;
      }
    }
  
  ImageRegionConstIterator<TInputImage>  in(input, output->GetRequestedRegion());
  ImageRegionIterator<TOutputImage> out(output, output->GetRequestedRegion());

  while( ! out.IsAtEnd() )
    {
    out.Value() =  static_cast<PixelType>(in.Get());  
    ++in;
    ++out;
    }
}

template <class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::AllocateUpdateBuffer()
{
  itkDebugMacro( << "AllocateUpdateBuffer() called" ); 

  std::cerr << "AllocateUpdaeBuffer() " << std::endl;

  /* The update buffer looks just like the output and holds the change in 
   the pixel  */
  
  typename TOutputImage::Pointer output = this->GetOutput();

  m_UpdateBuffer->SetSpacing(output->GetSpacing());
  m_UpdateBuffer->SetOrigin(output->GetOrigin());
  m_UpdateBuffer->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_UpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_UpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_UpdateBuffer->Allocate();
}

template <class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::AllocateDiffusionTensorImage()
{
  itkDebugMacro( << "AllocateDiffusionTensorImage() called" ); 
  
  std::cerr << "AllocateDiffusionTensorImage() " << std::endl;
  

  /* The diffusionTensor image has the same size as the output and holds 
     the diffusion tensor matrix for each pixel */

  typename TOutputImage::Pointer output = this->GetOutput();

  m_DiffusionTensorImage->SetSpacing(output->GetSpacing());
  m_DiffusionTensorImage->SetOrigin(output->GetOrigin());
  m_DiffusionTensorImage->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_DiffusionTensorImage->SetRequestedRegion(output->GetRequestedRegion());
  m_DiffusionTensorImage->SetBufferedRegion(output->GetBufferedRegion());
  m_DiffusionTensorImage->Allocate();
}

template<class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ApplyUpdate(TimeStepType dt)
{
  itkDebugMacro( << "ApplyUpdate Invoked with time step size: " << dt ); 
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(this->ApplyUpdateThreaderCallback,
                                            &str);
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

template<class TInputImage, class TOutputImage>
ITK_THREAD_RETURN_TYPE
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ApplyUpdateThreaderCallback( void * arg )
{
  DenseFDThreadStruct * str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (DenseFDThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  ThreadDiffusionTensorImageRegionType    splitRegionDiffusionTensorImage;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegionDiffusionTensorImage);
  if (threadId < total)
    {
    str->Filter->ThreadedApplyUpdate(str->TimeStep,
                                     splitRegion, 
                                     splitRegionDiffusionTensorImage,
                                     threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TInputImage, class TOutputImage>
typename
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>::TimeStepType
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::CalculateChange()
{
  itkDebugMacro( << "CalculateChange called" );

  std::cerr << "CalculateChange called" << std::endl;

  int threadCount;
  TimeStepType dt;

  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = NumericTraits<TimeStepType>::Zero;  // Not used during the
  // calculate change step.
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(this->CalculateChangeThreaderCallback,
                                            &str);

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  threadCount = this->GetMultiThreader()->GetNumberOfThreads();  
  str.TimeStepList = new TimeStepType[threadCount];
  str.ValidTimeStepList = new bool[threadCount];
  for (int i =0; i < threadCount; ++i)
    {      str.ValidTimeStepList[i] = false;    } 

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Resolve the single value time step to return
  dt = this->ResolveTimeStep(str.TimeStepList, str.ValidTimeStepList, threadCount);
  delete [] str.TimeStepList;
  delete [] str.ValidTimeStepList;

  return  dt;
}

template <class TInputImage, class TOutputImage>
ITK_THREAD_RETURN_TYPE
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::CalculateChangeThreaderCallback( void * arg )
{
  DenseFDThreadStruct * str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (DenseFDThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;

  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  ThreadDiffusionTensorImageRegionType splitDiffusionimageRegion;

  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitDiffusionimageRegion);

  if (threadId < total)
    { 
    str->TimeStepList[threadId]
      = str->Filter->ThreadedCalculateChange(splitRegion, splitDiffusionimageRegion, threadId);
    str->ValidTimeStepList[threadId] = true;
    }

  return ITK_THREAD_RETURN_VALUE;  
}

template <class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType &regionToProcess,
                      const ThreadDiffusionTensorImageRegionType &,
                      int)
{
  ImageRegionIterator<UpdateBufferType> u(m_UpdateBuffer,    regionToProcess);
  ImageRegionIterator<OutputImageType>  o(this->GetOutput(), regionToProcess);

  u = u.Begin();
  o = o.Begin();

  while ( !u.IsAtEnd() )
    {

    o.Value() += static_cast<PixelType>(u.Value() * dt);  // no adaptor support here

    ++o;
    ++u;
    }
}

template <class TInputImage, class TOutputImage>
typename
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>::TimeStepType
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::ThreadedCalculateChange(const ThreadRegionType &regionToProcess, 
    const ThreadDiffusionTensorImageRegionType & diffusionRegionToProcess, int)
{
  typedef typename OutputImageType::RegionType      RegionType;
  typedef typename OutputImageType::SizeType        SizeType;
  typedef typename OutputImageType::SizeValueType   SizeValueType;
  typedef typename OutputImageType::IndexType       IndexType;
  typedef typename OutputImageType::IndexValueType  IndexValueType;

  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
                                           NeighborhoodIteratorType;
  
  typedef ImageRegionIterator<UpdateBufferType> UpdateIteratorType;

  typename OutputImageType::Pointer output = this->GetOutput();
  TimeStepType timeStep;
  void *globalData;

  // Get the FiniteDifferenceFunction to use in calculations.
  const typename FiniteDifferenceFunctionType::Pointer df = 
     dynamic_cast<AnisotropicDiffusionTensorFunction<UpdateBufferType> *>
     ( this->GetDifferenceFunction().GetPointer());

  const SizeType  radius = df->GetRadius();
  
  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because input has been copied to output.
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<OutputImageType>
                                                        FaceCalculatorType;

  typedef typename FaceCalculatorType::FaceListType FaceListType;

  FaceCalculatorType faceCalculator;
    
  FaceListType faceList = faceCalculator(output, regionToProcess, radius);
  typename FaceListType::iterator fIt = faceList.begin();


  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<DiffusionTensorImageType>
    DiffusionTensorFaceCalculatorType;

  typedef typename DiffusionTensorFaceCalculatorType::FaceListType 
                                                DiffusionTensorFaceListType;

  DiffusionTensorFaceCalculatorType diffusionTensorFaceCalculator;
  
  DiffusionTensorFaceListType diffusionTensorFaceList = 
     diffusionTensorFaceCalculator(m_DiffusionTensorImage, diffusionRegionToProcess, radius);

  typename DiffusionTensorFaceListType::iterator dfIt = diffusionTensorFaceList.begin();
  

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();

  // Process the non-boundary region.
  NeighborhoodIteratorType nD(radius, output, *fIt);
  UpdateIteratorType       nU(m_UpdateBuffer,  *fIt);
  DiffusionTensorNeighborhoodType  dTN(radius,m_DiffusionTensorImage, *dfIt); 

  nD.GoToBegin();
  nU.GoToBegin();
  dTN.GoToBegin();
  while( !nD.IsAtEnd() )
    {
    nU.Value() = df->ComputeUpdate(nD, dTN, globalData);
    ++nD;
    ++nU;
    ++dTN;
    }

  // Process each of the boundary faces.
  NeighborhoodIteratorType bD;
  DiffusionTensorNeighborhoodType bDD;
  UpdateIteratorType   bU;
  ++dfIt;
  for (++fIt; fIt != faceList.end(); ++fIt)
    {
    bD = NeighborhoodIteratorType(radius, output, *fIt);
    bDD = DiffusionTensorNeighborhoodType(radius, m_DiffusionTensorImage, *dfIt);
    bU = UpdateIteratorType  (m_UpdateBuffer, *fIt);
     
    bD.GoToBegin();
    bU.GoToBegin();
    bDD.GoToBegin();
    while ( !bD.IsAtEnd() )
      {
      bU.Value() = df->ComputeUpdate(bD,bDD,globalData);
      ++bD;
      ++bU;
      ++bDD;
      }
    ++dfIt;
    }

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  timeStep = df->ComputeGlobalTimeStep(globalData);
  df->ReleaseGlobalDataPointer(globalData);

  return timeStep;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  itkDebugMacro( << "GenerateData is called" );

  if (this->GetState() == Superclass::UNINITIALIZED)
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

  while ( ! this->Halt() )
    {
    std::cout << "Iteration:\t" << iter << std::endl;
    this->InitializeIteration(); // An optional method for precalculating
                                 // global values, or otherwise setting up
                                 // for the next iteration
    dt = this->CalculateChange();

    this->ApplyUpdate(dt);

    ++iter;

    this->SetElapsedIterations( iter );

    // Invoke the iteration event.
    this->InvokeEvent( IterationEvent() );
    if( this->GetAbortGenerateData() )
      {
      this->InvokeEvent( IterationEvent() );
      this->ResetPipeline(); 
      throw ProcessAborted(__FILE__,__LINE__);
      }
    }
} 
 
template <class TInputImage, class TOutputImage>
typename AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::DiffusionTensorImagePointerType
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::GetDiffusionTensorImage()
{
  return m_DiffusionTensorImage;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "TimeStep: " << m_TimeStep  << std::endl;
}

}// end namespace itk

#endif
