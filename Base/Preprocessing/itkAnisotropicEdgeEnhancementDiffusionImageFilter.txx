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
#ifndef __itkAnisotropicEdgeEnhancementDiffusionImageFilter_txx
#define __itkAnisotropicEdgeEnhancementDiffusionImageFilter_txx

#include "itkAnisotropicEdgeEnhancementDiffusionImageFilter.h"
#include "itkAnisotropicEdgeEnhancementDiffusionFunction.h"

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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
::AnisotropicEdgeEnhancementDiffusionImageFilter()
{
  m_UpdateBuffer = UpdateBufferType::New(); 

  m_DiffusionTensorImage  = DiffusionTensorImageType::New();
 
  this->SetNumberOfIterations(1);

  m_TimeStep = 10e-3;

  //set the function
  typename AnisotropicEdgeEnhancementDiffusionFunction<UpdateBufferType>::Pointer q
      = AnisotropicEdgeEnhancementDiffusionFunction<UpdateBufferType>::New();
  this->SetDifferenceFunction(q);

  // Vesselness guided vesselness function algorithm parameter
  m_WStrength  = 25.0;
  m_Sensitivity  = 5.0;
  m_Epsilon = 10e-2;
}

/** Prepare for the iteration process. */
 template <class TInputImage, class TOutputImage>
 void
 AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
 ::InitializeIteration()
{
  //itkDebugMacro( << "InitializeIteration() called " );
  std::cerr << "InitalizeIteration" << std::endl;

  AnisotropicEdgeEnhancementDiffusionFunction<UpdateBufferType> *f = 
     dynamic_cast<AnisotropicEdgeEnhancementDiffusionFunction<UpdateBufferType> *>
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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

template <class TInputImage, class TOutputImage>
void
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
::UpdateDiffusionTensorImage()
{
  itkDebugMacro( << "UpdateDiffusionTensorImage() called" ); 

  std::cerr << "UpdateDiffusionTensorImage()" << std::endl;

  /* IN THIS METHOD, the following items will be implemented
   - Compute the structure tensor ( Multiscale version structure tensor )
   - Compute its eigen vectors
   - Compute eigen values corresponding to the diffusion matrix tensor 
     ( Here is where all the magic happens for EED, CED and hybrid switch )
  */

  //Step 1: Compute the structure tensor and identify the eigen vectors 
  //Step 1.1: Compute the structure tensor
  m_StructureTensorFilter->SetInput( this->GetOutput() );
  m_StructureTensorFilter->Update();

  // Step 1.2: Identify the eigen vectors of the structure tensor 
  typedef  Matrix< double, 3, 3>                            EigenVectorMatrixType;
  typedef  Image< EigenVectorMatrixType, 3>                 EigenVectorImageType;
  typedef  itk::FixedArray< double, 3>                      EigenValueArrayType;
  typedef  itk::Image< EigenValueArrayType, 3>              EigenValueImageType;

  typedef  typename StructureTensorFilterType::OutputImageType  SymmetricSecondRankTensorImageType;
  typedef itk::
   SymmetricEigenVectorAnalysisImageFilter<SymmetricSecondRankTensorImageType, 
                                           EigenValueImageType, EigenVectorImageType> 
                    EigenVectorAnalysisFilterType;

  typename EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter = 
                                  EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( 3 );
  eigenVectorAnalysisFilter->OrderEigenValuesBy( 
      EigenVectorAnalysisFilterType::FunctorType::OrderByValue );
  
  eigenVectorAnalysisFilter->SetInput( m_StructureTensorFilter->GetOutput() );
  eigenVectorAnalysisFilter->Modified();
  eigenVectorAnalysisFilter->Update();

  //Step 1.3: Compute the eigen values 
  typedef itk::
    SymmetricEigenAnalysisImageFilter<SymmetricSecondRankTensorImageType, EigenValueImageType> 
                               EigenAnalysisFilterType;

  typename EigenAnalysisFilterType::Pointer eigenAnalysisFilter = EigenAnalysisFilterType::New();
  eigenAnalysisFilter->SetDimension( 3 );
  eigenAnalysisFilter->OrderEigenValuesBy( 
      EigenAnalysisFilterType::FunctorType::OrderByValue );
  
  eigenAnalysisFilter->SetInput( m_StructureTensorFilter->GetOutput() );
  eigenAnalysisFilter->Update();

  /* Step 2: Generate the diffusion tensor matrix 
      D = [v1 v2 v3] [DiagonalMatrixContainingLambdas] [v1 v2 v3]^t
  */

  //Setup the iterators
  //
  //Iterator for the eigenvector matrix image
  EigenVectorImageType::ConstPointer eigenVectorImage = 
                    eigenVectorAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenVectorImageType> eigenVectorImageIterator;
  eigenVectorImageIterator = itk::ImageRegionConstIterator<EigenVectorImageType>(
      eigenVectorImage, eigenVectorImage->GetRequestedRegion());
  eigenVectorImageIterator.GoToBegin();

  //Iterator for the diffusion tensor image
  typedef itk::ImageRegionIterator< DiffusionTensorImageType > DiffusionTensorIteratorType;
  DiffusionTensorIteratorType 
      it( m_DiffusionTensorImage, m_DiffusionTensorImage->GetLargestPossibleRegion() );

  //Iterator for the eigen value image
  EigenValueImageType::ConstPointer eigenImage = eigenAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenValueImageType> eigenValueImageIterator;
  eigenValueImageIterator = itk::ImageRegionConstIterator<EigenValueImageType>(
      eigenImage, eigenImage->GetRequestedRegion());
  eigenValueImageIterator.GoToBegin();

  it.GoToBegin();
  eigenVectorImageIterator.GoToBegin();
  eigenValueImageIterator.GoToBegin();

  MatrixType  eigenValueMatrix;
  while( !it.IsAtEnd() )
    {
    // Generate the diagonal matrix with the eigen values
    eigenValueMatrix.SetIdentity();

    //Set the lambda's appropriately. For now, set them to be equal to the eigen values
    double Lambda1;
    double Lambda2;
    double Lambda3;
  
    // Get the eigen value
    EigenValueArrayType eigenValue;
    eigenValue = eigenValueImageIterator.Get();
      
    Lambda1 = eigenValue[0];
    Lambda2 = eigenValue[1];
    Lambda3 = eigenValue[2];

    eigenValueMatrix(0,0) = Lambda1;
    eigenValueMatrix(1,1) = Lambda2;
    eigenValueMatrix(2,2) = Lambda3;

    //Get the eigenVector matrix
    EigenVectorMatrixType eigenVectorMatrix;
    eigenVectorMatrix = eigenVectorImageIterator.Get();

    EigenVectorMatrixType  eigenVectorMatrixTranspose;
    eigenVectorMatrixTranspose = eigenVectorMatrix.GetTranspose();

    // Generate the tensor matrix
    EigenVectorMatrixType  productMatrix;
    productMatrix = eigenVectorMatrix * eigenValueMatrix * eigenVectorMatrixTranspose;

    //Copy the ITK::Matrix to the tensor...there should be a better way of doing this TODO
    typename DiffusionTensorImageType::PixelType        tensor;

    tensor(0,0) = productMatrix(0,0);
    tensor(0,1) = productMatrix(0,1);
    tensor(0,2) = productMatrix(0,2);

    tensor(1,0) = productMatrix(1,0);
    tensor(1,1) = productMatrix(1,1);
    tensor(1,2) = productMatrix(1,2);

    tensor(2,0) = productMatrix(2,0);
    tensor(2,1) = productMatrix(2,1);
    tensor(2,2) = productMatrix(2,2);
    it.Set( tensor );

    ++it;
    ++eigenValueImageIterator;
    ++eigenVectorImageIterator;
    }
}

template<class TInputImage, class TOutputImage>
void
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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

  ThreadDiffusionImageRegionType    splitRegionDiffusionImage;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegionDiffusionImage);
  if (threadId < total)
    {
    str->Filter->ThreadedApplyUpdate(str->TimeStep,
                                     splitRegion, 
                                     splitRegionDiffusionImage,
                                     threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TInputImage, class TOutputImage>
typename
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>::TimeStepType
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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

  ThreadDiffusionImageRegionType splitDiffusionimageRegion;

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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType &regionToProcess,
                      const ThreadDiffusionImageRegionType &,
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>::TimeStepType
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
::ThreadedCalculateChange(const ThreadRegionType &regionToProcess, 
    const ThreadDiffusionImageRegionType & diffusionRegionToProcess, int)
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
     dynamic_cast<AnisotropicEdgeEnhancementDiffusionFunction<UpdateBufferType> *>
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

   // Process the non-boundary region.
  NeighborhoodIteratorType nD(radius, output, *fIt);

  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<DiffusionTensorImageType>
    DiffusionTensorFaceCalculatorType;

  typedef typename DiffusionTensorFaceCalculatorType::FaceListType 
                                                DiffusionTensorFaceListType;

  DiffusionTensorFaceCalculatorType diffusionTensorFaceCalculator;
  
  DiffusionTensorFaceListType diffusionTensorFaceList = 
     diffusionTensorFaceCalculator(m_DiffusionTensorImage, diffusionRegionToProcess, radius);

  typename DiffusionTensorFaceListType::iterator dfIt = diffusionTensorFaceList.begin();
  
  DiffusionTensorNeighborhoodType  dTN(radius,m_DiffusionTensorImage, *dfIt); 

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();


  UpdateIteratorType       nU(m_UpdateBuffer,  *fIt);
  nD.GoToBegin();
  while( !nD.IsAtEnd() )
    {
    nU.Value() = df->ComputeUpdate(nD, dTN, globalData);
    ++nD;
    ++nU;
    }

  // Process each of the boundary faces.

  NeighborhoodIteratorType bD;
  
  DiffusionTensorNeighborhoodType bDD;

  UpdateIteratorType   bU;
  for (++fIt; fIt != faceList.end(); ++fIt)
    {
    bD = NeighborhoodIteratorType(radius, output, *fIt);
    bDD = DiffusionTensorNeighborhoodType(radius, m_DiffusionTensorImage, *dfIt);
    bU = UpdateIteratorType  (m_UpdateBuffer, *fIt);
     
    bD.GoToBegin();
    bU.GoToBegin();
    while ( !bD.IsAtEnd() )
      {
      bU.Value() = df->ComputeUpdate(bD,bDD,globalData);
      ++bD;
      ++bU;
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
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
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
void
AnisotropicEdgeEnhancementDiffusionImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "TimeStep: " << m_TimeStep  << std::endl;
  os << indent << "Epsilon : " << m_Epsilon << std::endl;
  os << indent << "WStrength: " << m_WStrength  << std::endl;
  os << indent << "Sensitivity: " << m_Sensitivity  << std::endl;
}

}// end namespace itk

#endif
