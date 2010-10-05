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
#ifndef _itkImageToImageDiffusiveDeformableRegistrationFilter_txx
#define _itkImageToImageDiffusiveDeformableRegistrationFilter_txx

#include "itkImageToImageDiffusiveDeformableRegistrationFilter.h"

#include "itkVectorIndexSelectionCastImageFilter.h"


namespace itk
{ 
    
/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ImageToImageDiffusiveDeformableRegistrationFilter()
{
  m_UpdateBuffer = UpdateBufferType::New();

  m_Normals.Fill( 0 );
  m_TangentalDiffusionTensorImage = DiffusionTensorImageType::New();

  typename RegistrationFunctionType::Pointer registrationFunction =
                                                RegistrationFunctionType::New();
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                          registrationFunction.GetPointer() ) );

  // Setup the component extractor to extract the components from the deformation
  // field
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_ComponentExtractor[i] = SelectionCastImageFilterType::New();
    m_ComponentExtractor[i]->SetInput( this->GetOutput() );
    m_ComponentExtractor[i]->SetIndex( i );
    }

  // Setup the deformation field component images
  for (unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_DeformationFieldComponents[i] = DeformationFieldComponentImageType::New();
    m_DeformationFieldComponents[i] = m_ComponentExtractor[i]->GetOutput();
    }
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Normals: ";
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    std::cout << m_Normals[i] << " ";
    }
}

/**
 * Set/Get the border normals
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::NormalVectorType&
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetNormals() const
{
  return m_Normals;
}

/**
 * Set/Get the border normals
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::SetNormals( NormalVectorType& normals )
{
  m_Normals = normals;
}

/**
 * Set/Get the image of tangental diffusion tensors
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::DiffusionTensorImagePointer
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetTangentalDiffusionTensorImage() const
{
  return m_TangentalDiffusionTensorImage;
}

// TODO create superclass with these methods for anisotropic diffusion registration filter

/**
 * Allocate space for the update buffer, and the diffusion tensor image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::AllocateUpdateBuffer()
{

  /* The update buffer looks just like the output and holds the change in
   the pixel  */

  typename OutputImageType::Pointer output = this->GetOutput();

  m_UpdateBuffer->SetSpacing(output->GetSpacing());
  m_UpdateBuffer->SetOrigin(output->GetOrigin());
  m_UpdateBuffer->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_UpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_UpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_UpdateBuffer->Allocate();

  // Allocate the deformation field component images
  // TODO necessary?
  this->AllocateDeformationFieldComponentImages();

  // Also allocate the diffusion tensor image
  this->AllocateDiffusionTensorImage();
}

/**
 * Allocate space for the deformation field components
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::AllocateDeformationFieldComponentImages()
{
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    DeformationFieldPointer deformationField = this->GetDeformationField();
    m_DeformationFieldComponents[i]->SetOrigin( deformationField->GetOrigin() );
    m_DeformationFieldComponents[i]->SetSpacing( deformationField->GetSpacing() );
    m_DeformationFieldComponents[i]->SetLargestPossibleRegion(
                                  deformationField->GetLargestPossibleRegion() );
    m_DeformationFieldComponents[i]->SetRequestedRegion(
                                  deformationField->GetRequestedRegion() );
    m_DeformationFieldComponents[i]->SetBufferedRegion(
                                  deformationField->GetBufferedRegion() );
    m_DeformationFieldComponents[i]->Allocate();
    }
}

/**
 * Allocate space for the diffusion tensor image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::AllocateDiffusionTensorImage()
{
  std::cout << "AllocatDiffusionTensorImage for FILTER"<< std::endl;

  // The diffusion tensor image has the same size as the deformation field and
  // holds the diffusion tensor matrix at each pixel
  DeformationFieldPointer deformationField = this->GetDeformationField();
  m_TangentalDiffusionTensorImage->SetOrigin( deformationField->GetOrigin() );
  m_TangentalDiffusionTensorImage->SetSpacing( deformationField->GetSpacing() );
  m_TangentalDiffusionTensorImage->SetLargestPossibleRegion(
                                deformationField->GetLargestPossibleRegion() );
  m_TangentalDiffusionTensorImage->SetRequestedRegion(
                                deformationField->GetRequestedRegion() );
  m_TangentalDiffusionTensorImage->SetBufferedRegion(
                                deformationField->GetBufferedRegion() );
  m_TangentalDiffusionTensorImage->Allocate();
}

/**
 * Initialize the state of the filter and equation before each iteration.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::InitializeIteration()
{
  std::cout << "InitializeIteration for FILTER" << std::endl;

  // Call the superclass implementation
  Superclass::InitializeIteration();

  if ( !this->GetFixedImage() || !this->GetMovingImage()
        || !this->GetDeformationField() )
    {
    itkExceptionMacro( << "FixedImage, MovingImage and/or DeformationField not set");
    }

  // TODO checking the timestep for stability as in the anisotropic filter

  // Update the diffusion tensor image
  this->UpdateDiffusionTensorImage();

  // Extract the components
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_ComponentExtractor[i]->Update();
    }
}

/**
 * Updates the diffusion tensor image before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::UpdateDiffusionTensorImage()
{
  typedef itk::Matrix< DeformationFieldScalarType,
                       ImageDimension, ImageDimension > MatrixType;

  // 1. Compute tangental diffusion tensor image

  // P = I - wnn^T
  NormalVectorType                n;
  DeformationFieldScalarType      w;
  MatrixType                      nMatrix;
  MatrixType                      P;

  // result = D = P^TP
  MatrixType                      result;
  DiffusionTensorImagePixelType   D;

  typedef itk::ImageRegionIterator< DiffusionTensorImageType > IteratorType;
  IteratorType it( m_TangentalDiffusionTensorImage,
                   m_TangentalDiffusionTensorImage->GetLargestPossibleRegion() );
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    // TODO calculate w and n here
    n = this->GetNormals();
    w = (DeformationFieldScalarType) 1.0;

    // Create the nMatrix used to calculate nn^T
    // The first column is filled with the values of n, the rest are 0s
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      nMatrix(i,0) = n[i];
      for ( unsigned int j = 1; j < ImageDimension; j++ )
        {
        nMatrix(i,j) = 0;
        }
      }

    nMatrix = nMatrix * nMatrix.GetTranspose(); // nn^T
    nMatrix = nMatrix * w; // wnn^T
    P.SetIdentity();
    P = P - nMatrix; // I - wnn^T
    result = P.GetTranspose();
    result = result * P; // P^TP

    // Copy the itk::Matrix to the tensor - there should be a better way to do
    // this
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        D(i,j) = result(i,j);
        }
      }


    // TODO test with more complicated n, then take me out
//    std::cout << "*****" << std::endl;
//    for (unsigned int i = 0; i < ImageDimension; i++)
//      {
//      for (unsigned int j = 0; j < ImageDimension; j++)
//        {
//        std::cout << D(i,j) << " ";
//        }
//      std::cout << std::endl;
//      }
//    std::cout << "*****" << std::endl;


    // Copy the diffusion tensor to m_TangentalDiffusionTensorImage
    it.Set(D);
    }
}

/**
 * Populates the update buffer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >::TimeStepType
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::CalculateChange()
{
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

  // Explicitely call Modified on m_UpdateBuffer here
  // since ThreadedCalculateChange changes this buffer
  // through iterators which do not increment the
  // update buffer timestamp
  //this->m_UpdateBuffer->Modified();

  return  dt;
}

/**
 * Calls ThreadedCalculateChange for processing
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_TYPE
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
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

  ThreadDiffusionTensorImageRegionType splitDiffusionImageRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitDiffusionImageRegion);

  ThreadDeformationFieldComponentImageRegionType
                                      splitDeformationFieldComponentImageRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                      splitDeformationFieldComponentImageRegion);

  if (threadId < total)
    {
    str->TimeStepList[threadId]
      = str->Filter->ThreadedCalculateChange(
                                      splitRegion,
                                      splitDiffusionImageRegion,
                                      splitDeformationFieldComponentImageRegion,
                                      threadId);
    str->ValidTimeStepList[threadId] = true;
    }

  return ITK_THREAD_RETURN_VALUE;
}

/**
 * Does the actual work of calculating change over a region supplied by the
 * multithreading mechanism
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >::TimeStepType
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ThreadedCalculateChange(
          const ThreadRegionType &regionToProcess,
          const ThreadDiffusionTensorImageRegionType &diffusionRegionToProcess,
          const ThreadDeformationFieldComponentImageRegionType
                                          &deformationComponentRegionToProcess,
          int)
{

  typedef typename OutputImageType::RegionType        RegionType;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::SizeValueType     SizeValueType;
  typedef typename OutputImageType::IndexType         IndexType;
  typedef typename OutputImageType::IndexValueType    IndexValueType;
  typedef typename
    FiniteDifferenceFunctionType::NeighborhoodType    NeighborhoodIteratorType;
  typedef ImageRegionIterator<UpdateBufferType>       UpdateIteratorType;

  typename OutputImageType::Pointer output = this->GetOutput();

  TimeStepType timeStep;
  void *globalData;

  // Get the FiniteDifferenceFunction to use in calculations.
  const typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );

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

  // Setup the boundary faces for the deformation field component images
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
                                    < DeformationFieldComponentImageType >
                                    DeformationFieldComponentFaceCalculatorType;

  typedef typename DeformationFieldComponentFaceCalculatorType::FaceListType
                                    DeformationFieldComponentFaceListType;

  DeformationFieldComponentFaceCalculatorType
                                    deformationFieldComponentFaceCalculator;

  itk::FixedArray< DeformationFieldComponentFaceListType, ImageDimension >
                                    deformationFieldComponentFaceList;

  typedef typename DeformationFieldComponentFaceListType::iterator
                                    DeformationFieldFaceListIterator;
  itk::FixedArray< DeformationFieldFaceListIterator, ImageDimension >
                                    deformationFieldComponentFaceListIterator;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    deformationFieldComponentFaceList[i]
        = deformationFieldComponentFaceCalculator(
                                          m_DeformationFieldComponents[i],
                                          deformationComponentRegionToProcess,
                                          radius );
    deformationFieldComponentFaceListIterator[i] =
                                  deformationFieldComponentFaceList[i].begin();

    }

  // Setup the boundary faces calculator for the diffusion tensor image
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
                                        < DiffusionTensorImageType >
                                        DiffusionTensorFaceCalculatorType;

  typedef typename DiffusionTensorFaceCalculatorType::FaceListType
                                        DiffusionTensorFaceListType;

  DiffusionTensorFaceCalculatorType diffusionTensorFaceCalculator;

  DiffusionTensorFaceListType diffusionTensorFaceList
      = diffusionTensorFaceCalculator( m_TangentalDiffusionTensorImage,
                                       diffusionRegionToProcess, radius );

  typename DiffusionTensorFaceListType::iterator dfIt
                                              = diffusionTensorFaceList.begin();

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();

  // Process the non-boundary region.
  NeighborhoodIteratorType        nD(radius, output, *fIt);
  UpdateIteratorType              nU(m_UpdateBuffer,  *fIt);
  DiffusionTensorNeighborhoodType dTN(radius, m_TangentalDiffusionTensorImage,
                                      *dfIt);
  DeformationFieldComponentNeighborhoodArrayType dFC;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    dFC[i] = DeformationFieldComponentNeighborhoodType(
                                radius,
                                m_DeformationFieldComponents[i],
                                *deformationFieldComponentFaceListIterator[i]);
    dFC[i].GoToBegin();
    }

  nD.GoToBegin();
  nU.GoToBegin();
  dTN.GoToBegin();
  while( !nD.IsAtEnd() )
    {
    nU.Value() = df->ComputeUpdate(nD, dTN, dFC, globalData);
    ++nD;
    ++nU;
    ++dTN;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      ++dFC[i];
      }
    }

  // Process each of the boundary faces.

  // TODO

//  NeighborhoodIteratorType bD;
//  DiffusionTensorNeighborhoodType bDD;
//  UpdateIteratorType   bU;
//  ++dfIt;

//  for (++fIt; fIt != faceList.end(); ++fIt)
//    {
//    bD = NeighborhoodIteratorType(radius, output, *fIt);
//    bDD = DiffusionTensorNeighborhoodType(radius,
//                                          m_TangentalDiffusionTensorImage,
//                                          *dfIt);

//    bU = UpdateIteratorType  (m_UpdateBuffer, *fIt);

//    bD.GoToBegin();
//    bU.GoToBegin();
//    bDD.GoToBegin();
//    while ( !bD.IsAtEnd() )
//      {
//      bU.Value() = df->ComputeUpdate(bD, bDD, globalData);
//      ++bD;
//      ++bU;
//      ++bDD;
//      }
//    ++dfIt;
//    }

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  timeStep = df->ComputeGlobalTimeStep(globalData);
  df->ReleaseGlobalDataPointer(globalData);

  return timeStep;
}

/**
 * Applies changes from the update buffer to the output
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ApplyUpdate(TimeStepType dt)
{
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = dt;
  // TODO
  //this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetNumberOfThreads(1);
  this->GetMultiThreader()->SetSingleMethod(this->ApplyUpdateThreaderCallback,
                                            &str);
  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Explicitely call Modified on GetOutput here
  // since ThreadedApplyUpdate changes this buffer
  // through iterators which do not increment the
  // output timestamp
  this->GetOutput()->Modified();
}

/**
 * Calls ThreadedApplyUpdate, need to reimplement here to also split the
 * diffusion tensor image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_TYPE
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
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

  ThreadDiffusionTensorImageRegionType splitRegionDiffusionTensorImage;
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

/**
 * Does the actual work of updating the output from the UpdateContainer over an
 * output region supplied by the multithreading mechanism.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType &regionToProcess,
                      const ThreadDiffusionTensorImageRegionType &,
                      int threadId)
{
  ImageRegionIterator<UpdateBufferType> u(m_UpdateBuffer,    regionToProcess);
  ImageRegionIterator<OutputImageType>  o(this->GetOutput(), regionToProcess);

  u = u.Begin();
  o = o.Begin();

  while ( !u.IsAtEnd() )
    {

    o.Value() += static_cast<DeformationFieldVectorType>(u.Value() * dt);
                                                    // no adaptor support here

    ++o;
    ++u;
    }
}






} // end namespace itk

#endif
