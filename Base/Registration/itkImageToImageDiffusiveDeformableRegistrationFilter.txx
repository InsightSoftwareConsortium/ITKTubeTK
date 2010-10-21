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

#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"

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

  // TODO take me out
  m_DummyVector.Fill(0);

  m_UseDiffusiveRegularization      = true;
  m_BorderSurface                   = 0;
  m_BorderNormalsSurface            = 0;
  m_NormalVectorImage               = NormalVectorImageType::New();
  m_OutputTangentialImage           = OutputImageType::New();
  m_OutputNormalImage               = OutputImageType::New();
  m_TangentialDiffusionTensorImage  = DiffusionTensorImageType::New();
  m_NormalDiffusionTensorImage      = DiffusionTensorImageType::New();

  typename RegistrationFunctionType::Pointer registrationFunction =
                                                RegistrationFunctionType::New();
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                          registrationFunction.GetPointer() ) );

  // Setup the component extractor to extract the components from the deformation
  // field
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialComponentExtractor[i] = SelectionCastImageFilterType::New();
    m_TangentialComponentExtractor[i]->SetInput( m_OutputTangentialImage );
    m_TangentialComponentExtractor[i]->SetIndex( i );

    m_NormalComponentExtractor[i] = SelectionCastImageFilterType::New();
    m_NormalComponentExtractor[i]->SetInput( m_OutputNormalImage );
    m_NormalComponentExtractor[i]->SetIndex( i );
    }

  // Setup the deformation field component images
  for (unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_DeformationFieldTangentialComponents[i]
                              = DeformationFieldComponentImageType::New();
    m_DeformationFieldTangentialComponents[i]
                              = m_TangentialComponentExtractor[i]->GetOutput();
    m_DeformationFieldNormalComponents[i]
                              = DeformationFieldComponentImageType::New();
    m_DeformationFieldNormalComponents[i]
                              = m_NormalComponentExtractor[i]->GetOutput();
    }

  // By default, compute the intensity distance and regularization terms
  this->SetComputeIntensityDistanceTerm( true );
  this->SetComputeRegularizationTerm( true );
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
  Superclass::PrintSelf( os, indent );

  // TODO PrintSelf not compiling because of os type mismatch
  os << indent << "Border Surface: " << m_BorderSurface;
  //m_BorderSurface->PrintSelf( os, indent )
  os << indent << "Border Normals Surface: " << m_BorderNormalsSurface;
  //m_BorderNormalsSurface->PrintSelf( os, indent );
  os << indent << "UseDiffusiveRegularization: " << m_UseDiffusiveRegularization;

}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::SetComputeRegularizationTerm( bool compute )
{
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  df->SetComputeRegularizationTerm( compute );

  // Do not smooth the deformation field with the
  // PDEDeformableRegistrationFilter method if we are using our own
  // regularization term
  if ( compute )
    {
    this->SmoothDeformationFieldOff();
    this->SmoothUpdateFieldOff();
    }
  // TODO
  // Do smooth the deformation field with the default method if we are not
  // using our own regularization term
  else
    {
    this->SmoothDeformationFieldOff();
    //this->SmoothDeformationFieldOn();
    this->SmoothUpdateFieldOff();
    }
}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
bool
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetComputeRegularizationTerm() const
{
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  return df->GetComputeRegularizationTerm();
}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::SetComputeIntensityDistanceTerm( bool compute )
{
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  df->SetComputeIntensityDistanceTerm( compute );
}

/**
 * Set/Get whether to compute terms
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
bool
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetComputeIntensityDistanceTerm() const
{
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  return df->GetComputeIntensityDistanceTerm();
}

/**
 * Set/Get the border surface
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::BorderSurfacePointer
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetBorderSurface() const
{
  return m_BorderSurface;
}

/**
 * Set/Get the border surface
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::SetBorderSurface( BorderSurfacePointer border )
{
  m_BorderSurface = border;

  // Update the polydata of border surface normals
  vtkSmartPointer< vtkPolyDataNormals > normalExtractor
                                                  = vtkPolyDataNormals::New();
  normalExtractor->SetInput( m_BorderSurface );
  //normalExtractor->SetFeatureAngle(30);
  // NOTE: default settings compute point normals, not cell normals
  normalExtractor->Update();

  // extract generic(double) point normals
  m_BorderNormalsSurface = normalExtractor->GetOutput();
}

/**
 * Get the surface of border normals
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::BorderSurfacePointer
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetBorderNormalsSurface() const
{
  return m_BorderNormalsSurface;
}

/**
 * Get the image of normal vectors
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::NormalVectorImagePointer
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetNormalVectorImage() const
{
  return m_NormalVectorImage;
}

/**
 * Get the image of tangential diffusion tensors
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::DiffusionTensorImagePointer
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetTangentialDiffusionTensorImage() const
{
  return m_TangentialDiffusionTensorImage;
}

/**
 * Set/Get the image of tangential diffusion tensors
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::DiffusionTensorImagePointer
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetNormalDiffusionTensorImage() const
{
  return m_NormalDiffusionTensorImage;
}
// TODO create superclass with these methods for anisotropic diffusion
// registration filter

/**
 * Helper function to allocate space for an image given a template image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
template < class UnallocatedImageType, class TemplateImageType >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::AllocateSpaceForImage( UnallocatedImageType& inputImage,
                         const TemplateImageType& templateImage )
{
  inputImage->SetSpacing( templateImage->GetSpacing() );
  inputImage->SetOrigin( templateImage->GetOrigin() );
  inputImage->SetLargestPossibleRegion( templateImage->GetLargestPossibleRegion() );
  inputImage->SetRequestedRegion( templateImage->GetRequestedRegion() );
  inputImage->SetBufferedRegion( templateImage->GetBufferedRegion() );
  inputImage->Allocate();
}

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
  this->AllocateSpaceForImage( m_UpdateBuffer, output );
}

/**
 * All other initialization done before the initialize / calculate change /
 * apply update loop
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::Initialize()
{
  Superclass::Initialize();

  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the deformation field component images
  // TODO necessary?
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {

    // Allocate the tangential components
    this->AllocateSpaceForImage( m_DeformationFieldTangentialComponents[i],
                                 output);

    // Allocate the normal components
    this->AllocateSpaceForImage( m_DeformationFieldNormalComponents[i],
                                 output );
    }

  // Also allocate the diffusion tensor image
  // The diffusion tensor image has the same size as the deformation field and
  // holds the diffusion tensor matrix at each pixel
  DeformationFieldPointer deformationField = this->GetDeformationField();

  // Allocate the image of normals
  this->AllocateSpaceForImage( m_NormalVectorImage,
                               output );

  // Allocate the output image's tangential and normal images
  this->AllocateSpaceForImage( m_OutputTangentialImage,
                               output );
  this->AllocateSpaceForImage( m_OutputNormalImage,
                               output );

  // Allocate the tangential and normal diffusion tensor images
  this->AllocateSpaceForImage( m_TangentialDiffusionTensorImage,
                               output );
  this->AllocateSpaceForImage( m_NormalDiffusionTensorImage,
                               output );

  // We can't fail here if we don't have a border surface, but we will fail
  // at InitializeIteration();
  if ( m_BorderNormalsSurface )
    {
    // Compute the border normals and the weighting factor w
    // Normals are dependent on the border geometry in the fixed image so this
    // only has to be computed once.
    this->ComputeNormalVectorImage();

    // Compute the diffusion tensor image
    // The diffusion tensors are dependent on the normals computed in the
    // previous line, so this only has to be computed once.
    this->ComputeDiffusionTensorImage();
    }
}

// TODO halting criteria?!?!
// TODO: SetUseImageSpacing() to on for this filter?  Would give derivates in physical space, default is off

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
  std::cout << "\tInitializeIteration for FILTER" << std::endl;

  if ( !this->GetFixedImage() || !this->GetMovingImage()
        || !this->GetDeformationField() || !this->GetBorderSurface()
        || !this->GetBorderNormalsSurface() )
    {
    itkExceptionMacro( << "FixedImage, MovingImage, DeformationField, border "
                       << "surface and/or border normals surface not set");
    }

  // TODO checking the timestep for stability as in the anisotropic filter

  // Update the deformation field component images
  // This depends on the current deformation field u, so it must be computed
  // on every iteration of the filter.
  this->UpdateDeformationFieldComponentImages();

  // Update the function's deformation field
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  df->SetDeformationField( this->GetDeformationField() );

  // Call the superclass implementation
  Superclass::InitializeIteration();

}

/**
 * Updates the border normals and the weighting factor w
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                  TDeformationField >
::ComputeNormalVectorImage()
{
  // TODO assert vs. if?
  assert( m_BorderNormalsSurface );

  // Get the normals
  vtkSmartPointer< vtkDataArray > normalData
                      = m_BorderNormalsSurface->GetPointData()->GetNormals();

  // Iterate over the normal vector image and insert the normal of the closest
  // point
  NormalVectorIteratorType normalVectorIt( m_NormalVectorImage,
                              m_NormalVectorImage->GetLargestPossibleRegion() );

  vtkPointLocator * locator = vtkPointLocator::New();
  locator->SetDataSet( m_BorderNormalsSurface );

  itk::Index< ImageDimension > index;
  typename NormalVectorImageType::SpacingType spacing
                                            = m_NormalVectorImage->GetSpacing();
  typename NormalVectorImageType::PointType origin
                                            = m_NormalVectorImage->GetOrigin();
  double coord[3];
  vtkIdType id;

  std::cout << "Computing normals... " << std::endl;

  for( normalVectorIt.GoToBegin(); !normalVectorIt.IsAtEnd(); ++normalVectorIt )
    {
    index = normalVectorIt.GetIndex();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      coord[i] = ( index[i] * spacing[i] ) + origin[i];
      }
    id = locator->FindClosestPoint(coord);
    normalVectorIt.Set( normalData->GetTuple( id ) );
    }

  std::cout << "Finished computing normals." << std::endl;

  // TODO calculate w here
  //DeformationFieldScalarType w = (DeformationFieldScalarType) 1.0;
  // TODO don't compute if m_UseDiffusiveRegularization is false
}

/**
 * Updates the diffusion tensor image before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ComputeDiffusionTensorImage()
{
  typedef itk::Matrix< DeformationFieldScalarType,
                       ImageDimension, ImageDimension > MatrixType;

  NormalVectorType                n;
  DeformationFieldScalarType      w;

  // Used to compute the tangential and normal diffusion tensor images

  // tangential:
  // P = I - wnn^T
  // tangentialMatrix = tangentialD = P^TP

  // normal:
  // normalMatrix = normalD = wnn^T

  MatrixType                      normalMatrix;
  MatrixType                      P;
  MatrixType                      tangentialMatrix;
  DiffusionTensorImagePixelType   tangentialD;
  DiffusionTensorImagePixelType   normalD;

  NormalVectorIteratorType normalVectorIt( m_NormalVectorImage,
                              m_NormalVectorImage->GetLargestPossibleRegion() );

  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
                                                    DiffusionTensorIteratorType;
  DiffusionTensorIteratorType tangentialIt( m_TangentialDiffusionTensorImage,
                 m_TangentialDiffusionTensorImage->GetLargestPossibleRegion() );
  DiffusionTensorIteratorType normalIt( m_NormalDiffusionTensorImage,
                 m_NormalDiffusionTensorImage->GetLargestPossibleRegion() );

  for( normalVectorIt.GoToBegin(), tangentialIt.GoToBegin(), normalIt.GoToBegin();
        !tangentialIt.IsAtEnd(); ++normalVectorIt, ++tangentialIt, ++normalIt )
    {

    // 1.  Get the border normal n and the weighting factor w
    n = normalVectorIt.Get();
    if ( !m_UseDiffusiveRegularization )
      {
      w = ( DeformationFieldScalarType ) 0.0;
      }
    else
      {
      // TODO get w here
      w = (DeformationFieldScalarType) 1.0;
      }

    // 2. Compute the tangential and normal diffusion tensor images

    // Create the nMatrix used to calculate nn^T
    // (The first column is filled with the values of n, the rest are 0s)
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      normalMatrix(i,0) = n[i];
      for ( unsigned int j = 1; j < ImageDimension; j++ )
        {
        normalMatrix(i,j) = 0;
        }
      }

    normalMatrix = normalMatrix * normalMatrix.GetTranspose(); // nn^T
    normalMatrix = normalMatrix * w; // wnn^T
    P.SetIdentity();
    P = P - normalMatrix; // I - wnn^T
    tangentialMatrix = P.GetTranspose();
    tangentialMatrix = tangentialMatrix * P; // P^TP

    // Copy the itk::Matrix to the tensor - there should be a better way to do
    // this
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        tangentialD( i,j ) = tangentialMatrix( i,j );
        normalD( i,j ) = normalMatrix( i,j );
        }
      }
    // Copy the diffusion tensors to m_TangentialDiffusionTensorImage and
    // m_NormalDiffusionTensorImage
    tangentialIt.Set( tangentialD );
    normalIt.Set( normalD );
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
::UpdateDeformationFieldComponentImages()
{
  // Get the border normals
  NormalVectorIteratorType normalVectorIterator( m_NormalVectorImage,
                              m_NormalVectorImage->GetLargestPossibleRegion() );

  // Get output (the current deformation field)
  typename OutputImageType::Pointer output = this->GetOutput();

  typedef itk::ImageRegionIterator< OutputImageType > IteratorType;
  IteratorType outputImageIterator(
                          output,
                          output->GetLargestPossibleRegion() );

  // Extract tangential and normal components from output
  IteratorType outputTangentialImageIterator(
                          m_OutputTangentialImage,
                          m_OutputTangentialImage->GetLargestPossibleRegion() );
  IteratorType outputNormalImageIterator(
                          m_OutputNormalImage,
                          m_OutputNormalImage->GetLargestPossibleRegion() );

  // Calculate the tangential and normal components of the deformation field
  NormalVectorType            n;
  DeformationFieldVectorType  u;
  DeformationFieldVectorType  normalU;
  DeformationFieldVectorType  tangentialU;

  normalVectorIterator.GoToBegin();
  outputTangentialImageIterator.GoToBegin();
  outputNormalImageIterator.GoToBegin();
  for( outputImageIterator.GoToBegin(); !outputImageIterator.IsAtEnd();
         ++outputImageIterator )
    {
    n = normalVectorIterator.Get();
    u = outputImageIterator.Get();

    // normal component = (u^Tn)n
    normalU = (u * n) * n;
    outputNormalImageIterator.Set( normalU );

    // tangential component = u - normal component
    tangentialU = u - normalU;
    outputTangentialImageIterator.Set( tangentialU );

    // We know normalU + tangentialU = u
    // Assertion to test that the normal and tangential components were computed
    // corectly - they should be orthogonal
    assert( normalU * tangentialU < 0.000005 );

    ++normalVectorIterator;
    ++outputTangentialImageIterator;
    ++outputNormalImageIterator;
    }

  // Feed tangential and normal components to the component extractor
  // Extract the components
  m_OutputTangentialImage->Modified();
  m_OutputNormalImage->Modified();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialComponentExtractor[i]->Update();
    m_NormalComponentExtractor[i]->Update();
    }

  // Little test to make sure that the component extractor was setup
  // properly - because m_DeformationFieldTangentialComponents and
  // m_DeformationFieldNormalComponents will be the input to the
  // registration function's ComputeUpdate()
  typename DeformationFieldComponentImageType::IndexType index;
  index.Fill(0);
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    assert( m_DeformationFieldTangentialComponents[i]->GetPixel( index )
        == m_TangentialComponentExtractor[i]->GetOutput()->GetPixel( index ) );
    assert( m_DeformationFieldNormalComponents[i]->GetPixel( index )
        == m_NormalComponentExtractor[i]->GetOutput()->GetPixel( index ) );
    }
}

/**
 * Populates the update buffer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >

::TimeStepType
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
  // TODO put back
  //this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetNumberOfThreads(1);
  this->GetMultiThreader()->SetSingleMethod(this->CalculateChangeThreaderCallback,
                                            &str);

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  // TODO put back
  //threadCount = this->GetMultiThreader()->GetNumberOfThreads();
  threadCount = 1;

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

  // TODO why did Andinet comment this out?
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
  // TODO put back
  //threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
  threadCount = 1;

  str = (DenseFDThreadStruct *)
                        (((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;

  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  ThreadNormalVectorImageRegionType splitNormalVectorImageRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitNormalVectorImageRegion);

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
                                      splitNormalVectorImageRegion,
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
typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ThreadedCalculateChange(
          const ThreadRegionType &regionToProcess,
          const ThreadNormalVectorImageRegionType &normalVectorRegionToProcess,
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

  // Setup the boundary faces for the normal vector images
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
                              < NormalVectorImageType >
                              NormalVectorFaceCalculatorType;

  typedef typename NormalVectorFaceCalculatorType::FaceListType
                              NormalVectorFaceListType;

  NormalVectorFaceCalculatorType
                              normalVectorFaceCalculator;

  NormalVectorFaceListType    normalVectorFaceList
                    = normalVectorFaceCalculator( m_NormalVectorImage,
                                                  normalVectorRegionToProcess,
                                                  radius );

  typename NormalVectorFaceListType::iterator
                  normalVectorFaceListIt = normalVectorFaceList.begin();

  // Setup the boundary faces for the deformation field component images
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
                              < DeformationFieldComponentImageType >
                              DeformationFieldComponentFaceCalculatorType;

  typedef typename DeformationFieldComponentFaceCalculatorType::FaceListType
                              DeformationFieldComponentFaceListType;

  DeformationFieldComponentFaceCalculatorType
                              deformationFieldTangentialComponentFaceCalculator;

  itk::FixedArray< DeformationFieldComponentFaceListType, ImageDimension >
                              deformationFieldTangentialComponentFaceList;

  typedef typename DeformationFieldComponentFaceListType::iterator
                              DeformationFieldFaceListIterator;
  itk::FixedArray< DeformationFieldFaceListIterator, ImageDimension >
                              deformationFieldTangentialComponentFaceListIterator;

  DeformationFieldComponentFaceCalculatorType
                              deformationFieldNormalComponentFaceCalculator;

  itk::FixedArray< DeformationFieldComponentFaceListType, ImageDimension >
                              deformationFieldNormalComponentFaceList;

  itk::FixedArray< DeformationFieldFaceListIterator, ImageDimension >
                              deformationFieldNormalComponentFaceListIterator;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    deformationFieldTangentialComponentFaceList[i]
        = deformationFieldTangentialComponentFaceCalculator(
                        m_DeformationFieldTangentialComponents[i],
                        deformationComponentRegionToProcess,
                        radius );
    deformationFieldTangentialComponentFaceListIterator[i] =
                        deformationFieldTangentialComponentFaceList[i].begin();

    deformationFieldNormalComponentFaceList[i]
        = deformationFieldNormalComponentFaceCalculator(
                        m_DeformationFieldNormalComponents[i],
                        deformationComponentRegionToProcess,
                        radius );
    deformationFieldNormalComponentFaceListIterator[i] =
                        deformationFieldNormalComponentFaceList[i].begin();

    }

  // Setup the boundary faces calculator for the diffusion tensor image
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
                                        < DiffusionTensorImageType >
                                        DiffusionTensorFaceCalculatorType;

  typedef typename DiffusionTensorFaceCalculatorType::FaceListType
                                        DiffusionTensorFaceListType;

  DiffusionTensorFaceCalculatorType diffusionTensorTangentialFaceCalculator;

  DiffusionTensorFaceListType diffusionTensorTangentialFaceList
      = diffusionTensorTangentialFaceCalculator(
                                        m_TangentialDiffusionTensorImage,
                                        diffusionRegionToProcess, radius );

  typename DiffusionTensorFaceListType::iterator tangentialDfIt
                                    = diffusionTensorTangentialFaceList.begin();

  DiffusionTensorFaceCalculatorType diffusionTensorNormalFaceCalculator;

  DiffusionTensorFaceListType diffusionTensorNormalFaceList
      = diffusionTensorNormalFaceCalculator(
                                        m_NormalDiffusionTensorImage,
                                        diffusionRegionToProcess, radius );

  typename DiffusionTensorFaceListType::iterator normalDfIt
                                    = diffusionTensorNormalFaceList.begin();

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();

  // Process the non-boundary region.
  NeighborhoodIteratorType        nD(radius, output, *fIt);
  UpdateIteratorType              nU(m_UpdateBuffer,  *fIt);

  NormalVectorImageNeighborhoodType normalVectorN(radius,
                                                  m_NormalVectorImage,
                                                  *normalVectorFaceListIt);

  DiffusionTensorNeighborhoodType tangentialDTN(radius,
                                                m_TangentialDiffusionTensorImage,
                                                *tangentialDfIt);
  DeformationFieldComponentNeighborhoodArrayType tangentialDFC;

  DiffusionTensorNeighborhoodType normalDTN(radius,
                                            m_NormalDiffusionTensorImage,
                                            *normalDfIt);
  DeformationFieldComponentNeighborhoodArrayType normalDFC;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    tangentialDFC[i] = DeformationFieldComponentNeighborhoodType(
                      radius,
                      m_DeformationFieldTangentialComponents[i],
                      *deformationFieldTangentialComponentFaceListIterator[i]);

    normalDFC[i] = DeformationFieldComponentNeighborhoodType(
                      radius,
                      m_DeformationFieldNormalComponents[i],
                      *deformationFieldNormalComponentFaceListIterator[i]);
    }

  nD.GoToBegin();
  nU.GoToBegin();
  normalVectorN.GoToBegin();
  tangentialDTN.GoToBegin();
  normalDTN.GoToBegin();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    tangentialDFC[i].GoToBegin();
    normalDFC[i].GoToBegin();
    }
  while( !nD.IsAtEnd() )
    {
    nU.Value() = df->ComputeUpdate(nD,
                                   normalVectorN,
                                   tangentialDTN,
                                   tangentialDFC,
                                   normalDTN,
                                   normalDFC,
                                   globalData);

    ++nD;
    ++nU;
    ++normalVectorN;
    ++tangentialDTN;
    ++normalDTN;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      ++tangentialDFC[i];
      ++normalDFC[i];
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
//                                          m_TangentialDiffusionTensorImage,
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
  // TODO put back
  //threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
  threadCount = 1;

  str = (DenseFDThreadStruct *)
                        (((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

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
