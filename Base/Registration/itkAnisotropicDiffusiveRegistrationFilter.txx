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
#ifndef _itkAnisotropicDiffusiveRegistrationFilter_txx
#define _itkAnisotropicDiffusiveRegistrationFilter_txx

#include "itkAnisotropicDiffusiveRegistrationFilter.h"

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"

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
  m_UpdateBuffer = UpdateBufferType::New();

  m_BorderSurface                   = 0;
  m_BorderNormalsSurface            = 0;
  m_BorderNormalsSurfaceFilter      = 0;
  m_NormalVectorImage               = 0;
  m_WeightImage                     = 0;
  m_NormalDeformationField          = 0;
  m_TangentialDiffusionTensorImage  = 0;
  m_NormalDiffusionTensorImage      = 0;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialComponentExtractor[i]           = 0;
    m_DeformationVectorTangentialComponents[i]  = 0;
    m_NormalComponentExtractor[i]               = 0;
    m_DeformationVectorNormalComponents[i]      = 0;
    }

  typename RegistrationFunctionType::Pointer registrationFunction =
      RegistrationFunctionType::New();
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
      registrationFunction.GetPointer() ) );

  // Lambda for exponential decay used to calculate weight from distance.
  m_lambda = -0.01;

  // By default, compute the intensity distance and regularization terms
  this->SetComputeIntensityDistanceTerm( true );
  this->SetComputeRegularizationTerm( true );
  this->SetUseAnisotropicRegularization( true );

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
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Border Surface: " << m_BorderSurface;
  //m_BorderSurface->PrintSelf( os, indent )
  os << indent << "Border Normals Surface: " << m_BorderNormalsSurface;
  //m_BorderNormalsSurface->PrintSelf( os, indent );
}

/**
 * Get the registration function pointer
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::RegistrationFunctionPointer
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetRegistrationFunctionPointer() const
{
  RegistrationFunctionPointer df = dynamic_cast< RegistrationFunctionType * >
       ( this->GetDifferenceFunction().GetPointer() );

  if ( !df )
    {
    itkExceptionMacro( << "Registration function pointer NULL" << std::endl );
    }

  return df;
}

/**
 * Helper function to allocate space for an image given a template image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
template < class UnallocatedImageType, class TemplateImageType >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateSpaceForImage( UnallocatedImageType& inputImage,
                         const TemplateImageType& templateImage )
{
  inputImage->SetSpacing( templateImage->GetSpacing() );
  inputImage->SetOrigin( templateImage->GetOrigin() );
  inputImage->SetLargestPossibleRegion(
      templateImage->GetLargestPossibleRegion() );
  inputImage->SetRequestedRegion( templateImage->GetRequestedRegion() );
  inputImage->SetBufferedRegion( templateImage->GetBufferedRegion() );
  inputImage->Allocate();
}

/**
 * Helper function to check whether the attributes of an image matches template
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
template < class CheckedImageType, class TemplateImageType >
bool
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CompareImageAttributes( const CheckedImageType& inputImage,
                          const TemplateImageType& templateImage )
{
  return inputImage->GetSpacing() == templateImage->GetSpacing()
      && inputImage->GetOrigin() == templateImage->GetOrigin()
      && inputImage->GetLargestPossibleRegion()
          == templateImage->GetLargestPossibleRegion()
      && inputImage->GetRequestedRegion() == templateImage->GetRequestedRegion()
      && inputImage->GetBufferedRegion() == templateImage->GetBufferedRegion();
}

/**
 * Allocate space for the update buffer, and the diffusion tensor image
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
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
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::Initialize()
{
  Superclass::Initialize();

  // Check the timestep for stability
  if( this->GetComputeRegularizationTerm() )
    {
    this->GetRegistrationFunctionPointer()->CheckTimeStepStability(
        this->GetInput(), this->GetUseImageSpacing() );
    }

  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the output image's normal images, tangential and normal diffusion
  // tensor images, and deformation field component images
  if( this->GetComputeRegularizationTerm() )
    {
    m_TangentialDiffusionTensorImage = DiffusionTensorImageType::New();
    this->AllocateSpaceForImage( m_TangentialDiffusionTensorImage, output );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      m_TangentialComponentExtractor[i] = SelectionCastImageFilterType::New();
      m_TangentialComponentExtractor[i]->SetInput( this->GetOutput() );
      m_TangentialComponentExtractor[i]->SetIndex( i );
      m_DeformationVectorTangentialComponents[i]
          = m_TangentialComponentExtractor[i]->GetOutput();
      }
    if( this->GetUseAnisotropicRegularization() )
      {
      m_NormalDeformationField = OutputImageType::New();
      this->AllocateSpaceForImage( m_NormalDeformationField, output );
      m_NormalDiffusionTensorImage = DiffusionTensorImageType::New();
      this->AllocateSpaceForImage( m_NormalDiffusionTensorImage, output );
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        m_NormalComponentExtractor[i] = SelectionCastImageFilterType::New();
        m_NormalComponentExtractor[i]->SetInput( m_NormalDeformationField );
        m_NormalComponentExtractor[i]->SetIndex( i );
        m_DeformationVectorNormalComponents[i]
            = m_NormalComponentExtractor[i]->GetOutput();
        }
      }
    }

  if( this->GetComputeRegularizationTerm() && this->GetUseAnisotropicRegularization() )
    {
    // Check the normal vector image and the weight image if one was supplied
    // by the user
    if( m_NormalVectorImage
          && !this->CompareImageAttributes( m_NormalVectorImage, output ) )
      {
      itkExceptionMacro( << "Normal vector image does not have the same "
                         << "attributes as the output deformation field"
                         << std::endl );
      }
    if( m_WeightImage
          && !this->CompareImageAttributes( m_WeightImage, output ) )
      {
      itkExceptionMacro( << "Weight image does not have the same attributes as "
                         << "the output deformation field" << std::endl );
      }

    // Compute the border normals and/or weighting factors if not supplied by the
    // user
    if( !m_NormalVectorImage || !m_WeightImage )
      {
      if ( !this->GetBorderSurface() )
        {
        itkExceptionMacro( << "Cannot perform registration without a border "
                           << "surface, or a normal vector image and a weight "
                           << "image" << std::endl );
        }

      // Setup the vtkPolyDataNormals to extract the normals from the surface
      m_BorderNormalsSurfaceFilter = BorderNormalsSurfaceFilterType::New();
      m_BorderNormalsSurfaceFilter->ComputePointNormalsOn();
      m_BorderNormalsSurfaceFilter->ComputeCellNormalsOff();
      //m_BorderNormalsSurfaceFilter->SetFeatureAngle(30);
      m_BorderNormalsSurfaceFilter->SetInput( m_BorderSurface );
      m_BorderNormalsSurfaceFilter->Update();
      m_BorderNormalsSurface = m_BorderNormalsSurfaceFilter->GetOutput();

      // Make sure we now have the normals
      if ( !this->GetBorderNormalsSurface() )
        {
        itkExceptionMacro( << "Error computing border normals" << std::endl );
        }

      bool computeNormalVectorImage = false;
      bool computeWeightImage = false;

      // Allocate the image of normals and weight image
      if( !m_NormalVectorImage )
        {
        m_NormalVectorImage = NormalVectorImageType::New();
        this->AllocateSpaceForImage( m_NormalVectorImage, output );
        computeNormalVectorImage = true;
        }
      if( !m_WeightImage )
        {
        m_WeightImage = WeightImageType::New();
        this->AllocateSpaceForImage( m_WeightImage, output );
        computeWeightImage = true;
        }

      // Compute the border normals and the weighting factor w
      // Normals are dependent on the border geometry in the fixed image so this
      // has to be completed only once.
      this->ComputeNormalVectorAndWeightImages( computeNormalVectorImage,
                                                computeWeightImage );
      }
    }

  // Compute the diffusion tensor image
  // The diffusion tensors are dependent on the normals computed in the
  // previous line, so this has to be completed only once.
  if( this->GetComputeRegularizationTerm() )
    {
    this->ComputeDiffusionTensorImage();
    }
}

/**
 * Initialize the state of the filter and equation before each iteration.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeIteration()
{
  std::cout << "Iteration #" << this->GetElapsedIterations() << std::endl;
  std::cout << "\tInitializeIteration for FILTER" << std::endl;

  if ( !this->GetFixedImage() || !this->GetMovingImage()
    || !this->GetDeformationField() )
    {
    itkExceptionMacro( << "FixedImage, MovingImage and/or DeformationField "
                       << "not set");
    }

  if( this->GetComputeRegularizationTerm() )
    {
    if ( this->GetUseAnisotropicRegularization()
         && ( !this->GetNormalVectorImage() || !this->GetWeightImage() ) )
      {
      itkExceptionMacro( << "NormalVector image and/or WeightImage not set");
      }

    // Update the deformation field component images
    // This depends on the current deformation field u, so it must be computed
    // on every iteration of the filter.
    this->UpdateDeformationVectorComponentImages();
    }

  // Update the function's deformation field
  this->GetRegistrationFunctionPointer()->SetDeformationField(
      this->GetDeformationField() );

  // Call the superclass implementation
  Superclass::InitializeIteration();

}

/**
 * Updates the border normals and the weighting factor w
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeNormalVectorAndWeightImages(
    bool computeNormalVectorImage, bool computeWeightImage )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetUseAnisotropicRegularization() );
  assert( m_BorderNormalsSurface );

  // Get the normals from the polydata
  vtkSmartPointer< vtkDataArray > normalData
      = m_BorderNormalsSurface->GetPointData()->GetNormals();
  vtkPointLocator* pointLocator = vtkPointLocator::New();
  pointLocator->SetDataSet( m_BorderNormalsSurface );

  // Iterate over the normal vector image and insert the normal of the closest
  // point
  NormalVectorImageIteratorType normalVectorIt(
      m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );

  // Iterate over the weight image and compute weight as a function of the
  // distance to the border
  WeightImageIteratorType weightIt(
      m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );

  itk::Index< ImageDimension >          imageIndex;
  itk::Point< double, ImageDimension >  imageCoordAsPoint;
  imageCoordAsPoint.Fill( 0 );
  double                                imageCoord[3] = {0, 0, 0};
  double                                borderCoord[3] = {0, 0, 0};
  vtkIdType                             id;
  WeightType                            distance;
  WeightType                            weight;
  NormalVectorType                      normal;

  std::cout << "Computing normals and weights... " << std::endl;

  // Determine the normals of and the distances to the nearest border
  for( normalVectorIt.GoToBegin(), weightIt.GoToBegin();
       !normalVectorIt.IsAtEnd();
       ++normalVectorIt, ++weightIt )
    {
    imageIndex = normalVectorIt.GetIndex();
    m_NormalVectorImage->TransformIndexToPhysicalPoint(
        imageIndex, imageCoordAsPoint );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      imageCoord[i] = imageCoordAsPoint[i];
      }
    id = pointLocator->FindClosestPoint( imageCoord );
    normal = normalData->GetTuple( id );

    if( computeNormalVectorImage )
      {
      normalVectorIt.Set( normal );
      }

    // Calculate distance between the current coord and the border surface coord
    m_BorderNormalsSurface->GetPoint( id, borderCoord );
    distance = 0.0;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      distance += pow( imageCoord[i] - borderCoord[i], 2 );
      }
    distance = sqrt( distance );

    // For now, put the distance in the weight image
    if( computeWeightImage )
      {
      weightIt.Set( distance );
      }
    }

  // Clean up memory
  pointLocator->Delete();

  if( computeNormalVectorImage )
    {
    //  // Smooth the normals to handle corners (because we are choosing the
    //  // closest point in the polydata
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

  if( computeWeightImage )
    {
    // Smooth the distance image to avoid "streaks" from faces of the polydata
    // (because we are choosing the closest point in the polydata)
    typedef itk::SmoothingRecursiveGaussianImageFilter
        < WeightImageType, WeightImageType > WeightSmoothingFilterType;
    typename WeightSmoothingFilterType::Pointer weightSmooth
        = WeightSmoothingFilterType::New();
    weightSmooth->SetInput( m_WeightImage );
    weightSmooth->SetSigma( 1.0 );
    weightSmooth->Update();
    m_WeightImage = weightSmooth->GetOutput();

    WeightImageIteratorType weightIt2(
        m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );

    // Iterate through the weight image and compute the weight from the distance
    for( weightIt2.GoToBegin(); !weightIt2.IsAtEnd(); ++weightIt2 )
      {
      weight = this->ComputeWeightFromDistance( weightIt2.Get() );
      weightIt2.Set( weight );
      }
    }

  std::cout << "Finished computing normals and weights." << std::endl;
}

/**
 * Updates the diffusion tensor image before each iteration
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
 * Updates the diffusion tensor image before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorImage()
{
  assert( this->GetComputeRegularizationTerm() );

  // Used to compute the tangential and normal diffusion tensor images
  // tangential:
  // P = I - wnn^T
  // tangentialMatrix = tangentialD = P^TP
  // normal:
  // normalMatrix = normalD = wnn^T

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
  NormalVectorImageIteratorType normalVectorIt;
  WeightImageIteratorType weightIt;
  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
      DiffusionTensorImageIteratorType;
  DiffusionTensorImageIteratorType tangentialDiffusionTensorIt;
  DiffusionTensorImageIteratorType normalDiffusionTensorIt;

  tangentialDiffusionTensorIt = DiffusionTensorImageIteratorType(
      m_TangentialDiffusionTensorImage,
      m_TangentialDiffusionTensorImage->GetLargestPossibleRegion() );
  if( this->GetUseAnisotropicRegularization() )
    {
    normalVectorIt = NormalVectorImageIteratorType(
        m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );
    weightIt = WeightImageIteratorType(
        m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );
    normalDiffusionTensorIt = DiffusionTensorImageIteratorType(
        m_NormalDiffusionTensorImage,
        m_NormalDiffusionTensorImage->GetLargestPossibleRegion() );
    }

  if( this->GetUseAnisotropicRegularization() )
    {
    normalVectorIt.GoToBegin();
    weightIt.GoToBegin();
    normalDiffusionTensorIt.GoToBegin();
    }

  for( tangentialDiffusionTensorIt.GoToBegin();
      !tangentialDiffusionTensorIt.IsAtEnd();
      ++ tangentialDiffusionTensorIt )
    {
    // Compute the tangential and normal diffusion tensor images
    if( !this->GetUseAnisotropicRegularization() )
      {
        // This is the diffusive (Gaussian) regularization
        tangentialMatrix.SetIdentity();
      }
    else
      {
      n = normalVectorIt.Get();
      w = weightIt.Get();

      // Create the normalMatrix used to calculate nn^T
      // (The first column is filled with the values of n, the rest are 0s)
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        normalMatrix( i,0 ) = n[i];
        for ( unsigned int j = 1; j < ImageDimension; j++ )
          {
          normalMatrix( i,j ) = 0;
          }
        }

      normalMatrix = normalMatrix * normalMatrix.GetTranspose(); // nn^T
      normalMatrix = normalMatrix * w; // wnn^T
      P.SetIdentity();
      P = P - normalMatrix; // I - wnn^T
      tangentialMatrix = P.GetTranspose();
      tangentialMatrix = tangentialMatrix * P; // P^TP
      }

    // Copy the itk::Matrix to the tensor
    // TODO there should be a better way to do this
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        tangentialDiffusionTensor( i,j ) = tangentialMatrix( i,j );
        }
      }
    if( this->GetUseAnisotropicRegularization() )
      {
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          normalDiffusionTensor( i,j ) = normalMatrix( i,j );
          }
        }
      }
    // Copy the diffusion tensors to m_TangentialDiffusionTensorImage and
    // m_NormalDiffusionTensorImage
    tangentialDiffusionTensorIt.Set( tangentialDiffusionTensor );
    if( this->GetUseAnisotropicRegularization() )
      {
      normalDiffusionTensorIt.Set( normalDiffusionTensor );
      }

    if( this->GetUseAnisotropicRegularization() )
      {
      ++normalVectorIt;
      ++weightIt;
      ++normalDiffusionTensorIt;
      }
    }
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

  if( this->GetUseAnisotropicRegularization() )
    {

    // Get the border normals
    NormalVectorImageIteratorType normalVectorIt(
        m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );

    // Get output (the current deformation field)
    typename OutputImageType::Pointer output = this->GetOutput();

    typedef itk::ImageRegionIterator< OutputImageType > OutputImageIteratorType;
    OutputImageIteratorType outputImageIt(
        output, output->GetLargestPossibleRegion() );

    // Extract normal components from output
    OutputImageIteratorType outputNormalImageIt(
        m_NormalDeformationField,
        m_NormalDeformationField->GetLargestPossibleRegion() );

    // Calculate the tangential and normal components of the deformation field
    NormalVectorType       n;
    DeformationVectorType  u; // deformation vector
    DeformationVectorType  normalDeformationVector;
    DeformationVectorType  tangentialDeformationVector;

    for( outputImageIt.GoToBegin(), normalVectorIt.GoToBegin(),
         outputNormalImageIt.GoToBegin();
         !outputImageIt.IsAtEnd(); ++outputImageIt )
      {
      n = normalVectorIt.Get();
      u = outputImageIt.Get();

      // normal component = (u^Tn)n
      normalDeformationVector = ( u * n ) * n;
      outputNormalImageIt.Set( normalDeformationVector );

      // tangential component = u - normal component
      tangentialDeformationVector = u - normalDeformationVector;

      // Assertion to test that the normal and tangential components were computed
      // corectly - they should be orthogonal
      if( normalDeformationVector * tangentialDeformationVector > 0.005 )
        {
        itkExceptionMacro( << "Normal and tangential deformation field "
                           << "components are not orthogonal" << std::endl
                           << "u = " << u[0] << " " << u[1] << " " << u[2]
                           << std::endl
                           << "n = " << n[0] << " " << n[1] << " " << n[2]
                           << std::endl
                           << "normal = " << normalDeformationVector[0]
                           << " " << normalDeformationVector[1] << " "
                           << normalDeformationVector[2] << std::endl
                           << "tangential = " << tangentialDeformationVector[0]
                           << " " << tangentialDeformationVector[1] << " "
                           << tangentialDeformationVector[2] << std::endl
                           << " dot product "
                           << normalDeformationVector
                                  * tangentialDeformationVector
                           << std::endl );
        }

      ++normalVectorIt;
      ++outputNormalImageIt;
      }
      m_NormalDeformationField->Modified();
    }

  // Updated the extracted components
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialComponentExtractor[i]->Update();
    if( this->GetUseAnisotropicRegularization() )
      {
      m_NormalComponentExtractor[i]->Update();
      }
    }
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
  int threadCount;
  TimeStepType dt;

  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = NumericTraits< TimeStepType >::Zero; // Not used during the
                                                      // calculate change step.
  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(
      this->CalculateChangeThreaderCallback, & str );

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  threadCount = this->GetMultiThreader()->GetNumberOfThreads();

  str.TimeStepList = new TimeStepType[threadCount];
  str.ValidTimeStepList = new bool[threadCount];
  for ( int i = 0; i < threadCount; ++i )
    {
    str.ValidTimeStepList[i] = false;
    }

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Resolve the single value time step to return
  dt = this->ResolveTimeStep(
      str.TimeStepList, str.ValidTimeStepList, threadCount);
  delete [] str.TimeStepList;
  delete [] str.ValidTimeStepList;

  // Explicitely call Modified on m_UpdateBuffer here
  // since ThreadedCalculateChange changes this buffer
  // through iterators which do not increment the
  // update buffer timestamp
  this->m_UpdateBuffer->Modified();

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
  DenseFDThreadStruct * str;
  int total, threadId, threadCount;

  threadId = (( MultiThreader::ThreadInfoStruct * )( arg ))->ThreadID;
  threadCount = (( MultiThreader::ThreadInfoStruct * )( arg ))->NumberOfThreads;

  str = ( DenseFDThreadStruct * )
            ((( MultiThreader::ThreadInfoStruct * )( arg ))->UserData );

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

  ThreadDeformationVectorComponentImageRegionType
      splitDeformationVectorComponentImageRegion;
  total = str->Filter->SplitRequestedRegion(
      threadId, threadCount, splitDeformationVectorComponentImageRegion);

  if (threadId < total)
    {
    str->TimeStepList[threadId]
        = str->Filter->ThreadedCalculateChange(
            splitRegion,
            splitNormalVectorImageRegion,
            splitDiffusionImageRegion,
            splitDeformationVectorComponentImageRegion,
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
typename AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChange(
          const ThreadRegionType &,
          int)
{
  // This function should never be called!
  itkExceptionMacro( << "ThreadedCalculateChange(regionToProcess, threadId) "
                     << "should never be called.  Use the other "
                     << "ThreadedCalculateChange function" );
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
    const ThreadDiffusionTensorImageRegionType &  diffusionRegionToProcess,
    const ThreadDeformationVectorComponentImageRegionType &
        deformationComponentRegionToProcess,
    int)
{
  typename OutputImageType::Pointer output = this->GetOutput();
  TimeStepType timeStep;
  void *globalData;

  // Get the FiniteDifferenceFunction to use in calculations.
  const RegistrationFunctionPointer df = this->GetRegistrationFunctionPointer();

  const typename OutputImageType::SizeType radius = df->GetRadius();

  bool computeRegularization = this->GetComputeRegularizationTerm();
  bool useAnisotropic = this->GetUseAnisotropicRegularization();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because input has been copied to output.

  // Define the boundary faces typedefs
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < OutputImageType > OutputImageFaceCalculatorType;
  typedef typename OutputImageFaceCalculatorType::FaceListType
      OutputImageFaceListType;
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < NormalVectorImageType > NormalVectorImageFaceCalculatorType;
  typedef typename NormalVectorImageFaceCalculatorType::FaceListType
      NormalVectorImageFaceListType;
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < DiffusionTensorImageType > DiffusionTensorImageFaceCalculatorType;
  typedef typename DiffusionTensorImageFaceCalculatorType::FaceListType
      DiffusionTensorImageFaceListType;
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < DeformationVectorComponentImageType >
      DeformationVectorComponentImageFaceCalculatorType;
  typedef typename DeformationVectorComponentImageFaceCalculatorType
      ::FaceListType
      DeformationVectorComponentImageFaceListType;
  typedef typename DeformationVectorComponentImageFaceListType::iterator
      DeformationVectorComponentImageFaceListIterator;

  // Setup the boundary faces for the output deformation field
  OutputImageFaceCalculatorType               outputImageFaceCalculator;
  OutputImageFaceListType                     outputImageFaceList;
  typename OutputImageFaceListType::iterator  outputImagefIt;

  // Setup the boundary faces for the normal vector images
  NormalVectorImageFaceCalculatorType
      normalVectorImageFaceCalculator;
  NormalVectorImageFaceListType                     normalVectorImageFaceList;
  typename NormalVectorImageFaceListType::iterator  normalVectorImagefIt;

  // Setup the boundary faces for the diffusion tensor images
  DiffusionTensorImageFaceCalculatorType    diffusionTensorImageFaceCalculator;
  DiffusionTensorImageFaceListType          tangentialDiffusionTensorFaceList;
  DiffusionTensorImageFaceListType          normalDiffusionTensorFaceList;
  typename DiffusionTensorImageFaceListType::iterator
      tangentialDiffusionTensorfIt;
  typename DiffusionTensorImageFaceListType::iterator
      normalDiffusionTensorfIt;

  // Setup the boundary faces for the deformation field component images
  DeformationVectorComponentImageFaceCalculatorType
      deformationVectorComponentImageFaceCalculator;
  itk::FixedArray
      < DeformationVectorComponentImageFaceListType, ImageDimension >
      deformationVectorTangentialComponentImageFaceListArray;
  itk::FixedArray
      < DeformationVectorComponentImageFaceListIterator, ImageDimension >
      deformationVectorTangentialComponentImagefItArray;
  itk::FixedArray
      < DeformationVectorComponentImageFaceListType, ImageDimension >
      deformationVectorNormalComponentImageFaceListArray;
  itk::FixedArray
      < DeformationVectorComponentImageFaceListIterator, ImageDimension >
      deformationVectorNormalComponentImagefItArray;

  // Actually initialize the face calculators and face list iterators
  outputImageFaceList = outputImageFaceCalculator(
      output, regionToProcess, radius );
  outputImagefIt = outputImageFaceList.begin();

  if( computeRegularization )
    {
    tangentialDiffusionTensorFaceList = diffusionTensorImageFaceCalculator(
        m_TangentialDiffusionTensorImage, diffusionRegionToProcess, radius );
    tangentialDiffusionTensorfIt = tangentialDiffusionTensorFaceList.begin();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      deformationVectorTangentialComponentImageFaceListArray[i]
          = deformationVectorComponentImageFaceCalculator(
              m_DeformationVectorTangentialComponents[i],
              deformationComponentRegionToProcess,
              radius );
      deformationVectorTangentialComponentImagefItArray[i]
          = deformationVectorTangentialComponentImageFaceListArray[i].begin();
      }

    if( useAnisotropic )
      {
      normalVectorImageFaceList = normalVectorImageFaceCalculator(
          m_NormalVectorImage, normalVectorRegionToProcess, radius );
      normalVectorImagefIt = normalVectorImageFaceList.begin();
      normalDiffusionTensorFaceList = diffusionTensorImageFaceCalculator(
          m_NormalDiffusionTensorImage, diffusionRegionToProcess, radius );
      normalDiffusionTensorfIt = normalDiffusionTensorFaceList.begin();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationVectorNormalComponentImageFaceListArray[i]
            = deformationVectorComponentImageFaceCalculator(
                m_DeformationVectorNormalComponents[i],
                deformationComponentRegionToProcess,
                radius );
        deformationVectorNormalComponentImagefItArray[i]
            = deformationVectorNormalComponentImageFaceListArray[i].begin();
        }
      }
    }

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();

  // Define the neighborhood iterator typedefs
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodIteratorType;
  typedef ImageRegionIterator< UpdateBufferType > UpdateIteratorType;
  typedef typename RegistrationFunctionType::
      NormalVectorImageNeighborhoodIteratorType
      NormalVectorImageNeighborhoodIteratorType;
  typedef typename RegistrationFunctionType::
      DiffusionTensorNeighborhoodIteratorType
      DiffusionTensorNeighborhoodIteratorType;
  typedef typename RegistrationFunctionType::
      DeformationVectorComponentNeighborhoodIteratorType
      DeformationVectorComponentNeighborhoodIteratorType;
  typedef typename RegistrationFunctionType::
      DeformationVectorComponentNeighborhoodIteratorArrayType
      DeformationVectorComponentNeighborhoodIteratorArrayType;

  // Process the boundary and non-boundary regions
  NeighborhoodIteratorType                    outputImageNeighborhoodIt;
  UpdateIteratorType                          updateIt;
  NormalVectorImageNeighborhoodIteratorType   normalVectorImageNeighborhoodIt;
  DiffusionTensorNeighborhoodIteratorType
      tangentialDiffusionTensorImageNeighborhoodIt;
  DiffusionTensorNeighborhoodIteratorType
      normalDiffusionTensorImageNeighborhoodIt;
  DeformationVectorComponentNeighborhoodIteratorArrayType
      deformationVectorTangentialComponentNeighborhoodItArray;
  DeformationVectorComponentNeighborhoodIteratorArrayType
      deformationVectorNormalComponentNeighborhoodItArray;

  for( ; outputImagefIt != outputImageFaceList.end(); ++outputImagefIt )
    {
    // Set the neighborhood iterators to the current face
    outputImageNeighborhoodIt = NeighborhoodIteratorType(
        radius, output, *outputImagefIt );
    updateIt = UpdateIteratorType( m_UpdateBuffer, *outputImagefIt );
    if( computeRegularization )
      {
      tangentialDiffusionTensorImageNeighborhoodIt
          = DiffusionTensorNeighborhoodIteratorType(
              radius, m_TangentialDiffusionTensorImage,
              *tangentialDiffusionTensorfIt );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationVectorTangentialComponentNeighborhoodItArray[i]
            = DeformationVectorComponentNeighborhoodIteratorType(
                radius, m_DeformationVectorTangentialComponents[i],
                *deformationVectorTangentialComponentImagefItArray[i] );
        }

      if( useAnisotropic )
        {
        normalVectorImageNeighborhoodIt
            = NormalVectorImageNeighborhoodIteratorType(
            radius, m_NormalVectorImage, *normalVectorImagefIt );
        normalDiffusionTensorImageNeighborhoodIt
            = DiffusionTensorNeighborhoodIteratorType(
                radius, m_NormalDiffusionTensorImage,
                *normalDiffusionTensorfIt );
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          deformationVectorNormalComponentNeighborhoodItArray[i]
              = DeformationVectorComponentNeighborhoodIteratorType(
                  radius, m_DeformationVectorNormalComponents[i],
                  *deformationVectorNormalComponentImagefItArray[i] );
          }
        }
      }

    // Go to the beginning of the neighborhood for this face
    outputImageNeighborhoodIt.GoToBegin();
    updateIt.GoToBegin();
    if( computeRegularization )
      {
      tangentialDiffusionTensorImageNeighborhoodIt.GoToBegin();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationVectorTangentialComponentNeighborhoodItArray[i].GoToBegin();
        }
      if( useAnisotropic )
        {
        normalVectorImageNeighborhoodIt.GoToBegin();
        normalDiffusionTensorImageNeighborhoodIt.GoToBegin();
        for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          deformationVectorNormalComponentNeighborhoodItArray[i].GoToBegin();
          }
        }
      }

    // Iterate through the neighborhood for this face and compute updates
    while ( !outputImageNeighborhoodIt.IsAtEnd() )
      {
      updateIt.Value() = df->ComputeUpdate(
          outputImageNeighborhoodIt,
          normalVectorImageNeighborhoodIt,
          tangentialDiffusionTensorImageNeighborhoodIt,
          deformationVectorTangentialComponentNeighborhoodItArray,
          normalDiffusionTensorImageNeighborhoodIt,
          deformationVectorNormalComponentNeighborhoodItArray,
          globalData);
      ++outputImageNeighborhoodIt;
      ++updateIt;
      if( computeRegularization )
        {
        ++tangentialDiffusionTensorImageNeighborhoodIt;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          ++deformationVectorTangentialComponentNeighborhoodItArray[i];
          }
        if( useAnisotropic )
          {
          ++normalVectorImageNeighborhoodIt;
          ++normalDiffusionTensorImageNeighborhoodIt;
          for( unsigned int i = 0; i < ImageDimension; i++ )
            {
            ++deformationVectorNormalComponentNeighborhoodItArray[i];
            }
          }
        }
      }

    // Go to the next face
    if( computeRegularization )
      {
      ++tangentialDiffusionTensorfIt;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        ++deformationVectorTangentialComponentImagefItArray[i];
        }
      if( useAnisotropic )
        {
        ++normalVectorImagefIt;
        ++normalDiffusionTensorfIt;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          ++deformationVectorNormalComponentImagefItArray[i];
          }
        }
      }
    }

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
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdate(TimeStepType dt)
{
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(this->ApplyUpdateThreaderCallback,
                                            &str);
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
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdateThreaderCallback( void * arg )
{
  DenseFDThreadStruct * str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (DenseFDThreadStruct *)
            (((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
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
 * Does the actual work of updating the output from the UpdateContainer over an
 * output region supplied by the multithreading mechanism.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::ThreadedApplyUpdate(TimeStepType dt,
                      const ThreadRegionType &regionToProcess,
                      int)
{
  ImageRegionIterator< UpdateBufferType > u(m_UpdateBuffer,   regionToProcess );
  ImageRegionIterator< OutputImageType > o(this->GetOutput(), regionToProcess );

  u = u.Begin();
  o = o.Begin();

  while ( !u.IsAtEnd() )
    {
    o.Value() += static_cast< DeformationVectorType >( u.Value() * dt );
                                                    // no adaptor support here

    ++o;
    ++u;
    }
}

} // end namespace itk

#endif
