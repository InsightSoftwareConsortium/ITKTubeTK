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
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "vtkDataArray.h"
#include "vtkPointData.h"

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

  m_UseDiffusiveRegularization      = true;
  m_BorderSurface                   = 0;
  m_BorderNormalsSurface            = 0;
  m_NormalVectorImage               = NormalVectorImageType::New();
  m_WeightImage                     = WeightImageType::New();
  m_TangentialDeformationField           = OutputImageType::New();
  m_NormalDeformationField               = OutputImageType::New();
  m_TangentialDiffusionTensorImage  = DiffusionTensorImageType::New();
  m_NormalDiffusionTensorImage      = DiffusionTensorImageType::New();

  typename RegistrationFunctionType::Pointer registrationFunction =
                                                RegistrationFunctionType::New();
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                          registrationFunction.GetPointer() ) );


  // Lambda for exponential decay used to calculate weight from distance.
  m_lambda = -0.01;

  // Setup the vtkPolyDataNormals to extract the normals from the surface
  m_BorderNormalsSurfaceFilter = BorderNormalsSurfaceFilterType::New();
  m_BorderNormalsSurfaceFilter->ComputePointNormalsOn();
  m_BorderNormalsSurfaceFilter->ComputeCellNormalsOff();
  //normalExtractor->SetFeatureAngle(30);
  // NOTE: default settings compute point normals, not cell normals
  m_BorderNormalsSurface = m_BorderNormalsSurfaceFilter->GetOutput();

  // Setup the point locator to find closest point on the surface
  m_PointLocator = PointLocatorType::New();

  // Setup the component extractor to extract the components from the deformation
  // field
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialComponentExtractor[i] = SelectionCastImageFilterType::New();
    m_TangentialComponentExtractor[i]->SetInput( this->GetOutput() );
    m_TangentialComponentExtractor[i]->SetIndex( i );

    m_NormalComponentExtractor[i] = SelectionCastImageFilterType::New();
    m_NormalComponentExtractor[i]->SetInput( m_NormalDeformationField );
    m_NormalComponentExtractor[i]->SetIndex( i );
    }

  // Setup the deformation field component images
  for (unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_DeformationVectorTangentialComponents[i]
                              = DeformationVectorComponentImageType::New();
    m_DeformationVectorTangentialComponents[i]
                              = m_TangentialComponentExtractor[i]->GetOutput();
    m_DeformationVectorNormalComponents[i]
                              = DeformationVectorComponentImageType::New();
    m_DeformationVectorNormalComponents[i]
                              = m_NormalComponentExtractor[i]->GetOutput();
    }

  // By default, compute the intensity distance and regularization terms
  this->SetComputeIntensityDistanceTerm( true );
  this->SetComputeRegularizationTerm( true );

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
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Border Surface: " << m_BorderSurface;
  //m_BorderSurface->PrintSelf( os, indent )
  os << indent << "Border Normals Surface: " << m_BorderNormalsSurface;
  //m_BorderNormalsSurface->PrintSelf( os, indent );
  os << indent << "UseDiffusiveRegularization: " << m_UseDiffusiveRegularization;
}

/**
 * Set/Get the timestep
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::SetTimeStep( const TimeStepType &t )
{
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  df->SetTimeStep( t );
}

/**
 * Set/Get the timestep
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
const typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >

::TimeStepType&
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::GetTimeStep() const
{
  typename RegistrationFunctionType::Pointer df
                              = dynamic_cast< RegistrationFunctionType * >
                                ( this->GetDifferenceFunction().GetPointer() );
  return df->GetTimeStep();
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

  // Do not smooth the deformation field with the
  // PDEDeformableRegistrationFilter method if we are using our own
  // regularization term
  df->SetComputeRegularizationTerm( compute );
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
 * Set the border surface
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
  m_BorderNormalsSurfaceFilter->SetInput( m_BorderSurface );
}

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

  // Check the timestep for stability
  double minSpacing;
  if ( this->GetUseImageSpacing() )
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

  TimeStepType ts = this->GetTimeStep();

  if ( ts > ratio )
    {
    itkWarningMacro(<< std::endl << "Anisotropic diffusion unstable time step:"
                    << ts << std::endl << "Minimum stable time step "
                    << "for this image is "
                    << ratio );
    }

  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the deformation field component images
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {

    // Allocate the tangential components
    this->AllocateSpaceForImage( m_DeformationVectorTangentialComponents[i],
                                 output);

    // Allocate the normal components
    this->AllocateSpaceForImage( m_DeformationVectorNormalComponents[i],
                                 output );
    }

  // Also allocate the diffusion tensor image
  // The diffusion tensor image has the same size as the deformation field and
  // holds the diffusion tensor matrix at each pixel
  DeformationFieldPointer deformationField = this->GetDeformationField();

  // Allocate the image of normals and weight image
  this->AllocateSpaceForImage( m_NormalVectorImage,
                               output );
  this->AllocateSpaceForImage( m_WeightImage,
                               output );

  // Allocate the output image's tangential and normal images
  this->AllocateSpaceForImage( m_TangentialDeformationField,
                               output );
  this->AllocateSpaceForImage( m_NormalDeformationField,
                               output );

  // Allocate the tangential and normal diffusion tensor images
  this->AllocateSpaceForImage( m_TangentialDiffusionTensorImage,
                               output );
  this->AllocateSpaceForImage( m_NormalDiffusionTensorImage,
                               output );

  // Compute the border normals, weighting factors and diffusion tensors
  if( !this->GetBorderSurface() )
    {
    itkExceptionMacro( << "Cannot perform registration without a border surface"
                       << std::endl );
    }

  // Update the border normals
  m_BorderNormalsSurfaceFilter->Update();

  // Make sure we now have the normals
  if ( !this->GetBorderNormalsSurface() )
    {
    itkExceptionMacro( << "Error computing border normals" << std::endl );
    }

  // Compute the border normals and the weighting factor w
  // Normals are dependent on the border geometry in the fixed image so this
  // has to be completed only once.
  this->ComputeNormalVectorAndWeightImages();

  // Compute the diffusion tensor image
  // The diffusion tensors are dependent on the normals computed in the
  // previous line, so this has to be completed only once.
  this->ComputeDiffusionTensorImage();
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
  std::cout << "Iteration #" << this->GetElapsedIterations() << std::endl;
  std::cout << "\tInitializeIteration for FILTER" << std::endl;

  if ( !this->GetFixedImage() || !this->GetMovingImage()
        || !this->GetDeformationField() || !this->GetBorderSurface()
        || !this->GetBorderNormalsSurface() )
    {
    itkExceptionMacro( << "FixedImage, MovingImage, DeformationField, border "
                       << "surface and/or border normals surface not set");
    }

  // Update the deformation field component images
  // This depends on the current deformation field u, so it must be computed
  // on every iteration of the filter.
  this->UpdateDeformationVectorComponentImages();

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
::ComputeNormalVectorAndWeightImages()
{
  assert( m_BorderNormalsSurface );

  // Get the normals
  vtkSmartPointer< vtkDataArray > normalData
                      = m_BorderNormalsSurface->GetPointData()->GetNormals();

  // Iterate over the normal vector image and insert the normal of the closest
  // point
  NormalVectorImageIteratorType normalVectorIt(
                              m_NormalVectorImage,
                              m_NormalVectorImage->GetLargestPossibleRegion() );
  // Iterate over the weight image and compute weight as a function of the
  // distance to the border
  WeightImageIteratorType weightIt( m_WeightImage,
                               m_WeightImage->GetLargestPossibleRegion() );

  m_PointLocator->SetDataSet( m_BorderNormalsSurface );

  itk::Index< ImageDimension > imageIndex;
  typename NormalVectorImageType::SpacingType spacing
                                            = m_NormalVectorImage->GetSpacing();
  typename NormalVectorImageType::PointType origin
                                            = m_NormalVectorImage->GetOrigin();
  itk::Point< double, ImageDimension > imageCoordAsPoint;
  imageCoordAsPoint.Fill( 0 );
  double imageCoord[3] = {0, 0, 0};
  double borderCoord[3] = {0, 0, 0};
  vtkIdType id;
  WeightType distance;
  WeightType weight;
  NormalVectorType normal;

  std::cout << "Computing normals and weights... " << std::endl;

  // Determine the normals of and the distances to the nearest border
  for( normalVectorIt.GoToBegin(), weightIt.GoToBegin();
       !normalVectorIt.IsAtEnd();
       ++normalVectorIt, ++weightIt )
    {
    imageIndex = normalVectorIt.GetIndex();

    m_NormalVectorImage->TransformIndexToPhysicalPoint( imageIndex,
                                                        imageCoordAsPoint );

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      imageCoord[i] = imageCoordAsPoint[i];
      }

    id = m_PointLocator->FindClosestPoint( imageCoord );

    normal = normalData->GetTuple( id );
//    for( int i = 0; i < ImageDimension; i++ )
//      {
//      normal[i] = fabs( normal[i] );
//      }

    normalVectorIt.Set( normal );

    // Calculate distance between the current coord and the border surface coord
    m_BorderNormalsSurface->GetPoint( id, borderCoord );
    distance = 0.0;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      distance += pow( imageCoord[i] - borderCoord[i], 2 );
      }
    distance = sqrt( distance );

    // For now, put the distance in the weight image
    weightIt.Set( distance );
    }

//  // Smooth the normals to handle corners (because we are choosing the closest
//  // point in the polydata
//  typedef itk::RecursiveGaussianImageFilter< NormalVectorImageType,
//                                                      NormalVectorImageType >
//                                                      NormalSmoothingFilterType;
//  typename NormalSmoothingFilterType::Pointer normalSmooth
//                                            = NormalSmoothingFilterType::New();
//  normalSmooth->SetInput( m_NormalVectorImage );
//  double normalSigma = 3.0;
//  normalSmooth->SetSigma( normalSigma );
//  normalSmooth->Update();
//  m_NormalVectorImage = normalSmooth->GetOutput();

  // Smooth the distance image to avoid "streaks" from faces of the polydata
  // (because we are choosing the closest point in the polydata)
  typedef itk::SmoothingRecursiveGaussianImageFilter< WeightImageType,
                                                      WeightImageType >
                                                      WeightSmoothingFilterType;
  typename WeightSmoothingFilterType::Pointer weightSmooth
                                            = WeightSmoothingFilterType::New();
  double weightSigma = 1.0;
  weightSmooth->SetInput( m_WeightImage );
  weightSmooth->SetSigma( weightSigma );
  weightSmooth->Update();
  m_WeightImage = weightSmooth->GetOutput();

  WeightImageIteratorType weightIt2( m_WeightImage,
                                m_WeightImage->GetLargestPossibleRegion() );

  // Iterate through the weight image and compute the weight from the distance
  for( weightIt2.GoToBegin(); !weightIt2.IsAtEnd(); ++weightIt2 )
    {
    weight = this->ComputeWeightFromDistance( weightIt2.Get() );
    weightIt2.Set( weight );
    }

  std::cout << "Finished computing normals and weights." << std::endl;
}

/**
 * Updates the diffusion tensor image before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename ImageToImageDiffusiveDeformableRegistrationFilter
                                < TFixedImage, TMovingImage, TDeformationField >
::WeightType
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
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
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
::ComputeDiffusionTensorImage()
{
  typedef itk::Matrix< DeformationVectorComponentType,
                       ImageDimension, ImageDimension > MatrixType;

  NormalVectorType                n;
  WeightType                      w;

  // Used to compute the tangential and normal diffusion tensor images

  // tangential:
  // P = I - wnn^T
  // tangentialMatrix = tangentialD = P^TP

  // normal:
  // normalMatrix = normalD = wnn^T

  MatrixType                      normalMatrix;
  MatrixType                      P;
  MatrixType                      tangentialMatrix;
  typedef typename DiffusionTensorImageType::PixelType
                                              DiffusionTensorImagePixelType;
  DiffusionTensorImagePixelType   tangentialD;
  DiffusionTensorImagePixelType   normalD;



  NormalVectorImageIteratorType normalVectorIt(
                              m_NormalVectorImage,
                              m_NormalVectorImage->GetLargestPossibleRegion() );
  WeightImageIteratorType weightIt( m_WeightImage,
                               m_WeightImage->GetLargestPossibleRegion() );

  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
                                                    DiffusionTensorIteratorType;
  DiffusionTensorIteratorType tangentialIt( m_TangentialDiffusionTensorImage,
                 m_TangentialDiffusionTensorImage->GetLargestPossibleRegion() );
  DiffusionTensorIteratorType normalIt( m_NormalDiffusionTensorImage,
                 m_NormalDiffusionTensorImage->GetLargestPossibleRegion() );

  for( normalVectorIt.GoToBegin(), tangentialIt.GoToBegin(),
       normalIt.GoToBegin(), weightIt.GoToBegin();
        !tangentialIt.IsAtEnd();
        ++normalVectorIt, ++tangentialIt, ++normalIt, ++weightIt )
    {

    // 1.  Get the border normal n and the weighting factor w
    n = normalVectorIt.Get();
    if ( !m_UseDiffusiveRegularization )
      {
      w = ( DeformationVectorComponentType ) 0.0;
      }
    else
      {
      // Get w here
      w = weightIt.Get();
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
::UpdateDeformationVectorComponentImages()
{
  // Get the border normals
  NormalVectorImageIteratorType normalVectorIterator( m_NormalVectorImage,
                              m_NormalVectorImage->GetLargestPossibleRegion() );

  // Get output (the current deformation field)
  typename OutputImageType::Pointer output = this->GetOutput();

  typedef itk::ImageRegionIterator< OutputImageType > IteratorType;
  IteratorType outputImageIterator(
                          output,
                          output->GetLargestPossibleRegion() );

  // Extract tangential and normal components from output
  IteratorType outputTangentialImageIterator(
                          m_TangentialDeformationField,
                          m_TangentialDeformationField->GetLargestPossibleRegion() );
  IteratorType outputNormalImageIterator(
                          m_NormalDeformationField,
                          m_NormalDeformationField->GetLargestPossibleRegion() );

  // Calculate the tangential and normal components of the deformation field
  NormalVectorType            n;
  DeformationVectorType  u;
  DeformationVectorType  normalU;
  DeformationVectorType  tangentialU;

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
    if( normalU * tangentialU > 0.005 )
      {
      itkExceptionMacro( << "Normal and tangential deformation field components"
                         << " are not orthogonal" << std::endl
                         << "u = " << u[0] << " " << u[1] << " " << u[2] << std::endl
                         << "n = " << n[0] << " " << n[1] << " " << n[2] << std::endl
                         << "normal = " << normalU[0] << " " << normalU[1]
                         << " " << normalU[2] << std::endl
                         << "tangential = "
                         << tangentialU[0] << " " << tangentialU[1] << " "
                         << tangentialU[2]
                         << " dot product "
                         << normalU * tangentialU << std::endl );
      }

    ++normalVectorIterator;
    ++outputTangentialImageIterator;
    ++outputNormalImageIterator;
    }

  // Feed tangential and normal components to the component extractor
  // Extract the components
  m_TangentialDeformationField->Modified();
  m_NormalDeformationField->Modified();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialComponentExtractor[i]->Update();
    m_NormalComponentExtractor[i]->Update();
    }

  // Little test to make sure that the component extractor was setup
  // properly - because m_DeformationVectorTangentialComponents and
  // m_DeformationVectorNormalComponents will be the input to the
  // registration function's ComputeUpdate()
  typename DeformationVectorComponentImageType::IndexType index;
  index.Fill(0);
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    assert( m_DeformationVectorTangentialComponents[i]->GetPixel( index )
        == m_TangentialComponentExtractor[i]->GetOutput()->GetPixel( index ) );
    assert( m_DeformationVectorNormalComponents[i]->GetPixel( index )
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

  typedef typename RegistrationFunctionType::DefaultBoundaryConditionType
      DefaultBoundaryConditionType;
  typedef typename RegistrationFunctionType::NormalVectorImageNeighborhoodType
                                NormalVectorImageNeighborhoodType;
  typedef ConstNeighborhoodIterator< DiffusionTensorImageType,
                                DefaultBoundaryConditionType >
                                DiffusionTensorNeighborhoodType;
  typedef typename
        RegistrationFunctionType::DeformationFieldComponentNeighborhoodType
                                DeformationFieldComponentNeighborhoodType;
  typedef typename
        RegistrationFunctionType::DeformationFieldComponentNeighborhoodArrayType
                                DeformationFieldComponentNeighborhoodArrayType;



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

  // Setup the boundary faces for the deformation field

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
                    = normalVectorFaceCalculator(
                              m_NormalVectorImage,
                              normalVectorRegionToProcess,
                              radius );

  typename NormalVectorFaceListType::iterator
                  normalVectorFaceListIt = normalVectorFaceList.begin();

  // Setup the boundary faces for the deformation field component images
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
                              < DeformationVectorComponentImageType >
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
                        m_DeformationVectorTangentialComponents[i],
                        deformationComponentRegionToProcess,
                        radius );
    deformationFieldTangentialComponentFaceListIterator[i] =
                        deformationFieldTangentialComponentFaceList[i].begin();

    deformationFieldNormalComponentFaceList[i]
        = deformationFieldNormalComponentFaceCalculator(
                        m_DeformationVectorNormalComponents[i],
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
                      m_DeformationVectorTangentialComponents[i],
                      *deformationFieldTangentialComponentFaceListIterator[i]);

    normalDFC[i] = DeformationFieldComponentNeighborhoodType(
                      radius,
                      m_DeformationVectorNormalComponents[i],
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
    nU.Value() = df->ComputeUpdate(nD,            // output (deformation field)
                                   normalVectorN, // m_NormalVectorImage
                                   tangentialDTN, // m_TangentialDiffusionTensorImage
                                   tangentialDFC, // m_DeformationVectorTangentialComponents
                                   normalDTN,     // m_NormalDiffusionTensorImage
                                   normalDFC,     // m_DeformationVectorTangentialComponents
                                   globalData);   // global data

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

  NeighborhoodIteratorType bD;
  UpdateIteratorType   bU;
  NormalVectorImageNeighborhoodType bNormalVectorN;
  DiffusionTensorNeighborhoodType bTangentialDTN;
  DiffusionTensorNeighborhoodType bNormalDTN;
  DeformationFieldComponentNeighborhoodArrayType bTangentialDFC;
  DeformationFieldComponentNeighborhoodArrayType bNormalDFC;

  ++normalVectorFaceListIt;
  ++tangentialDfIt;
  ++normalDfIt;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    ++bTangentialDFC[i];
    ++bNormalDFC[i];
    }
  for ( ++fIt; fIt != faceList.end(); ++fIt )
    {
    bD = NeighborhoodIteratorType( radius, output, *fIt );
    bU = UpdateIteratorType( m_UpdateBuffer, *fIt );
    bNormalVectorN = NormalVectorImageNeighborhoodType( radius, m_NormalVectorImage,
                                                        *normalVectorFaceListIt );
    bTangentialDTN = DiffusionTensorNeighborhoodType( radius,
                                                     m_TangentialDiffusionTensorImage,
                                                     *tangentialDfIt );
    bNormalDTN = DiffusionTensorNeighborhoodType( radius,
                                                  m_NormalDiffusionTensorImage,
                                                  *normalDfIt);
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      bTangentialDFC[i] = DeformationFieldComponentNeighborhoodType(
          radius, m_DeformationVectorTangentialComponents[i],
          *deformationFieldTangentialComponentFaceListIterator[i] );
      bNormalDFC[i] = DeformationFieldComponentNeighborhoodType(
          radius, m_DeformationVectorNormalComponents[i],
          *deformationFieldNormalComponentFaceListIterator[i] );
      }

    bD.GoToBegin();
    bU.GoToBegin();
    bNormalVectorN.GoToBegin();
    bTangentialDTN.GoToBegin();
    bNormalDTN.GoToBegin();
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      bTangentialDFC[i].GoToBegin();
      bNormalDFC[i].GoToBegin();
      }
    while ( !bD.IsAtEnd() )
      {
      nU.Value() = df->ComputeUpdate(bD,            // output (deformation field)
                                     bNormalVectorN, // m_NormalVectorImage
                                     bTangentialDTN, // m_TangentialDiffusionTensorImage
                                     bTangentialDFC, // m_DeformationVectorTangentialComponents
                                     bNormalDTN,     // m_NormalDiffusionTensorImage
                                     bNormalDFC,     // m_DeformationVectorTangentialComponents
                                     globalData);   // global data
      ++bD;
      ++bU;
      ++bNormalVectorN;
      ++bTangentialDTN;
      ++bNormalDTN;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        ++bTangentialDFC[i];
        ++bNormalDFC[i];
        }
      }
    ++normalVectorFaceListIt;
    ++tangentialDfIt;
    ++normalDfIt;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      ++bTangentialDFC[i];
      ++bNormalDFC[i];
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
ImageToImageDiffusiveDeformableRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
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
    o.Value() += static_cast<DeformationVectorType>(u.Value() * dt);
                                                    // no adaptor support here

    ++o;
    ++u;
    }
}

} // end namespace itk

#endif
