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
  m_UpdateBuffer = UpdateBufferType::New();

  // Initialize attributes to NULL
  m_BorderSurface                               = 0;
  m_NormalVectorImage                           = 0;
  m_WeightImage                                 = 0;
  m_NormalDeformationField                      = 0;
  m_TangentialDiffusionTensorImage              = 0;
  m_NormalDiffusionTensorImage                  = 0;
  m_TangentialDiffusionTensorDerivativeImage    = 0;
  m_NormalDiffusionTensorDerivativeImage        = 0;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_TangentialDeformationComponentImages[i]   = 0;
    m_NormalDeformationComponentImages[i]       = 0;
    }

  // Create the registration function
  typename RegistrationFunctionType::Pointer registrationFunction =
      RegistrationFunctionType::New();
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
      registrationFunction.GetPointer() ) );

  // Lambda for exponential decay used to calculate weight from distance
  m_lambda = -0.01;

  // By default, compute the intensity distance and regularization terms
  this->SetComputeRegularizationTerm( true );
  this->SetComputeIntensityDistanceTerm( true );
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
AnisotropicDiffusiveRegistrationFilter
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
AnisotropicDiffusiveRegistrationFilter
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
AnisotropicDiffusiveRegistrationFilter
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
  RegistrationFunctionPointer df = this->GetRegistrationFunctionPointer();
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

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the tangential diffusion tensor image and its derivative
  m_TangentialDiffusionTensorImage = DiffusionTensorImageType::New();
  this->AllocateSpaceForImage( m_TangentialDiffusionTensorImage,
                               output );

  m_TangentialDiffusionTensorDerivativeImage = TensorDerivativeImageType::New();
  this->AllocateSpaceForImage( m_TangentialDiffusionTensorDerivativeImage,
                               output );

  // Allocate the images needed when using the anisotropic diffusive
  // regularization
  if( this->GetUseAnisotropicRegularization() )
    {
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

  // Compute the diffusion tensors and their derivatives
  this->ComputeDiffusionTensorImages();
  this->ComputeDiffusionTensorDerivativeImages();
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
 * Update x, y, z components of the tangential and/or normal deformation field
 * components
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ExtractXYZFromDeformationComponents( bool extractTangentialComponents,
                                       bool extractNormalComponents )
{
  typename VectorIndexSelectionFilterType::Pointer indexSelector;

  if( extractTangentialComponents )
    {
    assert( this->GetOutput() );

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      indexSelector = VectorIndexSelectionFilterType::New();
      indexSelector->SetInput( this->GetOutput() );
      indexSelector->SetIndex( i );
      m_TangentialDeformationComponentImages[i] = indexSelector->GetOutput();
      indexSelector->Update();
      }
    }
  if( extractNormalComponents )
    {
    assert( m_NormalDeformationField );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      indexSelector = VectorIndexSelectionFilterType::New();
      indexSelector->SetInput( m_NormalDeformationField );
      indexSelector->SetIndex( i );
      m_NormalDeformationComponentImages[i] = indexSelector->GetOutput();
      indexSelector->Update();
      }
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
  assert( this->GetUseAnisotropicRegularization() );
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
  assert( m_TangentialDiffusionTensorImage );

  // If we are not using the anisotropic diffusive regularization, then we only
  // need to set the tangential diffusion tensors to the identity
  typename DiffusionTensorImageType::PixelType tangentialDiffusionTensor;
  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
      DiffusionTensorImageRegionType;
  DiffusionTensorImageRegionType tangentialTensorIt
      = DiffusionTensorImageRegionType(
          m_TangentialDiffusionTensorImage,
          m_TangentialDiffusionTensorImage->GetLargestPossibleRegion() );

  if( !this->GetUseAnisotropicRegularization() )
    {
    for( tangentialTensorIt.GoToBegin();
        !tangentialTensorIt.IsAtEnd();
        ++tangentialTensorIt )
          {
      tangentialDiffusionTensor.SetIdentity();
      tangentialTensorIt.Set( tangentialDiffusionTensor );
      }
    return;
    }

  // If we are using the anisotropic diffusive regularization, then we need
  // to setup the tangential diffusion tensors and the normal diffusion tensors
  assert( m_NormalVectorImage );
  assert( m_WeightImage );
  assert( m_NormalDiffusionTensorImage );

  // Used to compute the tangential and normal diffusion tensor images
  // tangential:
  // P = I - wnn^T
  // tangentialMatrix = tangentialDiffusionTensor = P^TP
  // normal:
  // normalMatrix = normalDiffusionTensor = wnn^T

  typedef itk::Matrix
      < DeformationVectorComponentType, ImageDimension, ImageDimension >
      MatrixType;

  NormalVectorType                              n;
  WeightType                                    w;
  MatrixType                                    P;
  MatrixType                                    normalMatrix;
  MatrixType                                    tangentialMatrix;
  typename DiffusionTensorImageType::PixelType  normalDiffusionTensor;

  // Setup iterators
  NormalVectorImageRegionType normalIt = NormalVectorImageRegionType(
      m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );
  WeightImageRegionType weightIt = WeightImageRegionType(
      m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );
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
  assert( this->GetComputeRegularizationTerm() );

  // Compute the diffusion tensor derivative image for the tangential image
  this->ComputeDiffusionTensorDerivativeImage(
      m_TangentialDiffusionTensorImage,
      m_TangentialDiffusionTensorDerivativeImage );

  // Compute the diffusion tensor derivative image for the normal plane
  if( this->GetUseAnisotropicRegularization() )
    {
    this->ComputeDiffusionTensorDerivativeImage(
        m_NormalDiffusionTensorImage, m_NormalDiffusionTensorDerivativeImage );
    }
}

/**
 * Actually computes the diffusion tensor derivative images
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorDerivativeImage(
    DiffusionTensorImagePointer tensorImage,
    TensorDerivativeImagePointer tensorDerivativeImage )
{
  assert( tensorImage );
  assert( tensorDerivativeImage );

  // Get the FiniteDifferenceFunction and radius to use in calculations.
  const RegistrationFunctionPointer df = this->GetRegistrationFunctionPointer();
  assert( df );
  const RegularizationFunctionPointer reg
      = df->GetRegularizationFunctionPointer();
  assert( reg );

  // Get the radius
  const typename OutputImageType::SizeType radius = df->GetRadius();

  // Setup the structs for the face calculations and their iterators
  FaceStruct< DiffusionTensorImagePointer > tensorStruct( tensorImage, radius );
  FaceStruct< TensorDerivativeImagePointer > tensorDerivativeStruct(
      tensorDerivativeImage, radius );

  // Iterator over the current face
  DiffusionTensorNeighborhoodType tensorNeighborhood;
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
      reg->ComputeDiffusionTensorFirstDerivative(tensorNeighborhood,
                                                 tensorDerivativeRegion );
      }
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
  Superclass::InitializeIteration();

  // Update the deformation field component images
  // Since the components depend on the current tangential and normal
  // deformation fields, they must be computed on every registration iteration
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
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationVectorComponentImages()
{
  assert( this->GetComputeRegularizationTerm() );

  if( this->GetUseAnisotropicRegularization() )
    {
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
    }

  // Update the extracted components
  this->ExtractXYZFromDeformationComponents(
      true, this->GetUseAnisotropicRegularization() );
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

  ThreadTensorDerivativeImageRegionType splitTensorDerivativeImageRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitTensorDerivativeImageRegion);

  ThreadDeformationVectorComponentImageRegionType
      splitDeformationVectorComponentImageRegion;
  total = str->Filter->SplitRequestedRegion(
      threadId, threadCount, splitDeformationVectorComponentImageRegion);

  if (threadId < total)
    {
    str->TimeStepList[threadId] = str->Filter->ThreadedCalculateChange(
      splitRegion,
      splitNormalVectorImageRegion,
      splitDiffusionImageRegion,
      splitTensorDerivativeImageRegion,
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
    const ThreadDiffusionTensorImageRegionType &  tensorRegionToProcess,
    const ThreadTensorDerivativeImageRegionType &derivativeRegionToProcess,
    const ThreadDeformationVectorComponentImageRegionType &
      deformationComponentRegionToProcess,
    int)
{
  typename OutputImageType::Pointer output = this->GetOutput();
  TimeStepType timeStep;
  void *globalData;

  // Get the FiniteDifferenceFunction to use in calculations.
  const RegistrationFunctionPointer df = this->GetRegistrationFunctionPointer();
  assert( df );

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
      < TensorDerivativeImageType > TensorDerivativeImageFaceCalculatorType;
  typedef typename TensorDerivativeImageFaceCalculatorType::FaceListType
      TensorDerivativeImageFaceListType;
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
  TensorDerivativeImageFaceCalculatorType   TensorDerivativeImageFaceCalculator;
  DiffusionTensorImageFaceListType          tangentialDiffusionTensorFaceList;
  TensorDerivativeImageFaceListType
      tangentialDiffusionTensorDerivativeFaceList;
  DiffusionTensorImageFaceListType          normalDiffusionTensorFaceList;
  TensorDerivativeImageFaceListType
      normalDiffusionTensorDerivativeFaceList;
  typename DiffusionTensorImageFaceListType::iterator
      tangentialDiffusionTensorfIt;
  typename TensorDerivativeImageFaceListType::iterator
      tangentialDiffusionTensorDerivativefIt;
  typename DiffusionTensorImageFaceListType::iterator
      normalDiffusionTensorfIt;
  typename TensorDerivativeImageFaceListType::iterator
      normalDiffusionTensorDerivativefIt;

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
        m_TangentialDiffusionTensorImage, tensorRegionToProcess, radius );
    tangentialDiffusionTensorfIt = tangentialDiffusionTensorFaceList.begin();
    tangentialDiffusionTensorDerivativeFaceList
        = TensorDerivativeImageFaceCalculator(
            m_TangentialDiffusionTensorDerivativeImage,
            tensorRegionToProcess,
            radius );
    tangentialDiffusionTensorDerivativefIt
        = tangentialDiffusionTensorDerivativeFaceList.begin();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      deformationVectorTangentialComponentImageFaceListArray[i]
          = deformationVectorComponentImageFaceCalculator(
              m_TangentialDeformationComponentImages[i],
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
          m_NormalDiffusionTensorImage, tensorRegionToProcess, radius );
      normalDiffusionTensorfIt = normalDiffusionTensorFaceList.begin();
      normalDiffusionTensorDerivativeFaceList
          = TensorDerivativeImageFaceCalculator(
              m_NormalDiffusionTensorDerivativeImage,
              tensorRegionToProcess,
              radius );
      normalDiffusionTensorDerivativefIt
          = normalDiffusionTensorDerivativeFaceList.begin();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationVectorNormalComponentImageFaceListArray[i]
            = deformationVectorComponentImageFaceCalculator(
                m_NormalDeformationComponentImages[i],
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
    NeighborhoodType;
  typedef ImageRegionIterator< UpdateBufferType >
    UpdateIteratorType;
  typedef typename RegistrationFunctionType::
    NormalVectorNeighborhoodType
    NormalVectorNeighborhoodType;

  typedef typename RegistrationFunctionType::
    TensorDerivativeImageRegionType
    TensorDerivativeImageRegionType;
  typedef typename RegistrationFunctionType::
    DeformationVectorComponentNeighborhoodType
    DeformationVectorComponentNeighborhoodType;
  typedef typename RegistrationFunctionType::
    DeformationVectorComponentNeighborhoodArrayType
    DeformationVectorComponentNeighborhoodArrayType;

  // Process the boundary and non-boundary regions
  NeighborhoodType                  outputImageNeighborhoodIt;
  UpdateIteratorType                        updateIt;
  NormalVectorNeighborhoodType normalVectorImageNeighborhoodIt;
  DiffusionTensorNeighborhoodType
    tangentialDiffusionTensorImageNeighborhoodIt;
  TensorDerivativeImageRegionType
    tangentialDiffusionTensorImageDerivativeNeighborhoodIt;
  DiffusionTensorNeighborhoodType
    normalDiffusionTensorImageNeighborhoodIt;
  TensorDerivativeImageRegionType
    normalDiffusionTensorImageDerivativeNeighborhoodIt;
  DeformationVectorComponentNeighborhoodArrayType
    deformationVectorTangentialComponentNeighborhoodItArray;
  DeformationVectorComponentNeighborhoodArrayType
    deformationVectorNormalComponentNeighborhoodItArray;

  for(; outputImagefIt != outputImageFaceList.end(); ++outputImagefIt )
    {
    // Set the neighborhood iterators to the current face
    outputImageNeighborhoodIt = NeighborhoodType(
        radius, output, *outputImagefIt );
    updateIt = UpdateIteratorType( m_UpdateBuffer, *outputImagefIt );
    if( computeRegularization )
      {
      tangentialDiffusionTensorImageNeighborhoodIt
        = DiffusionTensorNeighborhoodType(
          radius, m_TangentialDiffusionTensorImage,
          *tangentialDiffusionTensorfIt );
      tangentialDiffusionTensorImageDerivativeNeighborhoodIt
        = TensorDerivativeImageRegionType(
          m_TangentialDiffusionTensorDerivativeImage,
          *tangentialDiffusionTensorDerivativefIt );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationVectorTangentialComponentNeighborhoodItArray[i]
          = DeformationVectorComponentNeighborhoodType(
            radius, m_TangentialDeformationComponentImages[i],
            * deformationVectorTangentialComponentImagefItArray[i] );
        }

      if( useAnisotropic )
        {
        normalVectorImageNeighborhoodIt
          = NormalVectorNeighborhoodType(
          radius, m_NormalVectorImage, *normalVectorImagefIt );
        normalDiffusionTensorImageNeighborhoodIt
          = DiffusionTensorNeighborhoodType(
          radius, m_NormalDiffusionTensorImage,
          *normalDiffusionTensorfIt );
        normalDiffusionTensorImageDerivativeNeighborhoodIt
          = TensorDerivativeImageRegionType(
          m_NormalDiffusionTensorDerivativeImage,
          *normalDiffusionTensorDerivativefIt );
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          deformationVectorNormalComponentNeighborhoodItArray[i]
            = DeformationVectorComponentNeighborhoodType(
            radius, m_NormalDeformationComponentImages[i],
            * deformationVectorNormalComponentImagefItArray[i] );
          }
        }
      }

    // Go to the beginning of the neighborhood for this face
    outputImageNeighborhoodIt.GoToBegin();
    updateIt.GoToBegin();
    if( computeRegularization )
      {
      tangentialDiffusionTensorImageNeighborhoodIt.GoToBegin();
      tangentialDiffusionTensorImageDerivativeNeighborhoodIt.GoToBegin();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        deformationVectorTangentialComponentNeighborhoodItArray[i].
          GoToBegin();
        }
      if( useAnisotropic )
        {
        normalVectorImageNeighborhoodIt.GoToBegin();
        normalDiffusionTensorImageNeighborhoodIt.GoToBegin();
        normalDiffusionTensorImageDerivativeNeighborhoodIt.GoToBegin();
        for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          deformationVectorNormalComponentNeighborhoodItArray[i].
            GoToBegin();
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
        tangentialDiffusionTensorImageDerivativeNeighborhoodIt,
        deformationVectorTangentialComponentNeighborhoodItArray,
        normalDiffusionTensorImageNeighborhoodIt,
        normalDiffusionTensorImageDerivativeNeighborhoodIt,
        deformationVectorNormalComponentNeighborhoodItArray,
        globalData);
      ++outputImageNeighborhoodIt;
      ++updateIt;
      if( computeRegularization )
        {
        ++tangentialDiffusionTensorImageNeighborhoodIt;
        ++tangentialDiffusionTensorImageDerivativeNeighborhoodIt;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          ++deformationVectorTangentialComponentNeighborhoodItArray[i];
          }
        if( useAnisotropic )
          {
          ++normalVectorImageNeighborhoodIt;
          ++normalDiffusionTensorImageNeighborhoodIt;
          ++normalDiffusionTensorImageDerivativeNeighborhoodIt;
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
      ++tangentialDiffusionTensorDerivativefIt;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        ++deformationVectorTangentialComponentImagefItArray[i];
        }
      if( useAnisotropic )
        {
        ++normalVectorImagefIt;
        ++normalDiffusionTensorfIt;
        ++normalDiffusionTensorDerivativefIt;
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
  this->GetMultiThreader()->SetSingleMethod(
    this->ApplyUpdateThreaderCallback, &str);
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
    str->Filter->ThreadedApplyUpdate(str->TimeStep, splitRegion,
      threadId );
    }

  return ITK_THREAD_RETURN_VALUE;
}

/**
 * Does the actual work of updating the output from the UpdateContainer
 * over an output region supplied by the multithreading mechanism.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::ThreadedApplyUpdate(TimeStepType dt,
  const ThreadRegionType &regionToProcess, int )
{
  UpdateBufferRegionType  u(m_UpdateBuffer, regionToProcess );
  OutputImageRegionType   o(this->GetOutput(), regionToProcess );

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
