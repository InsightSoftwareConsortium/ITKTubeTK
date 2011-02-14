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
  os << indent << "lambda: " << m_lambda << std::endl;
}

/**
 * Setup the pointers for the deformation component images
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeDeformationComponentImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the images needed when using the anisotropic diffusive
  // regularization
  this->SetDeformationComponentImage( TANGENTIAL, this->GetOutput() );
  DeformationFieldPointer normalDeformationField = DeformationFieldType::New();
  this->AllocateSpaceForImage( normalDeformationField, output );
  this->SetDeformationComponentImage( NORMAL, normalDeformationField );
}

/**
 * All other initialization done before the initialize / calculate change /
 * apply update loop
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::SetupNormalVectorAndWeightImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

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
    this->ComputeNormalVectorAndWeightImages( computeNormals, computeWeights );
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
::ComputeWeightFromDistance( const WeightType distance ) const
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

  // If required, allocate and compute the normal vector and weight images
  this->SetupNormalVectorAndWeightImages();
  assert( m_NormalVectorImage );
  assert( m_WeightImage );

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
          this->GetDiffusionTensorImage( TANGENTIAL ),
          this->GetDiffusionTensorImage( TANGENTIAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType normalTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( NORMAL ),
          this->GetDiffusionTensorImage( NORMAL )->GetLargestPossibleRegion() );

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

/** Computes the multiplication vectors that the div(Tensor /grad u) values
 *  are multiplied by.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeMultiplicationVectorImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );
  assert( this->GetNormalVectorImage() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  // Allocate the images needed when using the anisotropic diffusive
  // regularization
  // There is no multiplication vector for the tangential term, and the
  // superclass will init it to zeros for us.
  // The normal multiplication vector is n_l*n

  // Iterate over the normal vector image
  typedef itk::ImageRegionConstIterator< NormalVectorImageType >
    NormalVectorImageConstRegionType;
  NormalVectorImageConstRegionType normalIt = NormalVectorImageConstRegionType(
      this->GetNormalVectorImage(),
      this->GetNormalVectorImage()->GetLargestPossibleRegion() );

  DeformationVectorImageArrayType normalMultsArray;
  DeformationVectorType normalVector;
  normalVector.Fill( 0.0 );
  for( int i = 0; i < ImageDimension; i++ )
    {
    // Create the multiplication vector image
    DeformationFieldPointer normalMultsImage = DeformationFieldType::New();
    this->AllocateSpaceForImage( normalMultsImage, output );

    // Calculate n_l*n
    DeformationVectorImageRegionType multIt = DeformationVectorImageRegionType(
        normalMultsImage, normalMultsImage->GetLargestPossibleRegion() );
    for( normalIt.GoToBegin(), multIt.GoToBegin();
         !multIt.IsAtEnd();
         ++normalIt, ++multIt )
           {
      normalVector = normalIt.Get();
      multIt.Set( normalVector[i] * normalVector );
      }

    // Set the multiplication vector image to the array
    normalMultsArray[i] = normalMultsImage;
    }
  this->SetMultiplicationVectorImage( NORMAL, normalMultsArray );
}

/**
 * Updates the deformation vector component images before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationComponentImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetNormalVectorImage() );
  assert( this->GetWeightImage() );

  // Setup iterators
  NormalVectorImageRegionType normalVectorRegion(
      m_NormalVectorImage, m_NormalVectorImage->GetLargestPossibleRegion() );

  typename OutputImageType::Pointer output = this->GetOutput();
  OutputImageRegionType outputRegion(output,
                                     output->GetLargestPossibleRegion() );

  DeformationFieldPointer normalDeformationField
      = this->GetDeformationComponentImage( NORMAL );
  OutputImageRegionType normalDeformationRegion(
      normalDeformationField,
      normalDeformationField->GetLargestPossibleRegion() );

  // Calculate the tangential and normal components of the deformation field
  NormalVectorType       n;
  DeformationVectorType  u; // deformation vector
  DeformationVectorType  normalDeformationVector;
  DeformationVectorType  tangentialDeformationVector;

  for( normalVectorRegion.GoToBegin(), outputRegion.GoToBegin(),
       normalDeformationRegion.GoToBegin();
  !outputRegion.IsAtEnd();
  ++normalVectorRegion, ++outputRegion, ++normalDeformationRegion )
    {
    n = normalVectorRegion.Get();
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
  normalDeformationField->Modified();
}

} // end namespace itk

#endif
