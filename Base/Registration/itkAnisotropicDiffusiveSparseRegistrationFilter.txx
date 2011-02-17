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
#ifndef __itkAnisotropicDiffusiveSparseRegistrationFilter_txx
#define __itkAnisotropicDiffusiveSparseRegistrationFilter_txx

#include "itkAnisotropicDiffusiveSparseRegistrationFilter.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
AnisotropicDiffusiveSparseRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::AnisotropicDiffusiveSparseRegistrationFilter()
{
  // Initialize attributes to NULL
  m_BorderSurface                               = 0;
  m_NormalMatrixImage                           = 0;
  m_WeightImage                                 = 0;

  // Lambda for exponential decay used to calculate weight from distance
  m_lambda = -0.01;
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  if( m_BorderSurface )
    {
    os << indent << "Border surface:" << std::endl;
    m_BorderSurface->Print( os );
    }
  if( m_NormalMatrixImage )
    {
    os << indent << "Normal vector image:" << std::endl;
    m_NormalMatrixImage->Print( os, indent );
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
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeDeformationComponentAndDerivativeImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // Setup pointers to the deformation component images - we have the
  // TANGENTIAL component, which is the entire deformation field, and the
  // NORMAL component, which is the deformation vectors projected onto their
  // normals
  // There are four terms that share these two deformation component images
  this->SetDeformationComponentImage( SMOOTH_TANGENTIAL, this->GetOutput() );
  this->SetDeformationComponentImage( PROP_TANGENTIAL, this->GetOutput() );

  DeformationFieldPointer normalDeformationField = DeformationFieldType::New();
  this->AllocateSpaceForImage( normalDeformationField, output );
  this->SetDeformationComponentImage( SMOOTH_NORMAL, normalDeformationField );
  this->SetDeformationComponentImage( PROP_NORMAL, normalDeformationField );

  // Setup the first and second order deformation component image derivatives
  // The two TANGENTIAL and two NORMAL components share images.
  int termOrder[4] = { SMOOTH_TANGENTIAL,
                       PROP_TANGENTIAL,
                       SMOOTH_NORMAL,
                       PROP_NORMAL };
  int t = 0;
  ScalarDerivativeImagePointer firstOrder = 0;
  TensorDerivativeImagePointer secondOrder = 0;
  for( int i = 0; i < ImageDimension; i++ )
    {
    for( int j = 0; j < this->GetNumberOfTerms(); j++ )
      {
      t = termOrder[j];
      if( t == SMOOTH_TANGENTIAL || t == SMOOTH_NORMAL )
        {
        firstOrder = ScalarDerivativeImageType::New();
        this->AllocateSpaceForImage( firstOrder, output );
        secondOrder = TensorDerivativeImageType::New();
        this->AllocateSpaceForImage( secondOrder, output );
        }
      this->SetDeformationComponentFirstOrderDerivative( t, i, firstOrder );
      this->SetDeformationComponentSecondOrderDerivative( t, i, secondOrder );
      }
    }
}

/**
 * All other initialization done before the initialize / calculate change /
 * apply update loop
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::SetupNormalMatrixAndWeightImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // If a normal vector image or weight image was supplied by the user, check
  // that it matches the output
  if( ( m_NormalMatrixImage
        && !this->CompareImageAttributes( m_NormalMatrixImage, output ) )
    || ( m_WeightImage
         && !this->CompareImageAttributes( m_WeightImage, output ) ) )
      {
    itkExceptionMacro( << "Normal matrix image and/or weight image must have "
                       << "the same attributes as the output deformation "
                       << "field" );
    }

  // Whether or not we must compute the normal vector and/or weight images
  bool computeNormals = !m_NormalMatrixImage;
  bool computeWeights = !m_WeightImage;

  // Compute the normal vector and/or weight images if required
  if( computeNormals || computeWeights )
    {
    // Ensure we have a border surface to work with
    if( !this->GetBorderSurface() )
      {
      itkExceptionMacro( << "You must provide a border surface, or both a "
                         << "normal matrix image and a weight image" );
      }

    // Compute the normals for the surface
    this->ComputeBorderSurfaceNormals();

    // Allocate the normal vector and/or weight images
    if( computeNormals )
      {
      m_NormalMatrixImage = NormalMatrixImageType::New();
      this->AllocateSpaceForImage( m_NormalMatrixImage, output );
      }
    if( computeWeights )
      {
      m_WeightImage = WeightImageType::New();
      this->AllocateSpaceForImage( m_WeightImage, output );
      }

    // Actually compute the normal vectors and/or weights
    this->ComputeNormalMatrixAndWeightImages( computeNormals, computeWeights );
    }
}

/**
 * Compute the normals for the border surface
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
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
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeNormalMatrixAndWeightImages( bool computeNormals, bool computeWeights )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( m_BorderSurface->GetPointData()->GetNormals() );
  assert( m_NormalMatrixImage );
  assert( m_WeightImage );

  std::cout << "Computing normals and weights... " << std::endl;

  // Setup iterators over the normal vector and weight images
  NormalMatrixImageRegionType normalIt(
      m_NormalMatrixImage, m_NormalMatrixImage->GetLargestPossibleRegion() );
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
  WeightComponentType                   distance = 0;
  NormalVectorType                      normalVector;
  normalVector.Fill(0);
  NormalMatrixType                      normalMatrix;
  normalMatrix.Fill(0);
  WeightComponentType                   weight = 0;
  WeightType                            weightMatrix;
  weightMatrix.Fill(0);

  // Determine the normals of and the distances to the nearest border
  for( normalIt.GoToBegin(), weightIt.GoToBegin();
       !normalIt.IsAtEnd();
       ++normalIt, ++weightIt )
    {
    normalMatrix.Fill(0);
    weightMatrix.Fill(0);

    // Find the normal of the surface point that is closest to the current voxel
    m_NormalMatrixImage->TransformIndexToPhysicalPoint( normalIt.GetIndex(),
                                                        imageCoordAsPoint );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      imageCoord[i] = imageCoordAsPoint[i];
      }
    id = pointLocator->FindClosestPoint( imageCoord );
    normalVector = normalData->GetTuple( id );
    if( computeNormals )
      {
      for( int i = 0; i < ImageDimension; i++ )
        {
        normalMatrix(i,0) = normalVector[i];
        }
      normalIt.Set( normalMatrix );
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
      weightMatrix(0,0) = distance;
      weightIt.Set( weightMatrix );
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
//    double weightSmoothingSigma = 1.0;
//    typedef itk::SmoothingRecursiveGaussianImageFilter
//        < WeightImageType, WeightImageType > WeightSmoothingFilterType;
//    typename WeightSmoothingFilterType::Pointer weightSmooth
//        = WeightSmoothingFilterType::New();
//    weightSmooth->SetInput( m_WeightImage );
//    weightSmooth->SetSigma( weightSmoothingSigma );
//    weightSmooth->Update();
//    m_WeightImage = weightSmooth->GetOutput();

    // Iterate through the weight image and compute the weight from the
    WeightImageRegionType weightIt(
        m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );
    for( weightIt.GoToBegin(); !weightIt.IsAtEnd(); ++weightIt )
      {
      weightMatrix = weightIt.Get();
      weight = this->ComputeWeightFromDistance( weightMatrix(0,0) );
      weightMatrix(0,0) = weight;
      weightIt.Set( weightMatrix );
      }
    }

  std::cout << "Finished computing normals and weights." << std::endl;
}

/**
 * Calculates the weighting between the anisotropic diffusive and diffusive
 * regularizations, based on a given distance from a voxel to the border
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::WeightComponentType
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeWeightFromDistance( const WeightComponentType distance ) const
{
  return exp( m_lambda * distance );
}

/**
 * Updates the diffusion tensor image before each run of the registration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorImages()
{
  assert( this->GetComputeRegularizationTerm() );

  // If required, allocate and compute the normal vector and weight images
  this->SetupNormalMatrixAndWeightImages();
  assert( m_NormalMatrixImage );
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

  NormalMatrixType        N;
  NormalMatrixType        N_transpose;
  WeightType              A;
  WeightComponentType     w;
  MatrixType              P;
  MatrixType              wP;   // wP
  MatrixType              wP_transpose;
  MatrixType              I_wP; // I-wP
  MatrixType              I_wP_transpose;

  MatrixType              smoothTangentialMatrix;  // (I-wP)^T * (I-wP)
  MatrixType              smoothNormalMatrix;      // (I-wP)^T * wP
  MatrixType              propTangentialMatrix;    // (wP)^T * (I-wP)
  MatrixType              propNormalMatrix;        // (wP)^T * wP

  DiffusionTensorType     smoothTangentialDiffusionTensor;
  DiffusionTensorType     smoothNormalDiffusionTensor;
  DiffusionTensorType     propTangentialDiffusionTensor;
  DiffusionTensorType     propNormalDiffusionTensor;

  // Setup iterators
  NormalMatrixImageRegionType normalIt = NormalMatrixImageRegionType(
      m_NormalMatrixImage, m_NormalMatrixImage->GetLargestPossibleRegion() );
  WeightImageRegionType weightIt = WeightImageRegionType(
      m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType smoothTangentialTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( SMOOTH_TANGENTIAL ),
          this->GetDiffusionTensorImage( SMOOTH_TANGENTIAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType smoothNormalTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( SMOOTH_NORMAL ),
          this->GetDiffusionTensorImage( SMOOTH_NORMAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType propTangentialTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( PROP_TANGENTIAL ),
          this->GetDiffusionTensorImage( PROP_TANGENTIAL )
          ->GetLargestPossibleRegion() );
  DiffusionTensorImageRegionType propNormalTensorIt
      = DiffusionTensorImageRegionType(
          this->GetDiffusionTensorImage( PROP_NORMAL ),
          this->GetDiffusionTensorImage( PROP_NORMAL )
          ->GetLargestPossibleRegion() );

  for( normalIt.GoToBegin(), weightIt.GoToBegin(),
       smoothTangentialTensorIt.GoToBegin(), smoothNormalTensorIt.GoToBegin(),
       propTangentialTensorIt.GoToBegin(), propNormalTensorIt.GoToBegin();
       !smoothTangentialTensorIt.IsAtEnd();
       ++normalIt, ++weightIt,
       ++smoothTangentialTensorIt, ++smoothNormalTensorIt,
       ++propTangentialTensorIt, ++propNormalTensorIt )
    {
    N = normalIt.Get();
    A = weightIt.Get();
    w = A(0,0);

    // The matrices are used for calculations, and will be copied to the
    // diffusion tensors afterwards.  The matrices are guaranteed to be
    // symmetric.
    P = N * A * N.GetTranspose();
    wP = P * w;
    I_wP.SetIdentity();
    I_wP = I_wP - wP;
    I_wP_transpose = I_wP.GetTranspose();
    smoothTangentialMatrix = I_wP_transpose * I_wP;
    smoothNormalMatrix = I_wP_transpose * wP;
    wP_transpose = wP.GetTranspose();
    propTangentialMatrix = wP_transpose * I_wP;
    propNormalMatrix = wP_transpose * wP;

    // Copy the matrices to the diffusion tensor
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        smoothTangentialDiffusionTensor(i,j) = smoothTangentialMatrix(i,j);
        smoothNormalDiffusionTensor(i,j) = smoothNormalMatrix(i,j);
        propTangentialDiffusionTensor(i,j) = propTangentialMatrix(i,j);
        propNormalDiffusionTensor(i,j) = propNormalMatrix(i,j);
        }
      }

    // Copy the diffusion tensors to their images
    smoothTangentialTensorIt.Set( smoothTangentialDiffusionTensor );
    smoothNormalTensorIt.Set( smoothNormalDiffusionTensor );
    propTangentialTensorIt.Set( propTangentialDiffusionTensor );
    propNormalTensorIt.Set( propNormalDiffusionTensor );
    }
}

/** Computes the multiplication vectors that the div(Tensor /grad u) values
 *  are multiplied by.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeMultiplicationVectorImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );
  assert( this->GetNormalMatrixImage() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // Allocate the images needed when using the anisotropic diffusive
  // regularization
  // There are no multiplication vectors for the smooth terms, and the
  // superclass will init them to zeros for us.
  // The prop multiplication vector is NAN_l

  // Iterate over the normal matrix and weight images
  NormalMatrixImageRegionType normalIt = NormalMatrixImageRegionType(
      m_NormalMatrixImage, m_NormalMatrixImage->GetLargestPossibleRegion() );
  WeightImageRegionType weightIt = WeightImageRegionType(
      m_WeightImage, m_WeightImage->GetLargestPossibleRegion() );

  DeformationVectorType multVector;
  multVector.Fill( 0.0 );
  NormalMatrixType N;
  N.Fill( 0.0 );
  WeightType A;
  A.Fill( 0.0 );
  DeformationVectorType N_l;
  N_l.Fill( 0.0 );

  for( int i = 0; i < ImageDimension; i++ )
    {
    // Create the multiplication vector image that is shared between the PROPs
    DeformationFieldPointer normalMultsImage = DeformationFieldType::New();
    this->AllocateSpaceForImage( normalMultsImage, output );

    // Calculate NAN_l
    DeformationVectorImageRegionType multIt = DeformationVectorImageRegionType(
        normalMultsImage, normalMultsImage->GetLargestPossibleRegion() );
    for( normalIt.GoToBegin(), weightIt.GoToBegin(), multIt.GoToBegin();
         !multIt.IsAtEnd();
         ++normalIt, ++weightIt, ++multIt )
           {
      multVector.Fill( 0.0 );
      N = normalIt.Get();
      A = weightIt.Get();
      for( int j = 0; j < ImageDimension; j++ )
        {
        N_l[j] = N[i][j];
        }
      multVector = N * A * N_l;
      multIt.Set( multVector );
      }

    // Set the multiplication vector image
    this->SetMultiplicationVectorImage( PROP_TANGENTIAL, i, normalMultsImage );
    this->SetMultiplicationVectorImage( PROP_NORMAL, i, normalMultsImage );
    }

  for( int i = 0; i < ImageDimension; i++ )
    {
    assert( this->GetMultiplicationVectorImage( SMOOTH_TANGENTIAL, i )
            == this->GetMultiplicationVectorImage( SMOOTH_NORMAL, i ) );
    assert( this->GetMultiplicationVectorImage( PROP_TANGENTIAL, i )
            == this->GetMultiplicationVectorImage( PROP_NORMAL, i ) );
    }
}

/**
 * Updates the deformation vector component images before each iteration
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationComponentImages()
{
  assert( this->GetComputeRegularizationTerm() );

  // The normal component of u_l is (NAN_l)^T * u
  // Conveniently, (NAN_l) is the prop multiplication vector (for both SMOOTH
  // and PROP), so we don't need to compute it again here
  DeformationVectorImageRegionArrayType NAN_lRegionArray;
  for( int i = 0; i < ImageDimension; i++ )
    {
    DeformationFieldPointer propMultImage
        = this->GetMultiplicationVectorImage( PROP_TANGENTIAL, i );
    assert( propMultImage );
    NAN_lRegionArray[i] = DeformationVectorImageRegionType(
        propMultImage, propMultImage->GetLargestPossibleRegion() );
    }

  // Setup the iterators for the deformation field
  OutputImagePointer output = this->GetOutput();
  OutputImageRegionType outputRegion(output,
                                     output->GetLargestPossibleRegion() );

  // We want to update the deformation component images for both SMOOTH_NORMAL
  // and PROP_NORMAL, but they point to the same image, so we will grab the
  // pointer from the first one to update both
  DeformationFieldPointer normalDeformationField
      = this->GetDeformationComponentImage( SMOOTH_NORMAL );
  OutputImageRegionType normalDeformationRegion(
      normalDeformationField,
      normalDeformationField->GetLargestPossibleRegion() );

  // Iterate over the deformation field and calculate the new normal components
  DeformationVectorType  u; // deformation vector
  u.Fill( 0.0 );
  DeformationVectorType  normalDeformationVector;
  normalDeformationVector.Fill( 0.0 );
  DeformationVectorType  tangentialDeformationVector;
  tangentialDeformationVector.Fill( 0.0 );

  outputRegion.GoToBegin();
  normalDeformationRegion.GoToBegin();
  for( int i = 0; i < ImageDimension; i++ )
    {
    NAN_lRegionArray[i].GoToBegin();
    }

  while( !outputRegion.IsAtEnd() )
    {
    u = outputRegion.Get();
    for( int i = 0; i < ImageDimension; i++ )
      {
      normalDeformationVector[i] = NAN_lRegionArray[i].Get() * u;
      }
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

    ++outputRegion;
    ++normalDeformationRegion;
    for( int i = 0; i < ImageDimension; i++ )
      {
      ++NAN_lRegionArray[i];
      }
    }
  normalDeformationField->Modified();

  assert( this->GetDeformationComponentImage( SMOOTH_TANGENTIAL )
          == this->GetDeformationComponentImage( PROP_TANGENTIAL ) );
  assert( this->GetDeformationComponentImage( SMOOTH_NORMAL )
          == this->GetDeformationComponentImage( PROP_NORMAL ) );
}

/**
 * Calculates the derivatives of the deformation vector derivatives after
 * each iteration.
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDeformationComponentDerivativeImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // Get the spacing and the radius
  SpacingType spacing = this->GetOutput()->GetSpacing();
  const RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  const typename OutputImageType::SizeType radius = df->GetRadius();

  // Calculate first and second order deformation component image derivatives.
  // The two TANGENTIAL and two NORMAL components share images.  We will
  // calculate derivatives only for the SMOOTH terms, and the PROP terms
  // will be automatically updated since they point to the same image data.
  DeformationComponentImageArrayType deformationComponentImageArray;
  deformationComponentImageArray.Fill( 0 );
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    if( i == SMOOTH_TANGENTIAL || i == SMOOTH_NORMAL )
      {
      this->ExtractXYZComponentsFromDeformationField(
          this->GetDeformationComponentImage(i),
          deformationComponentImageArray );

      for( int j = 0; j < ImageDimension; j++ )
        {
        this->ComputeDeformationComponentDerivativeImageHelper(
            deformationComponentImageArray[j], i, j, spacing, radius );
        }
      }
    }

  for( int i = 0; i < ImageDimension; i++ )
    {
    assert( this->GetDeformationComponentFirstOrderDerivative(
        SMOOTH_TANGENTIAL, i )
          == this->GetDeformationComponentFirstOrderDerivative(
              PROP_TANGENTIAL, i ) );
    assert( this->GetDeformationComponentFirstOrderDerivative(
        SMOOTH_NORMAL, i )
          == this->GetDeformationComponentFirstOrderDerivative(
              PROP_NORMAL, i ) );

    assert( this->GetDeformationComponentSecondOrderDerivative(
        SMOOTH_TANGENTIAL, i )
          == this->GetDeformationComponentSecondOrderDerivative(
              PROP_TANGENTIAL, i ) );
    assert( this->GetDeformationComponentSecondOrderDerivative(
        SMOOTH_NORMAL, i )
          == this->GetDeformationComponentSecondOrderDerivative(
              PROP_NORMAL, i ) );
    }
}


} // end namespace itk

#endif
