/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeAnisotropicHybridDiffusionImageFilter_hxx
#define __itktubeAnisotropicHybridDiffusionImageFilter_hxx


#include <itkFixedArray.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
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
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::AnisotropicHybridDiffusionImageFilter( void )
{
  m_ThresholdParameterC = 3.31488;
  m_ContrastParameterLambdaHybrid = 30.0;
  m_ContrastParameterLambdaCED = 30.0;
  m_ContrastParameterLambdaEED = 20.0;
  m_Sigma = 1.0;
  m_SigmaOuter = 1.0;
  m_Alpha = 0.001;
}

template< class TInputImage, class TOutputImage >
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::UpdateDiffusionTensorImage( void )
{
  itkDebugMacro( << "UpdateDiffusionTensorImage() called." );

  /* IN THIS METHOD, the following items will be implemented
   - Compute the structure tensor ( Multiscale version structure tensor )
   - Compute its eigenvectors
   - Compute eigenvalues corresponding to the diffusion matrix tensor
     ( Here is where all the magic happens for EED, CED and hybrid switch )
  */

  //Step 1: Compute the structure tensor and identify the eigenvectors
  //Step 1.1: Compute the structure tensor
  // Instantiate the structure tensor filter
  typename StructureTensorFilterType::Pointer
    structureTensorFilter  = StructureTensorFilterType::New();

  structureTensorFilter->SetInput( this->GetOutput() );
  structureTensorFilter->SetSigma( m_Sigma );
  structureTensorFilter->SetSigmaOuter( m_SigmaOuter );
  structureTensorFilter->Update();

  // Step 1.2: Identify the eigenvectors of the structure tensor
  typedef  Matrix< double, 3, 3>
    EigenVectorMatrixType;
  typedef  Image< EigenVectorMatrixType, 3>
    EigenVectorImageType;
  typedef  typename itk::Image< EigenValueArrayType, 3>
    EigenValueImageType;

  typedef  typename StructureTensorFilterType::OutputImageType
    SymmetricSecondRankTensorImageType;
  typedef  SymmetricEigenVectorAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType,
    EigenVectorImageType>
    EigenVectorAnalysisFilterType;

  typename EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter
    = EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( 3 );
  eigenVectorAnalysisFilter->OrderEigenValuesBy(
    EigenValueOrderEnum::OrderByValue );

  eigenVectorAnalysisFilter->SetInput( structureTensorFilter->GetOutput() );
  eigenVectorAnalysisFilter->Modified();
  eigenVectorAnalysisFilter->Update();

  //Step 1.3: Compute the eigenvalues
  typedef itk::SymmetricEigenAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType>
    EigenAnalysisFilterType;

  typename EigenAnalysisFilterType::Pointer eigenAnalysisFilter
    = EigenAnalysisFilterType::New();
  eigenAnalysisFilter->SetDimension( 3 );
  eigenAnalysisFilter->OrderEigenValuesBy(
    EigenValueOrderEnum::OrderByValue );

  eigenAnalysisFilter->SetInput( structureTensorFilter->GetOutput() );
  eigenAnalysisFilter->Update();

  /* Compute the gradient magnitude. This is required to set Lambda1 */

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<InputImageType>
    GradientMagnitudeFilterType;

  typename GradientMagnitudeFilterType::Pointer gradientMagnitudeFilter
    = GradientMagnitudeFilterType::New();
  gradientMagnitudeFilter->SetInput( this->GetInput() );
  gradientMagnitudeFilter->SetSigma( m_Sigma );
  gradientMagnitudeFilter->Update();

  /* Step 2: Generate the diffusion tensor matrix
      D = [v1 v2 v3] [DiagonalMatrixContainingLambdas] [v1 v2 v3]^t
  */

  //Setup the iterators
  //
  //Iterator for the eigenvector matrix image
  EigenVectorImageType::ConstPointer eigenVectorImage =
    eigenVectorAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenVectorImageType>
    eigenVectorImageIterator;
  eigenVectorImageIterator
    = itk::ImageRegionConstIterator<EigenVectorImageType>(
    eigenVectorImage, eigenVectorImage->GetRequestedRegion() );
  eigenVectorImageIterator.GoToBegin();

  //Iterator for the diffusion tensor image
  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
    DiffusionTensorIteratorType;
  DiffusionTensorIteratorType it( this->GetDiffusionTensorImage(),
    this->GetDiffusionTensorImage()->GetLargestPossibleRegion() );

  //Iterator for the eigenvalue image
  typename EigenValueImageType::ConstPointer eigenImage
    = eigenAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenValueImageType>
    eigenValueImageIterator;
  eigenValueImageIterator = itk::ImageRegionConstIterator<
    EigenValueImageType>( eigenImage, eigenImage->GetRequestedRegion() );

  //Iterator for the gradient magnitude image
  typedef typename GradientMagnitudeFilterType::OutputImageType
    GradientMagnitudeOutputImageType;
  typename GradientMagnitudeOutputImageType::Pointer
    gradientMagnitudeOutputImage = gradientMagnitudeFilter->GetOutput();

  itk::ImageRegionConstIterator<GradientMagnitudeOutputImageType>
    gradientMagnitudeImageIterator;
  gradientMagnitudeImageIterator = itk::ImageRegionConstIterator<
    GradientMagnitudeOutputImageType>(
    gradientMagnitudeOutputImage,
    gradientMagnitudeOutputImage->GetRequestedRegion() );

  it.GoToBegin();
  eigenVectorImageIterator.GoToBegin();
  eigenValueImageIterator.GoToBegin();
  gradientMagnitudeImageIterator.GoToBegin();

  MatrixType  eigenValueMatrix;
  while( !it.IsAtEnd() )
    {
    // Generate the diagonal matrix with the eigenvalues
    eigenValueMatrix.SetIdentity();

    // Get the eigenvalue
    EigenValueArrayType eigenValue;
    eigenValue = eigenValueImageIterator.Get();

    /* Assumption is that eigenvalue1 > eigenvalue2 > eigenvalue3 */

    // Find the smallest eigenvalue
    double smallest = std::fabs( eigenValue[0] );
    unsigned int smallestEigenValueIndex=0;

    for( unsigned int i=1; i <=2; i++ )
      {
      if( std::fabs( eigenValue[i] ) < smallest )
        {
        smallest = std::fabs( eigenValue[i] );
        smallestEigenValueIndex = i;
        }
      }

    // Find the largest eigenvalue
    double largest = std::fabs( eigenValue[0] );
    unsigned int largestEigenValueIndex=0;

    for( unsigned int i=1; i <=2; i++ )
      {
      if( std::fabs( eigenValue[i] ) > largest )
        {
        largestEigenValueIndex = i;
        }
      }

    unsigned int middleEigenValueIndex=0;

    for( unsigned int i=0; i <=2; i++ )
      {
      if( std::fabs( eigenValue[i] ) != smallest
        && std::fabs( eigenValue[i] ) != largest )
        {
        middleEigenValueIndex = i;
        break;
        }
      }

    //Set the lambda's appropriately.

    //Compute EED lambdas first

    double LambdaEED1;
    double LambdaEED2;
    double LambdaEED3;

    LambdaEED2 = 1.0;
    LambdaEED3 = 1.0;

    double zerovalueTolerance = 1e-15;

    double gradientMagnitude = gradientMagnitudeImageIterator.Get();

    if( gradientMagnitude < zerovalueTolerance )
      {
      LambdaEED1 = 1.0;
      }
    else
      {
      double gradientMagnitudeSquare = gradientMagnitude
        * gradientMagnitude;
      double ratio = ( gradientMagnitudeSquare )
        / ( m_ContrastParameterLambdaEED*m_ContrastParameterLambdaEED );
      double expVal = std::exp( ( -1.0 * m_ThresholdParameterC )
        / ( std::pow( ratio, 4.0 ) ) );
      LambdaEED1 = 1.0 - expVal;
      }

    /* std::cout << "LambdaEED1,LambdaEED2, LambdaEED3\t"
         << LambdaEED1 << "\t" << LambdaEED2 << "\t" << LambdaEED3
         << std::endl; */

    /*Next compute Lambda's for CED */

    double LambdaCED1;
    double LambdaCED2;
    double LambdaCED3;

    LambdaCED1 = m_Alpha;
    LambdaCED2 = m_Alpha;

    double zeroValueTolerance = 1.0e-20;

    if( ( std::fabs( eigenValue[middleEigenValueIndex] ) <
        zeroValueTolerance )  ||
      ( std::fabs( eigenValue[smallestEigenValueIndex] ) <
        zeroValueTolerance ) )
      {
      LambdaCED3 = 1.0;
      }
    else
      {
      double kappa =
       std::pow( ( ( float ) ( eigenValue[middleEigenValueIndex] ) /
                ( m_Alpha + eigenValue[smallestEigenValueIndex] ) ),
               4.0 );

      double contrastParameterLambdaCEDSquare
        = m_ContrastParameterLambdaCED * m_ContrastParameterLambdaCED;

      double expVal = std::exp( ( -1.0 * ( std::log( 2.0 )
        * contrastParameterLambdaCEDSquare )/kappa ) );
      LambdaCED3 = m_Alpha + ( 1.0 - m_Alpha )*expVal;
      }

    /* Compute the lambda's for the continuous switch */
    double Lambda1;
    double Lambda2;
    double Lambda3;

    double xi = ( eigenValue[largestEigenValueIndex]
      /( m_Alpha + eigenValue[middleEigenValueIndex] ) ) -
      ( eigenValue[middleEigenValueIndex]
      /( m_Alpha + eigenValue[smallestEigenValueIndex] ) );

    double numerator = eigenValue[middleEigenValueIndex] *
      ( ( m_ContrastParameterLambdaHybrid *
          m_ContrastParameterLambdaHybrid )
      * ( xi - std::fabs( xi ) ) - 2.0 *
          eigenValue[smallestEigenValueIndex] );


    double denominator = 2.0 * std::pow( m_ContrastParameterLambdaHybrid,
      4.0 );

    double epsilon = std::exp( numerator/denominator );

    Lambda1 = ( 1 - epsilon ) * LambdaCED1 + epsilon*LambdaEED1;
    Lambda2 = ( 1 - epsilon ) * LambdaCED2 + epsilon*LambdaEED2;
    Lambda3 = ( 1 - epsilon ) * LambdaCED3 + epsilon*LambdaEED3;

    eigenValueMatrix( 0, 0 ) = Lambda1;
    eigenValueMatrix( 1, 1 ) = Lambda2;
    eigenValueMatrix( 2, 2 ) = Lambda3;

    //Get the eigenVector matrix
    EigenVectorMatrixType eigenVectorMatrix =
      eigenVectorImageIterator.Get();

    unsigned int vectorLength = 3; // Eigenvector length

    itk::VariableLengthVector<double> firstEigenVector( vectorLength );
    itk::VariableLengthVector<double> secondEigenVector( vectorLength );
    itk::VariableLengthVector<double> thirdEigenVector( vectorLength );

    for( unsigned int i=0; i < vectorLength; i++ ) {
    // Get eigenvectors belonging to eigenvalue order
      firstEigenVector[i] = eigenVectorMatrix[largestEigenValueIndex][i];
      secondEigenVector[i] = eigenVectorMatrix[middleEigenValueIndex][i];
      thirdEigenVector[i] = eigenVectorMatrix[smallestEigenValueIndex][i];

      // Set eigenVectorMatrix in correct order
      eigenVectorMatrix[0][i] = firstEigenVector[i];
      eigenVectorMatrix[1][i] = secondEigenVector[i];
      eigenVectorMatrix[2][i] = thirdEigenVector[i];
      }

    EigenVectorMatrixType  eigenVectorMatrixTranspose;
    eigenVectorMatrixTranspose = eigenVectorMatrix.GetTranspose();

    // Generate the tensor matrix
    EigenVectorMatrixType  productMatrix;
    productMatrix = eigenVectorMatrix * eigenValueMatrix
      * eigenVectorMatrixTranspose;

    //Copy the ITK::Matrix to the tensor...there should be a better way of
    //doing this TODO
    typename DiffusionTensorImageType::PixelType        tensor;

    tensor( 0, 0 ) = productMatrix( 0, 0 );
    tensor( 0, 1 ) = productMatrix( 0, 1 );
    tensor( 0, 2 ) = productMatrix( 0, 2 );

    tensor( 1, 0 ) = productMatrix( 1, 0 );
    tensor( 1, 1 ) = productMatrix( 1, 1 );
    tensor( 1, 2 ) = productMatrix( 1, 2 );

    tensor( 2, 0 ) = productMatrix( 2, 0 );
    tensor( 2, 1 ) = productMatrix( 2, 1 );
    tensor( 2, 2 ) = productMatrix( 2, 2 );
    it.Set( tensor );

    ++it;
    ++eigenValueImageIterator;
    ++eigenVectorImageIterator;
    ++gradientMagnitudeImageIterator;
    }
}

template< class TInputImage, class TOutputImage >
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "CED Contrast parameter "
    << m_ContrastParameterLambdaCED << std::endl;
  os << indent << "EED Contrast parameter "
    << m_ContrastParameterLambdaEED << std::endl;
  os << indent << "Threshold C parameter "
    << m_ThresholdParameterC << std::endl;
  os << indent << "Hybrid Contrast parameter"
    << m_ContrastParameterLambdaHybrid << std::endl;
  os << indent << "Alpha " << m_Alpha << std::endl;
  os << indent << "Sigma " << m_Sigma << std::endl;
  os << indent << "Sigma outer " << m_SigmaOuter << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeAnisotropicHybridDiffusionImageFilter_hxx )
