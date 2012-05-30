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
#ifndef __itkAnisotropicHybridDiffusionImageFilter_txx
#define __itkAnisotropicHybridDiffusionImageFilter_txx

#include "itkAnisotropicHybridDiffusionImageFilter.h"

#include <list>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "itkFixedArray.h"

//#define INTERMEDIATE_OUTPUTS

namespace itk {

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::AnisotropicHybridDiffusionImageFilter()
{
  m_ThresholdParameterC = 3.31488;
  m_ContrastParameterLambdaHybrid = 30.0;
  m_ContrastParameterLambdaCED = 30.0;
  m_ContrastParameterLambdaEED = 20.0;
  m_Sigma = 1.0;
  m_Alpha = 0.001;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
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
  // Instantiate the structure tensor filter
  typename StructureTensorFilterType::Pointer
    StructureTensorFilter  = StructureTensorFilterType::New();

  StructureTensorFilter->SetInput( this->GetOutput() );
  StructureTensorFilter->SetSigma( m_Sigma );
  StructureTensorFilter->Update();

  // Step 1.2: Identify the eigen vectors of the structure tensor
  typedef  Matrix< double, 3, 3>
    EigenVectorMatrixType;
  typedef  Image< EigenVectorMatrixType, 3>
    EigenVectorImageType;
  typedef  typename itk::Image< EigenValueArrayType, 3>
    EigenValueImageType;

  typedef  typename StructureTensorFilterType::OutputImageType
    SymmetricSecondRankTensorImageType;
  typedef  typename itk::SymmetricEigenVectorAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType,
    EigenVectorImageType>
    EigenVectorAnalysisFilterType;

  typename EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter
    = EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( 3 );
  eigenVectorAnalysisFilter->OrderEigenValuesBy(
    EigenVectorAnalysisFilterType::FunctorType::OrderByValue );

  eigenVectorAnalysisFilter->SetInput( StructureTensorFilter->GetOutput() );
  eigenVectorAnalysisFilter->Modified();
  eigenVectorAnalysisFilter->Update();

  //Step 1.3: Compute the eigen values
  typedef itk::SymmetricEigenAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType>
    EigenAnalysisFilterType;

  typename EigenAnalysisFilterType::Pointer eigenAnalysisFilter
    = EigenAnalysisFilterType::New();
  eigenAnalysisFilter->SetDimension( 3 );
  eigenAnalysisFilter->OrderEigenValuesBy(
    EigenAnalysisFilterType::FunctorType::OrderByValue );

  eigenAnalysisFilter->SetInput( StructureTensorFilter->GetOutput() );
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
    eigenVectorImage, eigenVectorImage->GetRequestedRegion());
  eigenVectorImageIterator.GoToBegin();

  //Iterator for the diffusion tensor image
  typedef itk::ImageRegionIterator< DiffusionTensorImageType >
    DiffusionTensorIteratorType;
  DiffusionTensorIteratorType it( this->GetDiffusionTensorImage(),
    this->GetDiffusionTensorImage()->GetLargestPossibleRegion() );

  //Iterator for the eigen value image
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
    // Generate the diagonal matrix with the eigen values
    eigenValueMatrix.SetIdentity();

    // Get the eigen value
    EigenValueArrayType eigenValue;
    eigenValue = eigenValueImageIterator.Get();

    /* Assumption is that eigenvalue1 > eigenvalue2 > eigenvalue3 */

    // Find the smallest eigenvalue
    double smallest = vnl_math_abs( eigenValue[0] );
    unsigned int smallestEigenValueIndex=0;

    for ( unsigned int i=1; i <=2; i++ )
      {
      if ( vnl_math_abs( eigenValue[i] ) < smallest )
        {
        smallest = vnl_math_abs( eigenValue[i] );
        smallestEigenValueIndex = i;
        }
      }

    // Find the largest eigenvalue
    double largest = vnl_math_abs( eigenValue[0] );
    unsigned int largestEigenValueIndex=0;

    for ( unsigned int i=1; i <=2; i++ )
      {
      if (  vnl_math_abs( eigenValue[i] > largest ) )
        {
        largestEigenValueIndex = i;
        }
      }

    unsigned int middleEigenValueIndex=0;

    for ( unsigned int i=0; i <=2; i++ )
      {
      if ( eigenValue[i] != smallest && eigenValue[i] != largest )
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

    if ( gradientMagnitude < zerovalueTolerance )
      {
      LambdaEED1 = 1.0;
      }
    else
      {
      double gradientMagnitudeSquare = gradientMagnitude
        * gradientMagnitude;
      double ratio = (gradientMagnitudeSquare) /
               (m_ContrastParameterLambdaEED*m_ContrastParameterLambdaEED);
      double expVal = exp( (-1.0 * m_ThresholdParameterC)/(vcl_pow( ratio,
        4.0 )));
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

    if((fabs(eigenValue[middleEigenValueIndex]) < zeroValueTolerance)  ||
       (fabs(eigenValue[smallestEigenValueIndex]) < zeroValueTolerance) )
      {
      LambdaCED3 = 1.0;
      }
    else
      {
      double kappa =
       vcl_pow( ((float) (eigenValue[middleEigenValueIndex]) /
                ( m_Alpha + eigenValue[smallestEigenValueIndex])),
               4.0);

      double contrastParameterLambdaCEDSquare
        = m_ContrastParameterLambdaCED * m_ContrastParameterLambdaCED;

      double expVal = exp((-1.0 * (vcl_log( 2.0)
        * contrastParameterLambdaCEDSquare )/kappa ));
      LambdaCED3 = m_Alpha + (1.0 - m_Alpha)*expVal;
      }

    /* Compute the lambda's for the continous switch */
    double Lambda1;
    double Lambda2;
    double Lambda3;

    double xi = (eigenValue[largestEigenValueIndex]
      /(m_Alpha + eigenValue[middleEigenValueIndex])) -
      (eigenValue[middleEigenValueIndex]
      /(m_Alpha + eigenValue[smallestEigenValueIndex]));

    double numerator = eigenValue[middleEigenValueIndex] *
      ( (m_ContrastParameterLambdaHybrid * m_ContrastParameterLambdaHybrid)
      * (xi - fabs(xi)) - 2.0 * eigenValue[smallestEigenValueIndex] );


    double denominator = 2.0 * vcl_pow( m_ContrastParameterLambdaHybrid,
      4.0 );

    double epsilon = exp(numerator/denominator);

    Lambda1 = (1 - epsilon ) * LambdaCED1 + epsilon*LambdaEED1;
    Lambda2 = (1 - epsilon ) * LambdaCED2 + epsilon*LambdaEED2;
    Lambda3 = (1 - epsilon ) * LambdaCED3 + epsilon*LambdaEED3;

    eigenValueMatrix(0,0) = Lambda1;
    eigenValueMatrix(1,1) = Lambda2;
    eigenValueMatrix(2,2) = Lambda3;

    //Get the eigenVector matrix
    EigenVectorMatrixType eigenVectorMatrix;

    unsigned int vectorLength = 3; // Eigenvector length

    itk::VariableLengthVector<double> firstEigenVector( vectorLength );
    itk::VariableLengthVector<double> secondEigenVector( vectorLength );
    itk::VariableLengthVector<double> thirdEigenVector( vectorLength );

    for ( unsigned int i=0; i < vectorLength; i++ ) {
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
    ++gradientMagnitudeImageIterator;
    }
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::SetSigma( double sigma)
{
  m_Sigma = sigma;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::SetThresholdParameterC( double threshold)
{
  m_ThresholdParameterC = threshold;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::SetContrastParameterLambdaEED( double contrast)
{
  m_ContrastParameterLambdaEED = contrast;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::SetContrastParameterLambdaCED( double contrast)
{
  m_ContrastParameterLambdaCED = contrast;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::SetContrastParameterLambdaHybrid( double contrast)
{
  m_ContrastParameterLambdaHybrid = contrast;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::SetAlpha( double alpha)
{
  m_Alpha = alpha;
}

template <class TInputImage, class TOutputImage>
void
AnisotropicHybridDiffusionImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "CED Contrast parameter "
    << m_ContrastParameterLambdaCED << std::endl;
  os << indent << "EED Contrast parameter "
    << m_ContrastParameterLambdaEED << std::endl;
  os << indent << "Hybrid Contrast parameter"
    << m_ContrastParameterLambdaHybrid << std::endl;
  os << indent << "Alpha " << m_Alpha << std::endl;
  os << indent << "Sigma " << m_Sigma << std::endl;
}

} // end namespace itk

#endif
