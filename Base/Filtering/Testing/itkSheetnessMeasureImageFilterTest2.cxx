/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSheetnessMeasureImageFilter.h"
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkSymmetricEigenVectorAnalysisImageFilter.h>
#include <itkMatrix.h>
#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>


int itkSheetnessMeasureImageFilterTest2(int argc ,char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImage sheetnessImage primaryEigenVectorOutputImage [sheetnessThreshold]"<< std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const unsigned int Dimension = 3;

  // Define the pixel type
  typedef short PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>           InputImageType;

  // Declare the reader
  typedef itk::ImageFileReader< InputImageType > ReaderType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  // Declare the type for the Hessian filter
  typedef itk::HessianRecursiveGaussianImageFilter<
                                            InputImageType >  HessianFilterType;

  typedef HessianFilterType::OutputImageType HessianImageType;

  // Declare the type for the sheetness measure filter
  typedef itk::SheetnessMeasureImageFilter< float >  SheetnessFilterType;

  // Create a Hessian Filter
  HessianFilterType::Pointer filterHessian = HessianFilterType::New();

  // Create a sheetness Filter
  SheetnessFilterType::Pointer filterSheetness = SheetnessFilterType::New();

  // Connect the input images
  filterHessian->SetInput( reader->GetOutput() );
  filterSheetness->SetInput( filterHessian->GetOutput() );

  // Select the value of Sigma
  filterHessian->SetSigma( 0.5 );

  // Execute the filter
  std::cout << "Generate sheetness measure" << std::endl;
  filterSheetness->Update();

  //Write out the sheetness image
  //Define output type
  typedef SheetnessFilterType::OutputImageType SheetnessImageType;

  std::cout << "Write out the sheetness image" << std::endl;
  typedef itk::ImageFileWriter<SheetnessImageType>     SheetnessImageWriterType;
  SheetnessImageWriterType::Pointer writer= SheetnessImageWriterType::New();
  writer->SetFileName(argv[1]);
  writer->SetInput(filterSheetness->GetOutput());
  writer->Update();

  //Compute the eigen values

  typedef  SheetnessFilterType::InputImageType SymmetricSecondRankTensorImageType;
  typedef  itk::FixedArray< double, Dimension>    EigenValueArrayType;
  typedef  itk::Image< EigenValueArrayType, Dimension> EigenValueImageType;

  typedef itk::
    SymmetricEigenAnalysisImageFilter<SymmetricSecondRankTensorImageType, EigenValueImageType> EigenAnalysisFilterType;

  EigenAnalysisFilterType::Pointer eigenAnalysisFilter = EigenAnalysisFilterType::New();
  eigenAnalysisFilter->SetDimension( Dimension );
  eigenAnalysisFilter->OrderEigenValuesBy(
      EigenAnalysisFilterType::FunctorType::OrderByValue );

  eigenAnalysisFilter->SetInput( filterHessian->GetOutput() );
  eigenAnalysisFilter->Update();


  //Generate and write out the primary eigen vector image
  typedef  itk::Matrix< double, 3, 3>                         EigenVectorMatrixType;
  typedef  itk::Image< EigenVectorMatrixType, Dimension>      EigenVectorImageType;

    typedef  itk::SymmetricEigenVectorAnalysisImageFilter<SymmetricSecondRankTensorImageType, EigenValueImageType, EigenVectorImageType> EigenVectorAnalysisFilterType;

  EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter = EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( Dimension );
  eigenVectorAnalysisFilter->OrderEigenValuesBy(
      EigenVectorAnalysisFilterType::FunctorType::OrderByValue );

  eigenVectorAnalysisFilter->SetInput( filterHessian->GetOutput() );
  eigenVectorAnalysisFilter->Update();

  //Generate an image with eigen vector pixel that correspond to the largest eigen value
  EigenVectorImageType::ConstPointer eigenVectorImage =
                    eigenVectorAnalysisFilter->GetOutput();

  typedef itk::VectorImage< double, 3 >    VectorImageType;
  VectorImageType::Pointer primaryEigenVectorImage = VectorImageType::New();

  unsigned int vectorLength = 3; // Eigenvector length
  primaryEigenVectorImage->SetVectorLength ( vectorLength );

  VectorImageType::RegionType region;
  region.SetSize(eigenVectorImage->GetLargestPossibleRegion().GetSize());
  region.SetIndex(eigenVectorImage->GetLargestPossibleRegion().GetIndex());
  primaryEigenVectorImage->SetRegions( region );
  primaryEigenVectorImage->SetOrigin(eigenVectorImage->GetOrigin());
  primaryEigenVectorImage->SetSpacing(eigenVectorImage->GetSpacing());
  primaryEigenVectorImage->Allocate();

  //Fill up the buffer with null vector
  itk::VariableLengthVector< double > nullVector( vectorLength );
  for ( unsigned int i=0; i < vectorLength; i++ )
    {
    nullVector[i] = 0.0;
    }
  primaryEigenVectorImage->FillBuffer( nullVector );

  //Generate an image containing the largest eigen values
  typedef itk::Image< double, 3 >    PrimaryEigenValueImageType;
  PrimaryEigenValueImageType::Pointer primaryEigenValueImage = PrimaryEigenValueImageType::New();

  PrimaryEigenValueImageType::RegionType eigenValueImageRegion;
  eigenValueImageRegion.SetSize(eigenVectorImage->GetLargestPossibleRegion().GetSize());
  eigenValueImageRegion.SetIndex(eigenVectorImage->GetLargestPossibleRegion().GetIndex());
  primaryEigenValueImage->SetRegions( eigenValueImageRegion );
  primaryEigenValueImage->SetOrigin(eigenVectorImage->GetOrigin());
  primaryEigenValueImage->SetSpacing(eigenVectorImage->GetSpacing());
  primaryEigenValueImage->Allocate();
  primaryEigenValueImage->FillBuffer( 0.0 );


  //Setup the iterators
  //
  //Iterator for the eigenvector matrix image
  itk::ImageRegionConstIterator<EigenVectorImageType> eigenVectorImageIterator;
  eigenVectorImageIterator = itk::ImageRegionConstIterator<EigenVectorImageType>(
      eigenVectorImage, eigenVectorImage->GetRequestedRegion());
  eigenVectorImageIterator.GoToBegin();

  //Iterator for the output image with the largest eigenvector
  itk::ImageRegionIterator<VectorImageType> primaryEigenVectorImageIterator;
  primaryEigenVectorImageIterator = itk::ImageRegionIterator<VectorImageType>(
      primaryEigenVectorImage, primaryEigenVectorImage->GetRequestedRegion());
  primaryEigenVectorImageIterator.GoToBegin();

  //Iterator for the output image with the largest eigenvalue
  itk::ImageRegionIterator<PrimaryEigenValueImageType> primaryEigenValueImageIterator;
  primaryEigenValueImageIterator = itk::ImageRegionIterator<PrimaryEigenValueImageType>(
      primaryEigenValueImage, primaryEigenValueImage->GetRequestedRegion());
  primaryEigenValueImageIterator.GoToBegin();


  //Iterator for the eigen value image
  EigenValueImageType::ConstPointer eigenImage = eigenAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenValueImageType> eigenValueImageIterator;
  eigenValueImageIterator = itk::ImageRegionConstIterator<EigenValueImageType>(
      eigenImage, eigenImage->GetRequestedRegion());
  eigenValueImageIterator.GoToBegin();

  double toleranceEigenValues = 1e-4;

  while (!eigenValueImageIterator.IsAtEnd())
    {
    // Get the eigen value
    EigenValueArrayType eigenValue;
    eigenValue = eigenValueImageIterator.Get();

    // Find the smallest eigenvalue
    double smallest = vnl_math_abs( eigenValue[0] );
    double Lambda1 = eigenValue[0];
    unsigned int smallestEigenValueIndex=0;

    for ( unsigned int i=1; i <=2; i++ )
      {
      if ( vnl_math_abs( eigenValue[i] ) < smallest )
        {
        Lambda1 = eigenValue[i];
        smallest = vnl_math_abs( eigenValue[i] );
        smallestEigenValueIndex = i;
        }
      }

    // Find the largest eigenvalue
    double largest = vnl_math_abs( eigenValue[0] );
    double Lambda3 = eigenValue[0];
    unsigned int largestEigenValueIndex=0;

    for ( unsigned int i=1; i <=2; i++ )
      {
      if (  vnl_math_abs( eigenValue[i] > largest ) )
        {
        Lambda3 = eigenValue[i];
        largest = vnl_math_abs( eigenValue[i] );
        largestEigenValueIndex = i;
        }
      }

    // find Lambda2 so that |Lambda1| < |Lambda2| < |Lambda3|
    double Lambda2 = eigenValue[0];
    unsigned int middleEigenValueIndex=0;

    for ( unsigned int i=0; i <=2; i++ )
      {
      if ( eigenValue[i] != Lambda1 && eigenValue[i] != Lambda3 )
        {
        Lambda2 = eigenValue[i];
        middleEigenValueIndex = i;
        break;
        }
      }

    // Write out the largest eigen value
    primaryEigenValueImageIterator.Set( eigenValue[largestEigenValueIndex]);


    EigenValueImageType::IndexType pixelIndex;
    pixelIndex = eigenValueImageIterator.GetIndex();

    EigenVectorMatrixType   matrixPixel;
    matrixPixel = eigenVectorImageIterator.Get();

    if( fabs(largest) >  toleranceEigenValues  )
      {
      //Assuming eigenvectors are rows
      itk::VariableLengthVector<double> primaryEigenVector( vectorLength );
      for ( unsigned int i=0; i < vectorLength; i++ )
      {
      primaryEigenVector[i] = matrixPixel[largestEigenValueIndex][i];
      }

      primaryEigenVectorImageIterator.Set( primaryEigenVector );
      }

    ++eigenValueImageIterator;
    ++eigenVectorImageIterator;
    ++primaryEigenVectorImageIterator;
    ++primaryEigenValueImageIterator;
    }

  typedef itk::ImageFileWriter< VectorImageType > EigenVectorWriterType;
  EigenVectorWriterType::Pointer eigenVectorWriter = EigenVectorWriterType::New();
  eigenVectorWriter->SetFileName( argv[2] );
  eigenVectorWriter->SetInput( primaryEigenVectorImage );
  eigenVectorWriter->Update();

  return EXIT_SUCCESS;

}
