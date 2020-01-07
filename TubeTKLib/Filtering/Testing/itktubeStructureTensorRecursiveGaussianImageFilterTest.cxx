/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#include "itktubeStructureTensorRecursiveGaussianImageFilter.h"
#include "itktubeSymmetricEigenVectorAnalysisImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>

int itktubeStructureTensorRecursiveGaussianImageFilterTest( int argc,
  char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImage primaryEigenVectorOutputImage";
    std::cerr << " primaryEigenValueOutputImage [Sigma]"<< std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  enum { Dimension = 3 };

  // Define the pixel type
  typedef short PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  InputImageType;

  // Declare the reader
  typedef itk::ImageFileReader< InputImageType > ReaderType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  // Declare the type for the
  typedef itk::tube::StructureTensorRecursiveGaussianImageFilter <
    InputImageType >  StructureTensorFilterType;

  // Create a  Filter
  StructureTensorFilterType::Pointer filter =
    StructureTensorFilterType::New();

  // Connect the input images
  filter->SetInput( reader->GetOutput() );

  // Set the value of sigma if specified in command line
  if( argc > 4 )
    {
    double sigma = std::atof( argv[4] );
    filter->SetSigma( sigma );
    }

  // Execute the filter
  filter->Update();

  // Compute the eigenvectors and eigenvalues of the structure tensor matrix
  typedef  itk::FixedArray< double, Dimension>         EigenValueArrayType;
  typedef  itk::Image< EigenValueArrayType, Dimension> EigenValueImageType;

  typedef  StructureTensorFilterType::OutputImageType
    SymmetricSecondRankTensorImageType;

  typedef itk::SymmetricEigenAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType>
    EigenAnalysisFilterType;

  EigenAnalysisFilterType::Pointer eigenAnalysisFilter =
    EigenAnalysisFilterType::New();
  eigenAnalysisFilter->SetDimension( Dimension );
  eigenAnalysisFilter->OrderEigenValuesBy(
    EigenAnalysisFilterType::EigenValueOrderType::OrderByValue );

  eigenAnalysisFilter->SetInput( filter->GetOutput() );
  eigenAnalysisFilter->Update();

  // Generate eigenvector image
  typedef  itk::Matrix< double, 3, 3>
    EigenVectorMatrixType;
  typedef  itk::Image< EigenVectorMatrixType, Dimension>
    EigenVectorImageType;

  typedef itk::tube::SymmetricEigenVectorAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType,
    EigenVectorImageType> EigenVectorAnalysisFilterType;

  EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter =
    EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( Dimension );
  eigenVectorAnalysisFilter->OrderEigenValuesBy(
    EigenVectorAnalysisFilterType::EigenValueOrderType::OrderByValue );

  eigenVectorAnalysisFilter->SetInput( filter->GetOutput() );
  eigenVectorAnalysisFilter->Update();

  //Generate an image with eigenvector pixel that correspond to the
  //largest eigenvalue
  EigenVectorImageType::ConstPointer eigenVectorImage =
    eigenVectorAnalysisFilter->GetOutput();

  typedef itk::VectorImage< double, 3 >    VectorImageType;
  VectorImageType::Pointer primaryEigenVectorImage = VectorImageType::New();

  unsigned int vectorLength = 3; // Eigenvector length
  primaryEigenVectorImage->SetVectorLength ( vectorLength );

  VectorImageType::RegionType region;
  region.SetSize( eigenVectorImage->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( eigenVectorImage->GetLargestPossibleRegion().GetIndex() );
  primaryEigenVectorImage->SetRegions( region );
  primaryEigenVectorImage->SetOrigin( eigenVectorImage->GetOrigin() );
  primaryEigenVectorImage->SetSpacing( eigenVectorImage->GetSpacing() );
  primaryEigenVectorImage->Allocate();

  //Fill up the buffer with null vector
  itk::VariableLengthVector< double > nullVector( vectorLength );
  for( unsigned int i=0; i < vectorLength; i++ )
    {
    nullVector[i] = 0.0;
    }
  primaryEigenVectorImage->FillBuffer( nullVector );

  //Generate an image containing the largest eigenvalues
  typedef itk::Image< double, 3 >    PrimaryEigenValueImageType;
  PrimaryEigenValueImageType::Pointer primaryEigenValueImage =
    PrimaryEigenValueImageType::New();

  PrimaryEigenValueImageType::RegionType eigenValueImageRegion;
  eigenValueImageRegion.SetSize( eigenVectorImage->
    GetLargestPossibleRegion().GetSize() );
  eigenValueImageRegion.SetIndex( eigenVectorImage->
    GetLargestPossibleRegion().GetIndex() );
  primaryEigenValueImage->SetRegions( eigenValueImageRegion );
  primaryEigenValueImage->SetOrigin( eigenVectorImage->GetOrigin() );
  primaryEigenValueImage->SetSpacing( eigenVectorImage->GetSpacing() );
  primaryEigenValueImage->Allocate();
  primaryEigenValueImage->FillBuffer( 0.0 );

  //Setup the iterators
  //
  //Iterator for the eigenvector matrix image
  itk::ImageRegionConstIterator<EigenVectorImageType>
    eigenVectorImageIterator;
  eigenVectorImageIterator = itk::ImageRegionConstIterator<
    EigenVectorImageType>( eigenVectorImage,
    eigenVectorImage->GetRequestedRegion() );
  eigenVectorImageIterator.GoToBegin();

  //Iterator for the output image with the largest eigenvector
  itk::ImageRegionIterator<VectorImageType> primaryEigenVectorImageIterator;
  primaryEigenVectorImageIterator =
    itk::ImageRegionIterator<VectorImageType>( primaryEigenVectorImage,
    primaryEigenVectorImage->GetRequestedRegion() );
  primaryEigenVectorImageIterator.GoToBegin();

  //Iterator for the output image with the largest eigenvalue
  itk::ImageRegionIterator<PrimaryEigenValueImageType>
    primaryEigenValueImageIterator;
  primaryEigenValueImageIterator = itk::ImageRegionIterator<
    PrimaryEigenValueImageType>( primaryEigenValueImage,
    primaryEigenValueImage->GetRequestedRegion() );
  primaryEigenValueImageIterator.GoToBegin();

  //Iterator for the eigenvalue image
  EigenValueImageType::ConstPointer eigenImage =
    eigenAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenValueImageType>
    eigenValueImageIterator;
  eigenValueImageIterator = itk::ImageRegionConstIterator<
    EigenValueImageType>( eigenImage, eigenImage->GetRequestedRegion() );
  eigenValueImageIterator.GoToBegin();

  //Iterator for the structure tensor
  typedef StructureTensorFilterType::OutputImageType TensorImageType;
  TensorImageType::ConstPointer tensorImage = filter->GetOutput();
  itk::ImageRegionConstIterator<TensorImageType> tensorImageIterator;
  tensorImageIterator = itk::ImageRegionConstIterator<TensorImageType>(
    tensorImage, tensorImage->GetRequestedRegion() );
  tensorImageIterator.GoToBegin();


  double toleranceEigenValues = 1e-4;

  while( !eigenValueImageIterator.IsAtEnd() )
    {
    // Get the eigenvalue
    EigenValueArrayType eigenValue;
    eigenValue = eigenValueImageIterator.Get();

    // Find the largest eigenvalue
    double largest = std::fabs( eigenValue[0] );
    unsigned int largestEigenValueIndex=0;

    for( unsigned int i=1; i <=2; i++ )
      {
      if(  std::fabs( eigenValue[i] > largest ) )
        {
        largest = std::fabs( eigenValue[i] );
        largestEigenValueIndex = i;
        }
      }

    // Write out the largest eigenvalue
    primaryEigenValueImageIterator.Set( eigenValue[largestEigenValueIndex] );


    EigenValueImageType::IndexType pixelIndex;
    pixelIndex = eigenValueImageIterator.GetIndex();

    EigenVectorMatrixType   matrixPixel;
    matrixPixel = eigenVectorImageIterator.Get();

    //Tensor pixelType
    TensorImageType::PixelType tensorPixel;
    tensorPixel = tensorImageIterator.Get();

    /*
    std::cout << "[" << pixelIndex[0] << "," << pixelIndex[1] << ","
    << pixelIndex[2] << "]"
    << "\t" << eigenValue[0] << "\t" << eigenValue[1] << "\t"
    << eigenValue[2] << std::endl;
    std::cout << "[Smallest,Largest]" << "\t" << smallest << "\t"
    << largest << std::endl;
    std::cout <<"eigenvector" << std::endl;
    std::cout <<"\t" <<  matrixPixel << std::endl;
    std::cout <<"Tensor pixel" << std::endl;
    std::cout << "\t" << tensorPixel( 0,0 ) << "\t" << tensorPixel( 0,1 )
    << "\t" << tensorPixel( 0,2 ) << std::endl;
    std::cout << "\t" << tensorPixel( 1,0 ) << "\t" << tensorPixel( 1,1 )
    << "\t" << tensorPixel( 1,2 ) << std::endl;
    std::cout << "\t" << tensorPixel( 2,0 ) << "\t" << tensorPixel( 2,1 )
    << "\t" << tensorPixel( 2,2 ) << std::endl;
    */


    if( std::fabs( largest ) >  toleranceEigenValues  )
      {
      //Assuming eigenvectors are rows
      itk::VariableLengthVector<double> primaryEigenVector( vectorLength );
      for( unsigned int i=0; i < vectorLength; i++ )
      {
      primaryEigenVector[i] = matrixPixel[largestEigenValueIndex][i];
      }

      //std::cout << "\t" << "[" << primaryEigenVector[0] << ","
      //<< primaryEigenVector[1] << "," << primaryEigenVector[2] << "]"
      //<< std::endl;
      primaryEigenVectorImageIterator.Set( primaryEigenVector );
      }

    ++eigenValueImageIterator;
    ++eigenVectorImageIterator;
    ++primaryEigenVectorImageIterator;
    ++primaryEigenValueImageIterator;
    ++tensorImageIterator;

    }

  typedef itk::ImageFileWriter< VectorImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( primaryEigenVectorImage );
  writer->Update();

  typedef itk::ImageFileWriter< PrimaryEigenValueImageType >
    EigenValueImageWriterType;
  EigenValueImageWriterType::Pointer eigenValueWriter =
    EigenValueImageWriterType::New();
  eigenValueWriter->SetFileName( argv[3] );
  eigenValueWriter->SetInput( primaryEigenValueImage );
  eigenValueWriter->Update();


  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
