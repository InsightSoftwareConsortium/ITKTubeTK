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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include <itkImage.h>
#include <itkStructureTensorRecursiveGaussianImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkSymmetricEigenVectorAnalysisImageFilter.h>
#include <itkMatrix.h>
#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>


int itkStructureTensorRecursiveGaussianImageFilterTest(int argc, char* argv []  ) 
{
  if( argc < 3 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImage outputImage"<< std::endl;
    return EXIT_FAILURE;
    }
 
  // Define the dimension of the images
  const unsigned int Dimension = 3;

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
  typedef itk::StructureTensorRecursiveGaussianImageFilter < 
                                            InputImageType >  StructureTensorFilterType;
            
  typedef StructureTensorFilterType::OutputImageType myTensorImageType;

  // Create a  Filter                                
  StructureTensorFilterType::Pointer filter = StructureTensorFilterType::New();

  // Connect the input images
  filter->SetInput( reader->GetOutput() ); 

  // Select the value of Sigma
  filter->SetSigma( 2.5 ); 

  // Execute the filter
  filter->Update();

  // Compute the eigenvectors and eigenvalues of the structure tensor matrix 
  typedef  itk::FixedArray< double, Dimension>    EigenValueArrayType;
  typedef  itk::Image< EigenValueArrayType, Dimension> EigenValueImageType;

  typedef  StructureTensorFilterType::OutputImageType SymmetricSecondRankTensorImageType;

  typedef itk::
    SymmetricEigenAnalysisImageFilter<SymmetricSecondRankTensorImageType, EigenValueImageType> EigenAnalysisFilterType;

  EigenAnalysisFilterType::Pointer eigenAnalysisFilter = EigenAnalysisFilterType::New();
  eigenAnalysisFilter->SetDimension( Dimension );
  eigenAnalysisFilter->OrderEigenValuesBy( 
      EigenAnalysisFilterType::FunctorType::OrderByValue );
  
  eigenAnalysisFilter->SetInput( filter->GetOutput() );
  eigenAnalysisFilter->Update();

  // Generate eigen vector image
  typedef  itk::Matrix< double, 3, 3>                           EigenVectorMatrixType;
  typedef  itk::Image< EigenVectorMatrixType, Dimension>      EigenVectorImageType;

  typedef itk::
    SymmetricEigenVectorAnalysisImageFilter<SymmetricSecondRankTensorImageType, EigenValueImageType, EigenVectorImageType> EigenVectorAnalysisFilterType;

  EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter = EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( Dimension );
  eigenVectorAnalysisFilter->OrderEigenValuesBy( 
      EigenVectorAnalysisFilterType::FunctorType::OrderByValue );
  
  eigenVectorAnalysisFilter->SetInput( filter->GetOutput() );
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
  primaryEigenVectorImage->SetRegions( region );
  primaryEigenVectorImage->Allocate();

  //Fill up the buffer with null vector
  itk::VariableLengthVector< double > nullVector( vectorLength );
  for ( unsigned int i=0; i < vectorLength; i++ )
    {
    nullVector[i] = 0.0;
    }
  primaryEigenVectorImage->FillBuffer( nullVector );


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


    EigenValueImageType::IndexType pixelIndex;
    pixelIndex = eigenValueImageIterator.GetIndex();

    std::cout << "[" << pixelIndex[0] << "," << pixelIndex[1] << "," << pixelIndex[2] << "]" << "\t" << smallest << "\t" << largest << std::endl;

    if( fabs(largest) >  toleranceEigenValues  )
      {
      EigenVectorMatrixType   matrixPixel;
      matrixPixel = eigenVectorImageIterator.Get();
      //Assuming eigenvectors are columns
      itk::VariableLengthVector<double> primaryEigenVector( vectorLength );
      for ( unsigned int i=0; i < vectorLength; i++ )
      {
      primaryEigenVector[i] = matrixPixel[largestEigenValueIndex][i];
      }

      std::cout << "\t" << "[" << primaryEigenVector[0] << "," << primaryEigenVector[1] << "," << primaryEigenVector[2] << "]" << std::endl;
      primaryEigenVectorImageIterator.Set( primaryEigenVector );
      }

    ++eigenValueImageIterator;
    ++eigenVectorImageIterator;
    ++primaryEigenVectorImageIterator;
    }

  typedef itk::ImageFileWriter< VectorImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( primaryEigenVectorImage );
  writer->Update();

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;

}




