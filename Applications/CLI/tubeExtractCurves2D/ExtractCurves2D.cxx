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

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// needed ITK classes
#include "itkImage.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkBinaryMagnitudeImageFilter.h"
#include "itkEigenAnalysis2DImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageToParametricSpaceFilter.h"
#include "itkMesh.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkJoinImageFilter.h"
#include "itkFrustumSpatialFunction.h"
#include "itkInteriorExteriorMeshFilter.h"
#include "itkParametricSpaceToImageSpaceMeshFilter.h"
#include "itkRescaleIntensityImageFilter.h"

// The following three should be used in every CLI application
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "ExtractCurves2DCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  // of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "ExtractCurves2D",
                                                 CLPProcessInformation );
  progressReporter.Start();

  // Input Types
  typedef   double                                   PixelType;
  typedef   itk::Image< PixelType, 2 >               ImageType;

  // Output Types
  typedef   unsigned char                            OutputPixelType;
  typedef   itk::Image< OutputPixelType, 2 >         OutputImageType;

  // Vector Image Types
  typedef   itk::Vector< PixelType, 2 >              VectorType;
  typedef   itk::CovariantVector< PixelType, 2 >     CovariantVectorType;
  typedef   itk::Image< VectorType, 2 >              VectorImageType;
  typedef   itk::Image< CovariantVectorType, 2 >     CovariantVectorImageType;

  // Mesh Types
  typedef   itk::Point<float,2>                      MeshPointDataType;
  typedef   itk::Mesh< MeshPointDataType, 3 >        MeshType;

  typedef   itk::ImageFileReader< ImageType >        ReaderType;
  typedef   itk::ImageFileWriter< OutputImageType >  WriterType;

  typedef   itk::Mesh< MeshType::PointType, 2 >      ImageSpaceMeshType;

  typedef   itk::RecursiveGaussianImageFilter< ImageType, ImageType >
                                                     GaussianFilterType;

  typedef   itk::AddImageFilter< ImageType, ImageType, ImageType >
                                                     AddFilterType;

  typedef   itk::BinaryMagnitudeImageFilter< ImageType, ImageType, ImageType >
                                                     ModulusFilterType;

  typedef   itk::EigenAnalysis2DImageFilter< ImageType, ImageType, VectorImageType >
                                                     EigenFilterType;

  typedef   itk::GradientRecursiveGaussianImageFilter< ImageType, CovariantVectorImageType >
                                                     GradientFilterType;

  typedef   itk::MultiplyImageFilter< VectorImageType, VectorImageType, ImageType >
                                                     ScalarProductFilterType;

  typedef   itk::ImageToParametricSpaceFilter< ImageType, MeshType >
                                                     ParametricSpaceFilterType;

  typedef   itk::JoinImageFilter< ImageType, ImageType >
                                                     JoinFilterType;

  typedef   itk::RescaleIntensityImageFilter< ImageType, ImageType >
                                                     RescaleIntensityFilterType;

  typedef   itk::FrustumSpatialFunction< MeshType::PointDimension,
                                         MeshType::PointType >
                                                     FrustumSpatialFunctionType;

  typedef   FrustumSpatialFunctionType               SpatialFunctionType;

  typedef   itk::InteriorExteriorMeshFilter< MeshType, MeshType, SpatialFunctionType >
                                                     SpatialFunctionFilterType;

  typedef   itk::ParametricSpaceToImageSpaceMeshFilter< MeshType, ImageSpaceMeshType >
                                                     InverseParametricFilterType;

  typedef typename GaussianFilterType::RealType      RealType;

  typedef typename ImageSpaceMeshType::PointsContainer
                                                     OutputPointsContainer;
  typedef typename OutputPointsContainer::ConstPointer
                                                     OutputPointsContainerConstPointer;
  typedef typename OutputPointsContainer::ConstIterator
                                                     OutputPointsConstIterator;
  typedef itk::MinimumMaximumImageCalculator< ImageType >
                                                     CalculatorType;
  typedef itk::ImageRegionConstIterator< ImageType > InputIterator;
  typedef itk::ImageRegionIterator< OutputImageType >
                                                     OutputIterator;

  typename ReaderType::Pointer                       reader;
  typename WriterType::Pointer                       writer;
  typename GaussianFilterType::Pointer               hx;
  typename GaussianFilterType::Pointer               hy;
  typename GaussianFilterType::Pointer               hxy;
  typename GaussianFilterType::Pointer               h1xy;
  typename GaussianFilterType::Pointer               h1x;
  typename GaussianFilterType::Pointer               h1y;
  typename GaussianFilterType::Pointer               h2x;
  typename GaussianFilterType::Pointer               h2y;
  typename AddFilterType::Pointer                    add;
  typename ModulusFilterType::Pointer                modulus;
  typename EigenFilterType::Pointer                  eigen;
  typename GradientFilterType::Pointer               gradient;
  typename ScalarProductFilterType::Pointer          scalarProduct;
  typename JoinFilterType::Pointer                   join;
  typename RescaleIntensityFilterType::Pointer       rescaleIntensitySmoothed;
  typename RescaleIntensityFilterType::Pointer       rescaleIntensityMaxEigen;
  typename RescaleIntensityFilterType::Pointer       rescaleIntensityMedialness;
  typename ParametricSpaceFilterType::Pointer        parametricSpace;
  typename SpatialFunctionFilterType::Pointer        spatialFunctionFilter;
  typename InverseParametricFilterType::Pointer      inverseParametricFilter;

  // Get ready for some PROGRESS!!! Wooo!!
  double progress = 0;

  timeCollector.Start("Load Data");
  reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );

  // Load the input image with exception handling
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input image!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load Data");

  // Setup the string of filters
  rescaleIntensitySmoothed   = RescaleIntensityFilterType::New();
  rescaleIntensityMedialness = RescaleIntensityFilterType::New();
  rescaleIntensityMaxEigen   = RescaleIntensityFilterType::New();

  hx = GaussianFilterType::New();
  hy = GaussianFilterType::New();
  hx->SetDirection( 0 );
  hy->SetDirection( 1 );

  hxy  = GaussianFilterType::New();
  hxy->SetDirection( 1 );

  h1x  = GaussianFilterType::New();
  h1y  = GaussianFilterType::New();
  h1x->SetDirection( 0 );
  h1y->SetDirection( 1 );
  h1x->SetOrder( GaussianFilterType::FirstOrder );
  h1y->SetOrder( GaussianFilterType::FirstOrder );

  h1xy = GaussianFilterType::New();
  h1xy->SetDirection( 1 );
  h1xy->SetOrder( GaussianFilterType::FirstOrder );

  h2x  = GaussianFilterType::New();
  h2y  = GaussianFilterType::New();

  h2x->SetDirection( 0 );
  h2y->SetDirection( 1 );

  h2x->SetOrder( GaussianFilterType::SecondOrder );
  h2y->SetOrder( GaussianFilterType::SecondOrder );

  hx->SetInputImage( reader->GetOutput() );
  hy->SetInputImage( reader->GetOutput() );

  hxy->SetInputImage( hx->GetOutput() );

  h1x->SetInputImage( hy->GetOutput() );
  h1y->SetInputImage( hx->GetOutput() );

  h2x->SetInputImage( hy->GetOutput() );
  h2y->SetInputImage( hx->GetOutput() );

  h1xy->SetInputImage( h1x->GetOutput() );

  add = AddFilterType::New();

  add->SetInput1( h2x->GetOutput() );
  add->SetInput2( h2y->GetOutput() );

  modulus = ModulusFilterType::New();

  modulus->SetInput1( h1x->GetOutput() );
  modulus->SetInput2( h1y->GetOutput() );

  gradient = GradientFilterType::New();

  gradient->SetInput( reader->GetOutput() );

  eigen = EigenFilterType::New();

  eigen->SetInput1( h2x->GetOutput() );
  eigen->SetInput2( h1xy->GetOutput() );
  eigen->SetInput3( h2y->GetOutput() );

  join = JoinFilterType::New();

  join->SetInput1( h1x->GetOutput() );
  join->SetInput2( h1y->GetOutput() );

  scalarProduct = ScalarProductFilterType::New();

  scalarProduct->SetInput1( join->GetOutput() );
  scalarProduct->SetInput2( eigen->GetMaxEigenVector() );

  // Normalize the parametric space
  rescaleIntensitySmoothed->SetInput( hxy->GetOutput() );
  rescaleIntensitySmoothed->SetOutputMinimum( -1.0 );
  rescaleIntensitySmoothed->SetOutputMaximum(  1.0 );

  rescaleIntensityMedialness->SetInput( scalarProduct->GetOutput() );
  rescaleIntensityMedialness->SetOutputMinimum( -1.0 );
  rescaleIntensityMedialness->SetOutputMaximum(  1.0 );

  rescaleIntensityMaxEigen->SetInput( eigen->GetMaxEigenValue() );
  rescaleIntensityMaxEigen->SetOutputMinimum( 0.0 );
  rescaleIntensityMaxEigen->SetOutputMaximum( 1.0 );

  parametricSpace = ParametricSpaceFilterType::New();

  parametricSpace->SetInput( 0, rescaleIntensityMaxEigen->GetOutput() );
  parametricSpace->SetInput( 1, rescaleIntensityMedialness->GetOutput() );
  parametricSpace->SetInput( 2, rescaleIntensitySmoothed->GetOutput() );

  spatialFunctionFilter = SpatialFunctionFilterType::New();

  spatialFunctionFilter->SetInput( parametricSpace->GetOutput() );

  inverseParametricFilter = InverseParametricFilterType::New();

  inverseParametricFilter->SetInput( spatialFunctionFilter->GetOutput() );

  SpatialFunctionType::Pointer spatialFunction =
    spatialFunctionFilter->GetSpatialFunction();

  itk::Point<float,3> apex;
  apex[0] = 0.4;
  apex[1] = 0;
  apex[2] = 1;

  spatialFunction->SetAngleZ( 20.0f );
  spatialFunction->SetApertureAngleX( 12.0f );
  spatialFunction->SetApertureAngleY(  3.0f );
  spatialFunction->SetTopPlane( 0.1f );
  spatialFunction->SetBottomPlane( 2.0f );
  spatialFunction->SetApex( apex );

  // Set the sigma value to all the different filters
  hx->SetSigma( sigma );
  hy->SetSigma( sigma );
  hxy->SetSigma( sigma );
  h1x->SetSigma( sigma );
  h1y->SetSigma( sigma );
  h1xy->SetSigma( sigma );
  h2x->SetSigma( sigma );
  h2y->SetSigma( sigma );
  gradient->SetSigma( sigma );

  // Update the pipeline
  hxy->Update();
  h1xy->Update();
  add->Update();
  modulus->Update();
  inverseParametricFilter->Update();

  // Get the values into the output
  ImageType::ConstPointer inputImage = reader->GetOutput();

  ImageType::RegionType region = inputImage->GetLargestPossibleRegion();
  OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( region );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  ImageSpaceMeshType::ConstPointer outputMesh =
    inverseParametricFilter->GetOutput();

  OutputPointsContainerConstPointer outputPoints = outputMesh->GetPoints();

  OutputPointsConstIterator pointItr = outputPoints->Begin();
  OutputPointsConstIterator pointEnd = outputPoints->End();

  OutputImageType::IndexType    outputIndex;
  ImageSpaceMeshType::PointType outputPoint;

  while( pointItr != pointEnd )
    {
    outputPoint = pointItr.Value();
    outputImage->TransformPhysicalPointToIndex( outputPoint, outputIndex );
    outputImage->SetPixel( outputIndex, 255 );
    ++pointItr;
    }

  // Prepare to write the output image of scores to disk
  timeCollector.Start("Save data");
  writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( outputImage );

  // Write the output image to disk with exception handling
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");

  // Tell the progress reporter that we're done
  progress = 1.0;
  progressReporter.Report( progress );

  // Report the time each process step took
  timeCollector.Report();

  // Return with a good value
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
