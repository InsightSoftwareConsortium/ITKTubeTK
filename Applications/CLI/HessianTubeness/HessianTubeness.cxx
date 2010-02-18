
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTubeEnhancingDiffusion2DImageFilter.h"

// forward declaration
// - Required for tubeCLIHelperFunctions to work
template< class pixelT, unsigned int dimensionT > 
int DoIt( int argc, char **argv );

// CLP file must be included before tubeCLIHelperFunctions
#include "HessianTubenessCLP.h"
#include "tubeCLIHelperFunctions.h"

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char **argv )
  {
  PARSE_ARGS;

  typedef pixelT                                           PixelType;
  typedef itk::Image< PixelType,  dimensionT  >            ImageType;
  typedef itk::ImageFileReader< ImageType >                ReaderType;
  typedef itk::TubeEnhancingDiffusion2DImageFilter< PixelType,
                                                    dimensionT  >
                                                           FilterType;
  
  typename ReaderType::Pointer reader = ReaderType::New();

  //read input image  
  reader->SetFileName( inputVolume.c_str() );
  
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  
  typename FilterType::Pointer filter = NULL;
  filter = FilterType::New();
  
  // set input image
  filter->SetInput( reader->GetOutput() );
  filter->SetDefaultPars();

  // set paramters
  filter->SetTimeStep( timeStep );
  filter->SetIterations( numIterations );
  filter->SetRecalculateTubeness( recalculateTubeness );
  
  filter->SetBeta( beta );
  filter->SetGamma( gamma );
 
  filter->SetEpsilon( epsilon );
  filter->SetOmega( omega );
  filter->SetSensitivity( sensitivity );
 
  // Compute scales and then set them
  std::vector< float > scales( numSigmaSteps );
  double deltaSigma = maxSigma - minSigma;
  for( int i = 0; i <= numSigmaSteps; i++ )
    {
    scales[i] = minSigma + i * ( deltaSigma / numSigmaSteps );
    }
  filter->SetScales( scales );

  // compute vesselness image
  filter->Update();

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput ( filter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }

int main( int argc, char * argv[] )
  {   
  PARSE_ARGS;  

  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
  }

