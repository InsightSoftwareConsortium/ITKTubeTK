#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkVesselEnhancingDiffusion2DImageFilter.h"

#include "HessianVesselness2DCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef float                                                   PixelType;
  typedef itk::Image< PixelType,  2 >                             ImageType;
  typedef itk::ImageFileReader< ImageType >                       ReaderType;
  typedef itk::VesselEnhancingDiffusion2DImageFilter< PixelType > FilterType;
  
  ReaderType::Pointer reader = ReaderType::New();

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
  
  FilterType::Pointer filter = NULL;
  filter = FilterType::New();
  
  // set input image
  filter->SetInput( reader->GetOutput() );
  filter->SetDefaultPars();

  // set paramters
  filter->SetTimeStep( timeStep );
  filter->SetIterations( numIterations );
  filter->SetRecalculateVesselness( recalculateVesselness );
  
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
  ImageWriterType::Pointer writer = ImageWriterType::New();

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

