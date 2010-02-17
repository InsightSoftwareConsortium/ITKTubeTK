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
#include "HessianTubenessCLP.h"
#include "tubetkParseImageInformation.h"

template< class T, unsigned int dimension> int DoIt( int argc, char * argv[])
{
  PARSE_ARGS;

  typedef T                                                      PixelType;
  typedef itk::Image< PixelType,  dimension >                    ImageType;
  typedef itk::ImageFileReader< ImageType >                      ReaderType;
  typedef itk::TubeEnhancingDiffusion2DImageFilter< PixelType, dimension >  FilterType;
  
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

  itk::ImageIOBase::IOComponentType componentType;
  unsigned int dimension;

  try    
    {    
    GetImageInformation(inputVolume, componentType, dimension);   
    switch (componentType)      
    {      
    case itk::ImageIOBase::UCHAR:        
      if( dimension == 2 )
        {
        return DoIt<unsigned char, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<unsigned char, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break;      
    case itk::ImageIOBase::CHAR:        
      if( dimension == 2 )
        {
        return DoIt<char, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<char, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
    case itk::ImageIOBase::USHORT:        
      if( dimension == 2 )
        {
        return DoIt<unsigned short, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<unsigned short, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break;
    case itk::ImageIOBase::SHORT:        
      if( dimension == 2 )
        {
        return DoIt<short, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<short, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
     case itk::ImageIOBase::UINT:        
      if( dimension == 2 )
        {
        return DoIt<unsigned int, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<unsigned int, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
     case itk::ImageIOBase::INT:        
      if( dimension == 2 )
        {
        return DoIt<int, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<int, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
    case itk::ImageIOBase::ULONG:        
      if( dimension == 2 )
        {
        return DoIt<unsigned long, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<unsigned long, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
    case itk::ImageIOBase::LONG:        
      if( dimension == 2 )
        {
        return DoIt<long, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<long, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
    case itk::ImageIOBase::FLOAT:        
      if( dimension == 2 )
        {
        return DoIt<float, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<float, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
    case itk::ImageIOBase::DOUBLE:        
      if( dimension == 2 )
        {
        return DoIt<double, 2>( argc, argv);
        }
      else if( dimension == 3 ) 
        {
        return DoIt<double, 3>( argc, argv);
        }
      else  
        {
        std::cerr << "Unsupported image dimension" << std::endl;
        return EXIT_FAILURE;
        }
      break; 
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:     
    default:      
      std::cout << "unknown component type" << std::endl;   
      break;      
    }    
  }  
  catch( itk::ExceptionObject &excep)   
    {    
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;    
    return EXIT_FAILURE;    
    }  
  return EXIT_SUCCESS;
}

