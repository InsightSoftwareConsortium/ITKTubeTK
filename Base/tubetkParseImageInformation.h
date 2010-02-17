#ifndef __tubetkParseImageInformation_h
#define __tubetkParseImageInformation_h

#include "itkImage.h"
#include "itkImageIOBase.h"

// Description:
// Get the ComponentType and dimension of the image 
void GetImageInformation(std::string fileName,
                     itk::ImageIOBase::IOComponentType &componentType, 
                     unsigned int & dimension)
  {
  // Find out the component type of the image in file
  typedef itk::ImageIOBase::IOComponentType  PixelType;

  itk::ImageIOBase::Pointer imageIO = 
    itk::ImageIOFactory::CreateImageIO( fileName.c_str(), 
                                   itk::ImageIOFactory::ReadMode );
  if( !imageIO )
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to 
  // read the meta data from the image file.
  imageIO->SetFileName( fileName.c_str() );
  imageIO->ReadImageInformation();

  componentType = imageIO->GetComponentType();
  dimension = imageIO->GetNumberOfDimensions();
  }
#endif
