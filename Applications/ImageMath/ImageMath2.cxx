/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: main2.cxx,v $
  Language:  C++
  Date:      $Date: 2005/04/07 03:57:49 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <cmath>
#include <iostream>

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaImageIO.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNJetImageFunction.h"
#include "itkNormalVariateGenerator.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkErodeObjectMorphologyImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkDilateObjectMorphologyImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkCVTImageFilter.h"

int usage(void)
{
  std::cout << "ImageMath <infile> [options]  " << std::endl;
  std::cout << "  - note: options are applied in cmdline order" << std::endl;
  std::cout << "Options :" << std::endl;
  std::cout << "  -w <imageOutputFile>"
            << std::endl
            << "     : writes the current image to the designated file." 
            << std::endl
            << "       Image is stored as floats" 
            << std::endl;
  std::cout << "  -W <type> <imageOutputFile>"
            << std::endl
            << "     : writes the current image to the designated file." 
            << std::endl
            << "       Using pixel type (0=unsigned byte, 1=unsigned" 
            << std::endl
            << "       short, 2=short, 3=formatted for MIDAS)" 
            << std::endl;
  std::cout << "  -C <image2>"
            << std::endl
            << "     : concatenates iamge2 to the right of image" 
            << std::endl;
  return 0;
}

// Declare the image
typedef float                         PixelType;
typedef itk::Image<PixelType, 3>      ImageType;

typedef itk::Image<unsigned char, 3>  ImageTypeUChar;
typedef itk::Image<unsigned short, 3> ImageTypeUShort;
typedef itk::Image<short, 3>          ImageTypeShort;

int main(int argc, char **argv)
{
  typedef itk::Statistics::NormalVariateGenerator gaussGenType;
  gaussGenType::Pointer gaussGen = gaussGenType::New();

  if(argc < 4)
    {
    return usage();
    }

  typedef itk::ImageFileReader< ImageType >   VolumeReaderType;
  typedef itk::NJetImageFunction< ImageType > ImageFunctionType;


  // Declare a reader
  VolumeReaderType::Pointer reader = VolumeReaderType::New();
  reader->SetFileName( argv[1] );
  ImageType::Pointer imIn;
  imIn = reader->GetOutput();
  
  // See if the file can be read - "try" otherwise program will 
  //   mysteriously exit on failure in the Object factory
  std::cout << "Reading file (" << argv[1] <<")"<< std::endl;
  try
    {
    reader->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format" << std::endl;
    return 1;
    }

  int argNum = 2;
  while(argNum<argc && argv[argNum][0] == '-')
    {
    switch(argv[argNum][1])
      {
      case 'w' :
        {
        argNum++;
        const char* outFilename = argv[argNum++];
        std::cout << "Writing output (" << outFilename << ")" << std::endl;
        typedef itk::ImageFileWriter< ImageType > VolumeWriterType;
        VolumeWriterType::Pointer writer = VolumeWriterType::New();
        writer->SetFileName( outFilename );
        writer->SetInput( imIn );
        writer->SetUseCompression(true);
        writer->Write();
        break;
        } // end -w
      case 'W' :
        {
        argNum++;
        int type = (int)atof(argv[argNum++]);
        const char* outFilename = argv[argNum++];
        std::cout << "Writing output (" << outFilename << ")" << std::endl;
        switch(type)
          {
          case 0:
            {
            typedef itk::CastImageFilter< ImageType, ImageTypeUChar> 
                  CastFilterType;
            CastFilterType::Pointer castFilter = CastFilterType::New();
            castFilter->SetInput( imIn );

            typedef itk::ImageFileWriter< ImageTypeUChar > VolumeWriterType;
            VolumeWriterType::Pointer writer = VolumeWriterType::New();
            writer->SetFileName( outFilename );
            writer->SetInput( castFilter->GetOutput() );
            writer->SetUseCompression(true);
            writer->Write();
            break;
            }
          case 1:
            {
            typedef itk::CastImageFilter< ImageType, ImageTypeUShort> 
                  CastFilterType;
            CastFilterType::Pointer castFilter = CastFilterType::New();
            castFilter->SetInput( imIn );

            typedef itk::ImageFileWriter< ImageTypeUShort > VolumeWriterType;
            VolumeWriterType::Pointer writer = VolumeWriterType::New();
            writer->SetFileName( outFilename );
            writer->SetInput( castFilter->GetOutput() );
            writer->SetUseCompression(true);
            writer->Write();
            break;
            }
          case 2:
            {
            typedef itk::CastImageFilter< ImageType, ImageTypeShort> 
                  CastFilterType;
            CastFilterType::Pointer castFilter = CastFilterType::New();
            castFilter->SetInput( imIn );

            typedef itk::ImageFileWriter< ImageTypeShort > VolumeWriterType;
            VolumeWriterType::Pointer writer = VolumeWriterType::New();
            writer->SetFileName( outFilename );
            writer->SetInput( castFilter->GetOutput() );
            writer->SetUseCompression(true);
            writer->Write();
            break;
            }
          case 3:
            {
            typedef itk::CastImageFilter< ImageType, ImageTypeShort> 
                  CastFilterType;
            CastFilterType::Pointer castFilter = CastFilterType::New();
            castFilter->SetInput( imIn );

            typedef itk::ImageFileWriter< ImageTypeShort > VolumeWriterType;
            VolumeWriterType::Pointer writer = VolumeWriterType::New();

            itk::MetaImageIO::Pointer metaWriter = itk::MetaImageIO::New();
            writer->SetImageIO( metaWriter );

            writer->SetFileName( outFilename );
            writer->SetInput( castFilter->GetOutput() );
            writer->SetUseCompression(false);

            MetaImage * metaImage = metaWriter->GetMetaImagePointer();

            metaImage->ElementSize(0, imIn->GetSpacing()[0]);
            metaImage->ElementSize(1, imIn->GetSpacing()[1]);
            metaImage->ElementSize(2, imIn->GetSpacing()[2]);

            metaImage->AddUserField("ElementByteOrderMSB",
                                    MET_STRING, strlen("False"), "False");

            writer->Write();
            break;
            }
          }
        break;
        } // end -W
      case 'c' :
        {
        std::cout << "Intensity windowing" << std::endl;
        argNum++;
        VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
        reader2->SetFileName( argv[argNum++] );
        ImageType::Pointer imIn2;
        imIn2 = reader2->GetOutput();
        try
          {
          reader2->Update();
          }
        catch( ... )
          {
          std::cout << "Problems reading file format of inFile2." << std::endl;
          return 1;
          }
        itk::ImageRegionIterator< ImageType > it2(imIn2, 
                       imIn2->GetLargestPossibleRegion());
        ImageType::SizeType size;
        size = imIn->GetLargestPossibleRegion().GetSize();
        ImageType::SizeType size2;
        size2 = imIn2->GetLargestPossibleRegion().GetSize();
        ImageType::SizeType outSize;
        outSize[0] = size[0];
        outSize[1] = size[1];
        outSize[2] = size[2] + size2[2];
        ImageType::Pointer outIm = ImageType::New();
        outIm->SetRegions(outSize);
        outIm->SetSpacing(imIn->GetSpacing());
        outIm->Allocate();

        double tf;
        ImageType::IndexType indx;
        for(int i=0; i<(int)size[2]; i++)
          {
          indx[2] = i;
          for(int j=0; j<(int)size[1]; j++)
            {
            indx[1] = j;
            for(int k=0; k<(int)size[0]; k++)
              {
              indx[0] = k;
              outIm->SetPixel(indx, imIn->GetPixel(indx));
              }
            for(int k=0; k<(int)size2[0]; k++)
              {
              indx[0] = k;
              tf = imIn2->GetPixel(indx);
              indx[0] = k+size[0];
              outIm->SetPixel(indx, tf);
              }
            }
          }

        imIn = outIm;
        } // end -c
      default :
        {
        return usage();
        }
      } // end switch
    } // end while

  return 1;
}
