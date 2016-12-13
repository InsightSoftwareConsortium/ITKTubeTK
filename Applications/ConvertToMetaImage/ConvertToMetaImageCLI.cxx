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

#include "tubeMacro.h"

#include <fstream>
#include <list>

int main( int tubeNotUsed( argc ), char * tubeNotUsed( argv )[] )
{
  std::cout << "Importing your data into the MetaImage format"
    << std::endl
    << "  is merely of matter of making a text file that"
    << std::endl
    << "  contains some info about your data and that also"
    << std::endl
    << "  lists the file( s ) that contain your data."
    << std::endl
    << std::endl
    << "This program helps you create that text file - we"
    << std::endl
    << "  designate these text files with the suffix '.mhd'"
    << std::endl
    << std::endl;

  std::string filename;
  std::cout << "What is the name of the '.mhd' file that you"
    << std::endl
    << "  want to create?   Please include '.mhd' at the end."
    << std::endl
    << "  -=> ";
  std::cin >> filename;
  std::cout << std::endl;

  std::ofstream fp;
  fp.open( filename.c_str() );

  std::cout << "What is the dimensionality of your data ( 1, 2, 3, ... )?"
    << std::endl
    << "  -=> ";
  int nDims;
  std::cin >> nDims;
  if( nDims < 1 )
    {
    std::cout << "...have a nice day" << std::endl;
    return 1;
    }
  std::cout << std::endl;
  fp << "NDims = " << nDims << std::endl;

  int i;
  int * dimSize = new int [nDims];
  std::cout << "Assuming dimension 0 is the 'x' dimension..." << std::endl;
  for( i = 0; i < nDims; i++ )
    {
    std::cout << "  How large is dimension " << i << " -=> ";
    std::cin >> dimSize[i];
    if( dimSize[i]<1 )
      {
      std::cout << "    ...we'll assume you meant 1..." << std::endl;
      dimSize[i] = 1;
      }
    }
  std::cout << std::endl;
  fp << "DimSize =";
  for( i = 0; i < nDims; i++ )
    {
    fp << " " << dimSize[i];
    }
  fp << std::endl;


  float * pointSpacing = new float [nDims];
  std::cout << "Points in an image represent a measure in physical space."
    << std::endl
    << "   Measures in physical space are done using a particular"
    << std::endl
    << "   aperature spacing - that is, what is the spacing of the"
    << std::endl
    << "   points in your data...if you don't know/care, answer 1."
    << std::endl;
  for( i = 0; i < nDims; i++ )
    {
    std::cout << "  What is the point spacing in dimension " << i
      << "? -=> ";
    std::cin >> pointSpacing[i];
    if( pointSpacing[i]<0.0000000000000001 )
      {
      std::cout << "...we'll assume you meant 1..." << std::endl;
      pointSpacing[i] = 1;
      }
    }
  std::cout << std::endl;
  fp << "ElementSpacing =";
  for( i = 0; i < nDims; i++ )
    {
    fp << " " << pointSpacing[i];
    }
  fp << std::endl;

  float * origin = new float [nDims];
  std::cout << "Images are taken of a particular location in space."
    << std::endl
    << "   That is, they have an origin in space...if you don't"
    << std::endl
    << "   know/care, answer 0 for each dimension."
    << std::endl;
  for( i = 0; i < nDims; i++ )
    {
    std::cout << "  What is the origin in dimension " << i << "? -=> ";
    std::cin >> origin[i];
    }
  std::cout << std::endl;
  fp << "Position =";
  for( i = 0; i < nDims; i++ )
    {
    fp << " " << origin[i];
    }
  fp << std::endl;

  bool byteOrderMSB = false;
  std::cout << "Is your data stored using MSB byte ordering?"
    << std::endl
    << "  Data written on a 'PC' is in 'Little-endian' byte order."
    << std::endl
    << "  Data from almost every other machine is in 'Big-endian'"
    << std::endl
    << "  byte order.  If your data was written on a 'PC', enter 0."
    << std::endl
    << "  Otherwise, if you data was written on a Mac/Sun/Other, enter 1."
    << std::endl
    << "  -=> ";
  std::cin >> byteOrderMSB;
  std::cout << std::endl;
  if( byteOrderMSB )
    {
    fp << "ElementByteOrderMSB = True" << std::endl;
    }
  else
    {
    fp << "ElementByteOrderMSB = False" << std::endl;
    }

  int nComps = 1;
  std::cout << "How many channels are stored at each point in your image?"
    << std::endl
    << "  For example, RGB images have 3 channels, and most medical"
    << std::endl
    << "  images have a single channel."
    << std::endl
    << "  -=> ";
  std::cin >> nComps;
  std::cout << std::endl;
  if( nComps<1 )
    {
    nComps = 1;
    }
  if( nComps != 1 )
    {
    fp << "ElementNumberOfChannels = " << nComps << std::endl;
    }

  int elementType;
  std::cout << "What is the 'type' of the elements in your data :"
    << std::endl
    << "  0 - signed char ( one byte )" << std::endl
    << "  1 - unsigned char" << std::endl
    << "  2 - signed short ( two byte )" << std::endl
    << "  3 - unsigned short" << std::endl
    << "  4 - signed int ( four byte )" << std::endl
    << "  5 - unsigned int" << std::endl
    << "  6 - float ( four byte )" << std::endl
    << "  7 - double ( eight byte )" << std::endl
    << "  -=> ";
  std::cin >> elementType;
  std::cout << std::endl;
  switch( elementType )
    {
    case 0: fp << "ElementType = MET_CHAR" << std::endl;
      break;
    case 1: fp << "ElementType = MET_UCHAR" << std::endl;
      break;
    case 2: fp << "ElementType = MET_SHORT" << std::endl;
      break;
    case 3: fp << "ElementType = MET_USHORT" << std::endl;
      break;
    case 4: fp << "ElementType = MET_INT" << std::endl;
      break;
    case 5: fp << "ElementType = MET_UINT" << std::endl;
      break;
    case 6: fp << "ElementType = MET_FLOAT" << std::endl;
      break;
    case 7: fp << "ElementType = MET_DOUBLE" << std::endl;
      break;
    default: fp.close();
      std::cout << "...have a nice day." << std::endl;
    }

  int storage = 0;
  std::cout << "How is the data stored ?" << std::endl
    << "  0 - in one file" << std::endl
    << "  1 - in one file per slice ( e.g., like dicom )" << std::endl
    << "  -=> ";
  std::cin >> storage;
  std::cout << std::endl;

  int headerSize = -1;
  std::cout << "Size of the header that must be skipped to"
    << std::endl
    << "  reach the data in the file( s )? Enter 0 to not skip"
    << std::endl
    << "  a header, -1 to have MetaImageReader automatically"
    << std::endl
    << "  calculate the headersize assuming the data is at"
    << std::endl
    << "  the end of the file, or enter the headersize"
    << std::endl
    << "  -=> ";
  std::cin >> headerSize;
  std::cout << std::endl;
  if( headerSize > 0 || headerSize == -1 )
    {
    fp << "HeaderSize = " << headerSize << std::endl;
    }

  if( storage == 0 )
    {
    std::string fname;
    std::cout << "Name of the file containing the data -=> ";
    std::cin >> fname;
    fp << "ElementDataFile = " << fname << std::endl;
    fp.close();
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "You're done!" << std::endl;
    std::cout << std::endl;
    }
  else
    {
    std::cout << "You've got two options when doing one file per slice:"
      << std::endl
      << "  0 - Listing the name of each slice in the .mhd file"
      << std::endl
      << "  1 - Specifying a 'fprintf' style string and min, max,"
      << std::endl
      << "        and step values for integers to be substituted"
      << std::endl
      << "        into that string to specify the numbered files"
      << std::endl
      << "        that are the slices.  For example, use data.%03d to"
      << std::endl
      << "        match data.000, data.001, data.002,...If you don't"
      << std::endl
      << "        know what we are talking about, choose option 0."
      << std::endl
      << "  -=> ";
    int storageList = 0;
    std::cin >> storageList;
    std::cout << std::endl;

    if( storageList == 0 )
      {
      fp << "ElementDataFile = LIST ";
        std::cout << "What is the dimension of the data stored in files?"
        << std::endl
        << "It must be less than or equal to nDims ( " << nDims
        << " ) Enter 2 for typical image slice data."
        <<std::endl
        << "  -=> ";
        int fileImageDim;
      std::cin >>fileImageDim;
      std::cout << std::endl;
      if( ( fileImageDim < 0 ) || ( fileImageDim > nDims ) )
        {
        std::cout << "    ...we'll assume you meant " << ( nDims - 1 )
          << " ..." << std::endl;
        }

      fp << fileImageDim
        << std::endl;
      int totalFiles = 1;
      for( i = nDims; i > fileImageDim; i-- )
        {
        totalFiles *= dimSize[i-1];
        }

      std::string str;
      for( i=0; i<totalFiles; i++ )
        {
        std::cout << "Filename of slice " << i << " -=> ";
        std::cin >> str;
        fp << str << std::endl;
        str.erase();
        }
      fp.close();
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "You're done!" << std::endl;
      std::cout << std::endl;
      }
    else
      {
      std::string storageListString;
      float storageListStringMin;
      float storageListStringMax;
      float storageListStringStep;
      std::cout << "Filename string -=> ";
      std::cin >> storageListString;
      std::cout << "  Min value to substitute into string -=> ";
      std::cin >> storageListStringMin;
      std::cout << "  Max value to substitute into string -=> ";
      std::cin >> storageListStringMax;
      std::cout << "  Step size to use to go from min to max -=> ";
      std::cin >> storageListStringStep;
      std::cout << std::endl;
      fp << "ElementDataFile = " << storageListString << " "
        << storageListStringMin << " " << storageListStringMax << " "
        << storageListStringStep << std::endl;
      fp.close();
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "You're done!" << std::endl;
      std::cout << std::endl;
      }
    }

  // Memory cleanup
  delete [] origin;
  delete [] pointSpacing;
  delete [] dimSize;
  return 0;
}
