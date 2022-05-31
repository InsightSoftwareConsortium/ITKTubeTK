/*=========================================================================

Library:   TubeTK

Copyright 2017 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include <itkImage.h>
#include <itkImageFileWriter.h>

#include <string>

const int dimension = 2;
using PixelType = float;
using ImageType = itk::Image<PixelType, dimension>;
using VectorPixelType = itk::Vector<PixelType, dimension>;
using VectorImageType = itk::Image<VectorPixelType, dimension>;

// Create an image, with memory allocated
template <typename TImage>
typename TImage::Pointer CreateImage( typename TImage::SizeType size )
{
  auto im = TImage::New();
  im->SetRegions( size );
  im->Allocate();
  return im;
}

template <typename T>
void WriteImage( const itk::SmartPointer<T>& image, std::string file_name )
{
  auto writer = itk::ImageFileWriter<T>::New();
  writer->SetFileName( file_name );
  writer->SetInput( image );
  writer->Update();
}

int main()
{
  ImageType::SpacingType::ValueType spacing_a[] = { 1.3, 2.4 };
  ImageType::SpacingType spacing = spacing_a;

  auto in_im = CreateImage<ImageType>( {{ 15, 21 }} );
  in_im->SetSpacing( spacing );
  in_im->SetOrigin( ImageType::PointType( 0 ) );

  ImageType::SizeType out_size = {{ 3, 3 }};
  ImageType::SpacingType::ValueType out_spacing_a[] = { 6.5, 16.8 };
  ImageType::PointType::ValueType out_origin_a[] = { 2.6, 7.2 };
  ImageType::SpacingType out_spacing = out_spacing_a;
  ImageType::PointType out_origin = out_origin_a;

  auto in_vim = CreateImage<VectorImageType>( out_size );
  in_vim->SetSpacing( out_spacing );
  in_vim->SetOrigin( out_origin );

  auto out_im = CreateImage<ImageType>( out_size );
  out_im->SetSpacing( out_spacing );
  out_im->SetOrigin( out_origin );

  out_im->FillBuffer( 0 );
  out_im->SetPixel( {{ 0, 0 }}, -1 );
  out_im->SetPixel( {{ 2, 0 }}, -2 );
  out_im->SetPixel( {{ 1, 2 }}, -3 );

  for( int x = 0; x < 3; x++ )
    {
    for( int y = 0; y < 3; y++ )
      {
      PixelType v[] = { ( 6 * x + 1 ) * 1.3f, ( 6 * y + 4 ) * 2.4f };
      in_vim->SetPixel( {{ x, y }}, v );
      }
    }

  in_im->FillBuffer( 0 );
  in_im->SetPixel( {{ 1, 4 }}, -1 );
  in_im->SetPixel( {{ 13, 4 }}, -2 );
  in_im->SetPixel( {{ 7, 9 }}, 4 );
  in_im->SetPixel( {{ 7, 16 }}, -3 );

  WriteImage( in_im, "ShrinkImageTest5_In.mha" );
  WriteImage( in_vim, "ShrinkImageTest5_Points.mha" );
  WriteImage( out_im, "ShrinkImageTest5_Out.mha" );
}
