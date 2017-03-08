#include <itkImage.h>
#include <itkImageFileWriter.h>

#include <iostream>

const int d = 2;
using NumType = float;
using Image = itk::Image<NumType, d>;
using VImage = itk::Image<itk::Vector<NumType, d>, d>;

// Create an image, with memory allocated
template <typename T>
typename T::Pointer CreateImage( typename T::SizeType size )
{
  auto im = T::New();
  im->SetRegions( size );
  im->Allocate();
  return im;
}

template <typename T>
void WriteImage( typename itk::SmartPointer<T> image, const char* file_name )
{
  auto writer = itk::ImageFileWriter<T>::New();
  writer->SetFileName( file_name );
  writer->SetInput( image );
  writer->Update();
}

int main()
{
  Image::SpacingType::ValueType spacing_a[] = { 1.3, 2.4 };
  Image::SpacingType spacing = spacing_a;

  auto in_im = CreateImage<Image>( {{ 15, 21 }} );
  in_im->SetSpacing( spacing );
  in_im->SetOrigin( Image::PointType( 0 ) );

  Image::SizeType out_size = {{ 3, 3 }};
  Image::SpacingType::ValueType out_spacing_a[] = { 6.5, 16.8 };
  Image::PointType::ValueType out_origin_a[] = { 2.6, 7.2 };
  Image::SpacingType out_spacing = out_spacing_a;
  Image::PointType out_origin = out_origin_a;

  auto in_vim = CreateImage<VImage>( out_size );
  in_vim->SetSpacing( out_spacing );
  in_vim->SetOrigin( out_origin );

  auto out_im = CreateImage<Image>( out_size );
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
      NumType v[] = { ( 6 * x + 1 ) * 1.3f, ( 6 * y + 4 ) * 2.4f };
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
