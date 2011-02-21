#include "sitkIMTypeImageWriter.h"

// SimpleITK includes
#include "sitkCastImageFilter.h"

// ITK includes
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
TypeImageWriter::TypeImageWriter()
{
  m_Type = UCharType;
  m_OutputFilename = "";
}


//
// ToString
//
std::string TypeImageWriter::ToString() const
{
  std::stringstream out;
  out << "Type Image Writer" << std::endl
      << "  Type: " << m_Type << std::endl
      << "  Output Filename: " << m_OutputFilename
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetType / SetType
//
unsigned int TypeImageWriter::GetType()
{
  return this->m_Type;
}
TypeImageWriter::Self& TypeImageWriter::SetType(unsigned int t)
{
  this->m_Type = t;
  return *this;
}

//
// GetOutputFilename / SetOutputFilename
//
std::string TypeImageWriter::GetOutputFilename()
{
  return this->m_OutputFilename;
}
TypeImageWriter::Self& TypeImageWriter::SetOutputFilename(std::string f)
{
  this->m_OutputFilename = f;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
void TypeImageWriter::Execute(itk::simple::Image* image,
                              unsigned int type, std::string filename)
{
  this->SetType(type);
  this->SetOutputFilename(filename);
  this->Execute(image);
}

void TypeImageWriter::Execute(itk::simple::Image* image)
{
  // Dispatch ExecuteInternalDim with the right dimension
  if (image->GetDimension() == 2)
    {
    this->ExecuteInternalDim<2>( image );
    }
  else if (image->GetDimension() == 3)
    {
    this->ExecuteInternalDim<3>( image );
    }
  else
    {
    sitkExceptionMacro("Unknown Dimension " << image->GetDimension());
    }
}


//----------------------------------------------------------------------------
// ExecuteInternalDim
//
template <unsigned int TDimension>
void TypeImageWriter::ExecuteInternalDim( itk::simple::Image* inImage )
{
  // ITK typedefs
  typedef itk::Image< unsigned char, TDimension >  ImageTypeUChar;
  typedef itk::Image< unsigned short, TDimension > ImageTypeUShort;
  typedef itk::Image< short, TDimension >          ImageTypeShort;
  typedef itk::Image< float, TDimension >          ImageTypeFloat;

  // Set up the cast filter
  itk::simple::CastImageFilter castFilter;

  // Switch over type
  switch( this->m_Type )
    {
    case UCharType:
    case UCharUCMPType:
      {
      // Cast the sitk image to the right type
      castFilter.SetOutputPixelType( itk::simple::sitkUInt8 );
      std::auto_ptr<itk::simple::Image> castImage( castFilter.Execute( inImage ) );

      // Get the ITK image
      typename ImageTypeUChar::Pointer image =
        dynamic_cast <ImageTypeUChar*> ( castImage->GetImageBase().GetPointer() );
      if ( image.IsNull() )
        {
        sitkExceptionMacro( "Couldn't save image as unsigned char" );
        }

      // Write out the image
      typedef itk::ImageFileWriter< ImageTypeUChar >
        VolumeWriterType;
      typename VolumeWriterType::Pointer writer =
        VolumeWriterType::New();
      writer->SetFileName( this->m_OutputFilename.c_str() );
      writer->SetInput( image );

      // Only turn compression on if called with the compressed version flag
      if( this->m_Type == UCharType )
        {
        writer->SetUseCompression( true );
        }
      writer->Write();
      break;
      }
    case UShortType:
    case UShortUCMPType:
      {
      // Cast the sitk image to the right type
      castFilter.SetOutputPixelType( itk::simple::sitkUInt16 );
      std::auto_ptr<itk::simple::Image> castImage( castFilter.Execute( inImage ) );

      // Get the ITK image
      typename ImageTypeUShort::Pointer image =
        dynamic_cast <ImageTypeUShort*> ( castImage->GetImageBase().GetPointer() );
      if ( image.IsNull() )
        {
        sitkExceptionMacro( "Couldn't save image as unsigned short" );
        }

      // Write out the image
      typedef itk::ImageFileWriter< ImageTypeUShort >
        VolumeWriterType;
      typename VolumeWriterType::Pointer writer =
        VolumeWriterType::New();
      writer->SetFileName( this->m_OutputFilename.c_str() );
      writer->SetInput( image );

      // Only turn compression on if called with the compressed version flag
      if( this->m_Type == UShortType )
        {
        writer->SetUseCompression( true );
        }
      writer->Write();
      break;
      }
    case ShortType:
    case ShortUCMPType:
      {
      // Cast the sitk image to the right type
      castFilter.SetOutputPixelType( itk::simple::sitkInt16 );
      std::auto_ptr<itk::simple::Image> castImage( castFilter.Execute( inImage ) );

      // Get the ITK image
      typename ImageTypeShort::Pointer image =
        dynamic_cast <ImageTypeShort*> ( castImage->GetImageBase().GetPointer() );
      if ( image.IsNull() )
        {
        sitkExceptionMacro( "Couldn't save image as short" );
        }

      // Write out the image
      typedef itk::ImageFileWriter< ImageTypeShort >
        VolumeWriterType;
      typename VolumeWriterType::Pointer writer =
        VolumeWriterType::New();
      writer->SetFileName( this->m_OutputFilename.c_str() );
      writer->SetInput( image );

      // Only turn compression on if called with the compressed version flag
      if( this->m_Type == ShortType )
        {
        writer->SetUseCompression( true );
        }
      writer->Write();
      break;
      }
    case ShortOldType:
      {
      // Cast the sitk image to the right type
      castFilter.SetOutputPixelType( itk::simple::sitkInt16 );
      std::auto_ptr<itk::simple::Image> castImage( castFilter.Execute( inImage ) );

      // Get the ITK image
      typename ImageTypeShort::Pointer image =
        dynamic_cast <ImageTypeShort*> ( castImage->GetImageBase().GetPointer() );
      if ( image.IsNull() )
        {
        sitkExceptionMacro( "Couldn't save image as short (old)" );
        }

      // Write out the image
      typedef itk::ImageFileWriter< ImageTypeShort >
        VolumeWriterType;
      typename VolumeWriterType::Pointer writer =
        VolumeWriterType::New();

      itk::MetaImageIO::Pointer metaWriter = itk::MetaImageIO::New();
      writer->SetImageIO( metaWriter );

      writer->SetFileName( this->m_OutputFilename.c_str() );
      writer->SetInput( image );
      writer->SetUseCompression( false );

      MetaImage * metaImage = metaWriter->GetMetaImagePointer();

      metaImage->ElementSize( 0, image->GetSpacing()[0] );
      metaImage->ElementSize( 1, image->GetSpacing()[1] );
      metaImage->ElementSize( 2, image->GetSpacing()[2] );

      metaImage->AddUserField( "ElementByteOrderMSB",
                              MET_STRING, strlen( "False" ), "False" );

      writer->Write();
      break;
      }
    case FloatType:
      {
      // Cast the sitk image to the right type
      castFilter.SetOutputPixelType( itk::simple::sitkFloat32 );
      std::auto_ptr<itk::simple::Image> castImage( castFilter.Execute( inImage ) );

      // Get the ITK image
      typename ImageTypeFloat::Pointer image =
        dynamic_cast <ImageTypeFloat*> ( castImage->GetImageBase().GetPointer() );
      if ( image.IsNull() )
        {
        sitkExceptionMacro( "Couldn't save image as float" );
        }

      // Write out the image
      typedef itk::ImageFileWriter< ImageTypeFloat >
        VolumeWriterType;
      typename VolumeWriterType::Pointer writer =
        VolumeWriterType::New();
      writer->SetFileName( this->m_OutputFilename.c_str() );
      writer->SetInput( image );
      writer->Write();
      break;
      }
    }
}

} // end namespace sitkIM
