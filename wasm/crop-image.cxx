/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkImage.h"
#include "itkInputImage.h"
#include "itkOutputImage.h"
#include "itkPipeline.h"
#include "itkSupportInputImageTypes.h"
#include "itkTimeProbesCollectorBase.h"

#include "tubeMessage.h"
#include "tubeCropImage.h"

template<typename TImage>
class PipelineFunctor
{
public:
  int operator()(itk::wasm::Pipeline & pipeline)
  {
    itk::TimeProbesCollectorBase timeCollector;

    using ImageType = TImage;

    using InputImageType = itk::wasm::InputImage<ImageType>;
    InputImageType inputVolume;
    pipeline.add_option("input-volume", inputVolume, "Input volume")->required()->type_name("INPUT_IMAGE");

    InputImageType matchVolume;
    pipeline.add_option("-v,--match-volume", matchVolume, "Match crop to volume.")->type_name("INPUT_IMAGE");

    InputImageType matchMask;
    pipeline.add_option("-k,--match-mask", matchMask, "Match crop to the axis-aligned box that fits the mask.")->type_name("INPUT_IMAGE");

    using OutputImageType = itk::wasm::OutputImage<ImageType>;
    OutputImageType outputVolume;
    pipeline.add_option("output-volume", outputVolume, "Output volume")->required()->type_name("OUTPUT_IMAGE");

    std::vector<int> min;
    pipeline.add_option("-m,--min", min, "Minimum coordinate. One corner of the hyper-rectangle.");

    std::vector<int> max;
    pipeline.add_option("-M,--max", max, "Maximum coordinate (use instead of size). Another corner of the hyper-rectangle.");

    std::vector<int> size;
    pipeline.add_option("-s,--size", size, "Size of ROI (use instead of maximum). Distance to the adjacent corner of the hyper-rectangle.");

    std::vector<int> center;
    pipeline.add_option("-c,--center", center, "Center of ROI (use instead of minimum/maximum). Center of the hyper-rectangle.");

    std::vector<int> boundary;
    pipeline.add_option("-b,--boundary", boundary, "Additional boundary pixels. Add pixels beyond specified edges.");

    ITK_WASM_PARSE(pipeline);

    const unsigned int Dimension = inputVolume.Get()->GetImageDimension();

    typedef tube::CropImage< ImageType, ImageType > CropFilterType;
    typename CropFilterType::Pointer cropFilter = CropFilterType::New();

    cropFilter->SetInput( inputVolume.Get() );

    if( size.size() > 0 || max.size() > 0 || min.size() > 0 ||
      matchVolume.Get() || matchMask.Get() )
      {
      if( size.size() > 0 && max.size() > 0 )
        {
        tube::ErrorMessage(
          "You must specify either --size or --max options.  Not both." );
        return EXIT_FAILURE;
        }

      if( center.size() > 0 && min.size() > 0 )
        {
        tube::ErrorMessage(
          "You must specify either --center or --min options.  Not both." );
        return EXIT_FAILURE;
        }

      timeCollector.Start( "CropFilter" );

      if( matchVolume.Get() )
        {
        cropFilter->SetMatchVolume( matchVolume.Get() );
        }

      if( matchMask.Get() )
        {
        timeCollector.Start( "Mask Bounding Box" );

        cropFilter->SetMatchMask( matchMask.Get() );

        timeCollector.Stop( "Mask Bounding Box" );
        }

      if( min.size() > 0 )
        {
        typename ImageType::IndexType minI;
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          minI[i] = min[i];
          }
        cropFilter->SetMin( minI );
        }

      if( max.size() > 0 )
        {
        typename ImageType::IndexType maxI;
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          maxI[i] = max[i];
          }
        cropFilter->SetMax( maxI );
        }

      if( size.size() > 0 )
        {
        typename ImageType::SizeType sizeI;
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          sizeI[i] = size[i];
          }
        cropFilter->SetSize( sizeI );
        }

      if( center.size() > 0 )
        {
        typename ImageType::IndexType centerI;
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          centerI[i] = center[i];
          }
        cropFilter->SetCenter( centerI );
        }

      if( boundary.size() > 0 )
        {
        typename ImageType::IndexType boundaryI;
        for( unsigned int i = 0; i < Dimension; i++ )
          {
          boundaryI[i] = boundary[i];
          }
        cropFilter->SetBoundary( boundaryI );
        }

      try
        {
        cropFilter->Update();
        }
      catch( itk::ExceptionObject & e )
        {
        std::stringstream out;
        out << "Crop Filter: itk exception: ";
        out << e;
        tube::ErrorMessage( out.str() );
        timeCollector.Stop( "CropFilter" );
        throw( out.str() );
        }
      catch( const std::string & s )
        {
        std::cerr << "Error during crop filter: " << s << std::endl;
        timeCollector.Stop( "CropFilter" );
        return EXIT_FAILURE;
        }
      catch( ... )
        {
        std::cerr << "Error during crop filter" << std::endl;
        timeCollector.Stop( "CropFilter" );
        return EXIT_FAILURE;
        }

      timeCollector.Stop( "CropFilter" );
      timeCollector.Start( "Save data" );

      try
        {
        cropFilter->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Exception caught: " << err << std::endl;
        timeCollector.Stop( "Save data" );
        return EXIT_FAILURE;
        }
      timeCollector.Stop( "Save data" );
      }

    typename ImageType::ConstPointer output = cropFilter->GetOutput();
    outputVolume.Set( output );

    timeCollector.Report();

    return EXIT_SUCCESS;
  }
};

int main (int argc, char * argv[])
{
  itk::wasm::Pipeline pipeline("crop-image", "Extract a hyper-rectangular region from an image.", argc, argv);

  return itk::wasm::SupportInputImageTypes<PipelineFunctor,
   uint8_t,
   int8_t,
   // List desired pixel types to support
   float,
   double
   //itk::Vector<uint8_t, 3>,
   //itk::Vector<float, 3>
   >
  ::Dimensions<2U,3U>("input-volume", pipeline);
}
