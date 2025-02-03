/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
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

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include <tubeMessage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFileReader.h>

#include <tubeResampleImage.h>
#include "ResampleImageCLP.h"

template <class TPixel, unsigned int VDimension>
int
DoIt(int argc, char * argv[]);

// Must follow include of "...CLP.h" and forward declaration of
//   int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

template <class TPixel, unsigned int DimensionI>
int
DoIt(int argc, char * argv[])
{
  PARSE_ARGS;

  typedef itk::Image<TPixel, DimensionI> ImageType;

  typedef typename tube::ResampleImage<ImageType> FilterType;
  typename FilterType::Pointer                    filter = FilterType::New();

  typedef typename FilterType::TransformType TransformType;

  typedef itk::ImageFileReader<ImageType> InputReaderType;
  typedef itk::ImageFileWriter<ImageType> OutputWriterType;

  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter reporter("Resample", CLPProcessInformation);
  reporter.Start();

  timeCollector.Start("LoadData");

  typename InputReaderType::Pointer reader = InputReaderType::New();
  reader->SetFileName(inputVolume.c_str());
  try
  {
    reader->Update();
    filter->SetInput(reader->GetOutput());
  }
  catch (itk::ExceptionObject & err)
  {
    tube::ErrorMessage("Reading input image. Exception caught: " + std::string(err.GetDescription()));
    timeCollector.Report();
    return EXIT_FAILURE;
  }

  timeCollector.Stop("LoadData");
  reporter.Report(0.1);

  if (matchImage.size() > 1)
  {
    typename InputReaderType::Pointer matchImReader = InputReaderType::New();
    matchImReader->SetFileName(matchImage);
    matchImReader->Update();
    filter->SetMatchImage(matchImReader->GetOutput());
    reporter.Report(0.2);
  }

  if (spacing.size() > 0)
  {
    filter->SetSpacing(spacing);
  }

  if (origin.size() > 0)
  {
    filter->SetOrigin(origin);
  }

  if (index.size() > 0)
  {
    filter->SetIndex(index);
  }

  if (resampleFactor.size() > 0)
  {
    filter->SetResampleFactor(resampleFactor);
  }

  filter->SetMakeIsotropic(makeIsotropic);
  filter->SetMakeHighResIso(makeHighResIso);
  filter->SetInterpolator(interpolator);

  if (loadTransform.size() > 0)
  {
    typename itk::TransformFileReader::Pointer treader = itk::TransformFileReader::New();
    treader->SetFileName(loadTransform);
    treader->Update();

    typename TransformType::Pointer tfm =
      static_cast<TransformType *>(treader->GetTransformList()->front().GetPointer());

    filter->SetTransform(tfm);
  }

  typename ImageType::Pointer outIm;

  timeCollector.Start("Resample");
  reporter.Report(0.25);

  filter->Update();

  outIm = filter->GetOutput();

  reporter.Report(0.95);
  timeCollector.Stop("Resample");

  timeCollector.Start("Write");
  typename OutputWriterType::Pointer writer = OutputWriterType::New();
  writer->SetFileName(outputVolume.c_str());
  writer->SetInput(outIm);
  writer->SetUseCompression(true);
  writer->Update();
  timeCollector.Stop("Write");

  timeCollector.Report();
  reporter.End();

  return 0;
}

// Main
int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt(inputVolume, argc, argv);
}
