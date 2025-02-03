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

#include "tubeBrentOptimizer1D.h"
#include "tubeCompareImageWithPrior.h"
#include "tubeSplineApproximation1D.h"
#include "tubeSplineND.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include "SimulateAcquisitionArtifactsUsingPriorCLP.h"

template <class TPixel, unsigned int VDimension>
int
DoIt(int argc, char * argv[]);

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

template <class TPixel, unsigned int VDimension>
class MyMIWPFunc : public tube::UserFunction<vnl_vector<int>, double>
{
public:
  typedef tube::CompareImageWithPrior<TPixel, VDimension> ImageEvalType;
  MyMIWPFunc(ImageEvalType & eval)
    : m_Eval(eval)
    , m_GoF(0)
  {}

  const double &
  Value(const vnl_vector<int> & x)
  {
    m_Eval.SetErode(x[0]);
    m_Eval.SetDilate(x[1]);
    m_Eval.SetGaussianBlur(x[2]);
    m_Eval.Update();
    m_GoF = m_Eval.GetGoodnessOfFit();
    return m_GoF;
  }

private:
  ImageEvalType m_Eval;
  double        m_GoF;

}; // End class MyMIWPFunc

template <class TPixel, unsigned int VDimension>
int
DoIt(int argc, char * argv[])
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter progressReporter("MatchImageWithPrior", CLPProcessInformation);
  progressReporter.Start();

  typedef float                             PixelType;
  typedef itk::Image<PixelType, VDimension> ImageType;

  /** Read input images */
  typename ImageType::Pointer inVolume;
  typename ImageType::Pointer inMask;

  timeCollector.Start("Read");
  {
    typedef itk::ImageFileReader<ImageType> ReaderType;

    typename ReaderType::Pointer readerVolume = ReaderType::New();
    typename ReaderType::Pointer readerMask = ReaderType::New();

    // read input image
    readerVolume->SetFileName(inputVolume.c_str());
    readerMask->SetFileName(inputMask.c_str());

    try
    {
      readerVolume->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      tube::ErrorMessage("Reading volume. Exception caught: " + std::string(err.GetDescription()));
      timeCollector.Report();
      return EXIT_FAILURE;
    }

    try
    {
      readerMask->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      tube::ErrorMessage("Reading mask. Exception caught: " + std::string(err.GetDescription()));
      timeCollector.Report();
      return EXIT_FAILURE;
    }

    inVolume = readerVolume->GetOutput();
    inMask = readerMask->GetOutput();
  }
  progressReporter.Report(0.1);
  timeCollector.Stop("Read");

  typename ImageType::Pointer metricMaskImage = nullptr;
  if (metricMask.size() != 0)
  {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer            readerMetricMask = ReaderType::New();
    readerMetricMask->SetFileName(metricMask.c_str());
    readerMetricMask->Update();
    metricMaskImage = readerMetricMask->GetOutput();
  }
  progressReporter.Report(0.2);

  typedef tube::CompareImageWithPrior<TPixel, VDimension> ImageEvalType;
  ImageEvalType                                           eval;
  eval.SetInput(inVolume);
  eval.SetMaskImage(inMask);
  eval.SetForeground(foreground);
  eval.SetBackground(background);
  eval.SetBoundarySize(outputBoundary);
  eval.SetNormalize(true);
  eval.SetSamplingRate(samplingRate);
  if (seed > 0)
  {
    eval.SetSeed(seed);
  }
  eval.SetTimeCollector(&timeCollector);
  eval.SetProgressReporter(&progressReporter, 0.3, 0.1);

  if (metricMaskImage.IsNotNull())
  {
    eval.SetMetricMask(metricMaskImage);
  }

  typedef typename ImageEvalType::RegistrationMethodType::TransformType TransformType;
  typename TransformType::Pointer                                       regTfm;
  if (loadTransform.size() == 0)
  {
    eval.SetUseRegistrationTransform(false);
  }
  else
  {
    itk::TransformFileReader::Pointer treader = itk::TransformFileReader::New();
    treader->SetFileName(loadTransform);
    treader->Update();
    typename TransformType::Pointer transform;
    transform = static_cast<TransformType *>(treader->GetTransformList()->front().GetPointer());

    eval.SetRegistrationTransform(transform);
    eval.SetUseRegistrationTransform(true);
  }

  if (disableRegistrationOptimization)
  {
    eval.SetUseRegistrationOptimization(false);
    if (!eval.GetUseRegistrationTransform())
    {
      eval.SetUseRegistration(false);
    }
  }

  eval.SetErode(erode);
  eval.SetDilate(dilate);
  eval.SetGaussianBlur(gaussianBlur);
  eval.Update();

  int                         outErode = erode;
  int                         outDilate = dilate;
  double                      outGaussianBlur = gaussianBlur;
  double                      outGoF = eval.GetGoodnessOfFit();
  typename ImageType::Pointer outVolume = eval.GetOutput();
  typename ImageType::Pointer outMask = eval.GetOutputMaskImage();

  if (loadTransform.size() == 0)
  {
    regTfm = eval.GetRegistrationTransform();
  }

  if (saveTransform.size() > 0)
  {
    itk::TransformFileWriter::Pointer twriter = itk::TransformFileWriter::New();
    twriter->SetInput(regTfm);
    twriter->SetFileName(saveTransform);
    twriter->Update();
  }

  if (!disableParameterOptimization)
  {
    eval.SetUseRegistrationTransform(true);
    eval.SetRegistrationTransform(regTfm);
    eval.SetUseRegistrationOptimization(false);
    eval.SetNormalize(true);
    eval.SetProgressReporter(&progressReporter, 0.4, 0.5);

    MyMIWPFunc<TPixel, VDimension> * myFunc = new MyMIWPFunc<TPixel, VDimension>(eval);
    tube::SplineApproximation1D *    spline1D = new tube::SplineApproximation1D();
    tube::BrentOptimizer1D *         opt = new tube::BrentOptimizer1D();
    tube::SplineND                   spline(3, myFunc, spline1D, opt);

    vnl_vector<int> xMin(3);
    xMin.fill(1);
    vnl_vector<int> xMax(3);
    xMax.fill(12);
    spline.SetXMin(xMin);
    spline.SetXMax(xMax);

    vnl_vector<double> x(3);
    x[0] = erode;
    x[1] = dilate;
    x[2] = gaussianBlur;

    spline.Extreme(x, &outGoF);

    outErode = x[0];
    outDilate = x[1];
    outGaussianBlur = x[2];

    std::cout << "Opt erode best = " << outErode << std::endl;
    std::cout << "Opt dilate best = " << outDilate << std::endl;
    std::cout << "Opt gaussian best = " << outGaussianBlur << std::endl;
    std::cout << "Opt gof best = " << outGoF << std::endl;

    eval.SetInput(inVolume);
    eval.SetMaskImage(inMask);
    eval.SetForeground(foreground);
    eval.SetBackground(background);
    eval.SetBoundarySize(outputBoundary);
    eval.SetNormalize(true);
    eval.SetUseRegistration(true);
    eval.SetErode(outErode);
    eval.SetDilate(outDilate);
    eval.SetGaussianBlur(outGaussianBlur);
    eval.SetProgressReporter(&progressReporter, 0.3, 0.1);
    eval.Update();
    outVolume = eval.GetOutput();
    outMask = eval.GetOutputMaskImage();

    delete myFunc;
    delete spline1D;
    delete opt;
  }

  std::cout << "Erode = " << outErode << std::endl;
  std::cout << "Dilate = " << outDilate << std::endl;
  std::cout << "Blur = " << outGaussianBlur << std::endl;

  progressReporter.Report(0.9);

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  typename ImageWriterType::Pointer       writerVolume = ImageWriterType::New();
  typename ImageWriterType::Pointer       writerMask = ImageWriterType::New();

  writerVolume->SetFileName(outputVolume.c_str());
  writerVolume->SetInput(outVolume);
  writerVolume->SetUseCompression(true);

  writerMask->SetFileName(outputMask.c_str());
  writerMask->SetInput(outMask);
  writerMask->SetUseCompression(true);

  try
  {
    writerVolume->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    tube::ErrorMessage("Writing volume. Exception caught: " + std::string(err.GetDescription()));
    timeCollector.Report();
    return EXIT_FAILURE;
  }

  try
  {
    writerMask->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    tube::ErrorMessage("Writing mask. Exception caught: " + std::string(err.GetDescription()));
    timeCollector.Report();
    return EXIT_FAILURE;
  }

  progressReporter.Report(1.0);
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt(inputVolume, argc, argv);
}
