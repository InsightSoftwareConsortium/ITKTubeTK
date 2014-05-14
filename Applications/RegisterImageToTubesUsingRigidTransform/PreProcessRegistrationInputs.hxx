/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __PreProcessRegistrationInputs_hxx
#define __PreProcessRegistrationInputs_hxx

#include "PreProcessRegistrationInputs.h"

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"
#include "itktubeSubSampleTubeTreeSpatialObjectFilterSerializer.h"
#include "itktubeTubeAngleOfIncidenceWeightFunction.h"
#include "itktubeTubeAngleOfIncidenceWeightFunctionSerializer.h"
#include "itktubeTubeExponentialResolutionWeightFunction.h"
#include "itktubeTubePointWeightsCalculator.h"
#include "itkUltrasoundProbeGeometryCalculator.h"
#include "itkUltrasoundProbeGeometryCalculatorSerializer.h"

#include <itkImageFileReader.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include "tubeCLIProgressReporter.h"

template<
  unsigned int VDimension,
  typename TFloat,
  typename TTube,
  typename TTubeNet,
  typename TImage,
  typename TRegistrationMethod >
int
PreProcessRegistrationInputs( int argc,
  char * argv[],
  itk::TimeProbesCollectorBase & timeCollector,
  tube::CLIProgressReporter & progressReporter,
  typename TImage::Pointer & currentImage,
  typename TTubeNet::Pointer & tubeNet,
  typename TRegistrationMethod::FeatureWeightsType & pointWeights )
{
  const unsigned int Dimension = VDimension;
  typedef TFloat                                         FloatType;
  typedef TTube                                          TubeType;
  typedef TTubeNet                                       TubeNetType;
  typedef TImage                                         ImageType;
  typedef TRegistrationMethod                            RegistrationMethodType;
  typedef typename RegistrationMethodType::TransformType TransformType;

  PARSE_ARGS;

#ifdef SlicerExecutionModel_USE_SERIALIZER
  // If SlicerExecutionModel was built with Serializer support, there is
  // automatically a parametersDeSerialize argument.  This argument is a JSON
  // file that has values for the CLI parameters, but it can also hold other
  // entries without causing any issues.
  Json::Value parametersRoot;
  if( !parametersDeSerialize.empty() )
    {
    // Parse the Json.
    std::ifstream stream( parametersDeSerialize.c_str() );
    Json::Reader reader;
    reader.parse( stream, parametersRoot );
    stream.close();
    }
#endif

  timeCollector.Start("Load data");
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typedef itk::SpatialObjectReader< Dimension > TubeNetReaderType;
  typename TubeNetReaderType::Pointer vesselReader = TubeNetReaderType::New();
  vesselReader->SetFileName( inputVessel );
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading vessel: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start("Sub-sample data");
  typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter< TubeNetType, TubeType >
                                                         SubSampleTubeTreeFilterType;
  typename SubSampleTubeTreeFilterType::Pointer subSampleTubeTreeFilter =
    SubSampleTubeTreeFilterType::New();
  subSampleTubeTreeFilter->SetInput( vesselReader->GetGroup() );
  subSampleTubeTreeFilter->SetSampling( 100 );
#ifdef SlicerExecutionModel_USE_SERIALIZER
  if( !parametersDeSerialize.empty() )
    {
    // If the Json file has entries that describe the parameters for an
    // itk::tube::SubSampleTubeTreeSpatialObjectFilter, read them in, and set them on our
    // instance.
    if( parametersRoot.isMember( "SubSampleTubeTree" ) )
      {
      Json::Value & subSampleTubeTreeFilterValue = parametersRoot["SubSampleTubeTree"];
      typedef itk::tube::SubSampleTubeTreeSpatialObjectFilterSerializer<
        SubSampleTubeTreeFilterType > SerializerType;
      typename SerializerType::Pointer serializer = SerializerType::New();
      serializer->SetTargetObject( subSampleTubeTreeFilter );
      itk::JsonCppArchiver::Pointer archiver =
        dynamic_cast< itk::JsonCppArchiver * >( serializer->GetArchiver() );
      archiver->SetJsonValue( &subSampleTubeTreeFilterValue );
      serializer->DeSerialize();
      }
    }
#endif
  try
    {
    subSampleTubeTreeFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Sub-sampling vessel: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  tubeNet = subSampleTubeTreeFilter->GetOutput();

  timeCollector.Stop("Sub-sample data");
  progress = 0.2;
  progressReporter.Report( progress );


  currentImage = reader->GetOutput();
  if( gaussianBlurStdDev > 0.0 )
    {
    timeCollector.Start("Gaussian Blur");

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
      GaussianFilterType;
    typename GaussianFilterType::Pointer gaussianFilter;

    // Progress per iteration
    const double progressFraction = 0.1/Dimension;
    for( unsigned int ii = 0; ii < Dimension; ++ii )
      {
      gaussianFilter = GaussianFilterType::New();
      gaussianFilter->SetInput( currentImage );
      gaussianFilter->SetSigma( gaussianBlurStdDev );

      gaussianFilter->SetOrder(
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      gaussianFilter->SetDirection( ii );
      tube::CLIFilterWatcher watcher( gaussianFilter,
                                      "Blur Filter 1D",
                                      CLPProcessInformation,
                                      progressFraction,
                                      progress,
                                      true );

      gaussianFilter->Update();
      currentImage = gaussianFilter->GetOutput();
      }

    timeCollector.Stop("Gaussian Blur");
    }

  timeCollector.Start("Compute Model Feature Weights");
  typedef itk::tube::Function::TubeExponentialResolutionWeightFunction<
    typename TubeType::TubePointType, FloatType >    WeightFunctionType;
  typedef typename RegistrationMethodType::FeatureWeightsType
    PointWeightsType;
  typename WeightFunctionType::Pointer weightFunction =
    WeightFunctionType::New();
  typedef itk::tube::TubePointWeightsCalculator< Dimension,
    TubeType, WeightFunctionType,
    PointWeightsType > PointWeightsCalculatorType;
  typename PointWeightsCalculatorType::Pointer resolutionWeightsCalculator
    = PointWeightsCalculatorType::New();
  resolutionWeightsCalculator->SetTubeTreeSpatialObject(
    subSampleTubeTreeFilter->GetOutput() );
  resolutionWeightsCalculator->SetPointWeightFunction( weightFunction );
  resolutionWeightsCalculator->Compute();
  pointWeights = resolutionWeightsCalculator->GetPointWeights();

#ifdef SlicerExecutionModel_USE_SERIALIZER
  // Compute ultrasound probe geometry.
  if( !parametersDeSerialize.empty() )
    {
    // If the Json file has entries that describe the parameters for an
    // itk::tube::SubSampleTubeTreeSpatialObjectFilter, read them in, and set them on our
    // instance.
    if( parametersRoot.isMember( "UltrasoundProbeGeometryCalculator" ) )
      {
      timeCollector.Start("Compute probe geometry");
      Json::Value & probeGeometryCalculatorValue
        = parametersRoot["UltrasoundProbeGeometryCalculator"];

      typedef itk::tube::UltrasoundProbeGeometryCalculator< ImageType >
        GeometryCalculatorType;
      typename GeometryCalculatorType::Pointer geometryCalculator
        = GeometryCalculatorType::New();
      geometryCalculator->SetInput( reader->GetOutput() );

      typedef itk::tube::UltrasoundProbeGeometryCalculatorSerializer<
        GeometryCalculatorType > SerializerType;
      typename SerializerType::Pointer serializer =
        SerializerType::New();
      serializer->SetTargetObject( geometryCalculator );
      itk::JsonCppArchiver::Pointer archiver =
        dynamic_cast< itk::JsonCppArchiver * >( serializer->GetArchiver() );
      archiver->SetJsonValue( &probeGeometryCalculatorValue );
      serializer->DeSerialize();

      try
        {
        geometryCalculator->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        tube::ErrorMessage( "Computing probe geometry: Exception caught: "
                            + std::string(err.GetDescription()) );
        timeCollector.Report();
        return EXIT_FAILURE;
        }

      if( parametersRoot.isMember( "UltrasoundProbeGeometryFile" ) )
        {
        const char * outputFile
          = parametersRoot["UltrasoundProbeGeometryFile"].asCString();
        std::ofstream geometryOutput( outputFile );
        if( !geometryOutput.is_open() )
          {
          tube::ErrorMessage( "Could not open geometry output file: "
                              + std::string( outputFile ) );
          timeCollector.Report();
          return EXIT_FAILURE;
          }

        const typename GeometryCalculatorType::OriginType ultrasoundProbeOrigin =
          geometryCalculator->GetUltrasoundProbeOrigin();
        geometryOutput << "UltrasoundProbeOrigin:";
        for( unsigned int ii = 0; ii < Dimension; ++ii )
          {
          geometryOutput << " " << ultrasoundProbeOrigin[ii];
          }
        geometryOutput << std::endl;

        const typename GeometryCalculatorType::RadiusType startOfAcquisitionRadius =
          geometryCalculator->GetStartOfAcquisitionRadius();
        geometryOutput << "GetStartOfAcquisitionRadius: "
                       << startOfAcquisitionRadius
                       << std::endl;
        }

      if( parametersRoot.isMember( "AngleOfIncidenceWeightFunction" ) )
        {
        Json::Value & angleOfIncidenceWeightFunctionValue =
          parametersRoot["AngleOfIncidenceWeightFunction"];

        typedef itk::tube::Function::TubeAngleOfIncidenceWeightFunction<
          typename TubeType::TubePointType, FloatType > AngleOfIncidenceWeightFunctionType;
        typename AngleOfIncidenceWeightFunctionType::Pointer angleOfIncidenceWeightFunction =
          AngleOfIncidenceWeightFunctionType::New();

        typedef itk::tube::TubeAngleOfIncidenceWeightFunctionSerializer<
          AngleOfIncidenceWeightFunctionType >
            AngleOfIncidenceSerializerType;
        typename AngleOfIncidenceSerializerType::Pointer angleOfIncidenceSerializer =
          AngleOfIncidenceSerializerType::New();
        angleOfIncidenceSerializer->SetTargetObject( angleOfIncidenceWeightFunction );
        itk::JsonCppArchiver::Pointer angleOfIncidenceArchiver =
          dynamic_cast< itk::JsonCppArchiver * >(
            angleOfIncidenceSerializer->GetArchiver() );
        angleOfIncidenceArchiver->SetJsonValue( &angleOfIncidenceWeightFunctionValue );
        angleOfIncidenceSerializer->DeSerialize();

        angleOfIncidenceWeightFunction->SetUltrasoundProbeOrigin(
          geometryCalculator->GetUltrasoundProbeOrigin() );

        typedef itk::tube::TubePointWeightsCalculator< Dimension,
          TubeType, AngleOfIncidenceWeightFunctionType, PointWeightsType >
            AngleOfIncidenceWeightsCalculatorType;
        typename AngleOfIncidenceWeightsCalculatorType::Pointer angleOfIncidenceWeightsCalculator =
          AngleOfIncidenceWeightsCalculatorType::New();
        angleOfIncidenceWeightsCalculator->SetPointWeightFunction(
          angleOfIncidenceWeightFunction );
        angleOfIncidenceWeightsCalculator->SetTubeTreeSpatialObject(
          subSampleTubeTreeFilter->GetOutput() );
        angleOfIncidenceWeightsCalculator->Compute();
        const PointWeightsType & angleOfIncidenceWeights =
          angleOfIncidenceWeightsCalculator->GetPointWeights();
        for( itk::SizeValueType ii = 0; ii < pointWeights.GetSize(); ++ii )
          {
          pointWeights[ii] *= angleOfIncidenceWeights[ii];
          }
        }

      timeCollector.Stop("Compute probe geometry");
      }
    }

  if( parametersRoot.isMember( "TubePointWeightsFile" ) )
    {
    Json::Value weightsJSONRoot;
    weightsJSONRoot["TubePointWeights"] = Json::Value( Json::arrayValue );
    Json::Value & weightsJSON = weightsJSONRoot["TubePointWeights"];
    weightsJSON.resize( pointWeights.GetSize() );
    for( itk::SizeValueType ii = 0; ii < pointWeights.GetSize(); ++ii )
      {
      weightsJSON[static_cast<Json::ArrayIndex>(ii)] = pointWeights[ii];
      }

    Json::FastWriter writer;
    const std::string weightsString = writer.write( weightsJSONRoot );

    Json::Value & tubePointWeightsFileValue =
      parametersRoot["TubePointWeightsFile"];
    std::ofstream tubePointWeightsFile( tubePointWeightsFileValue.asCString() );
    if( !tubePointWeightsFile.is_open() )
      {
      tube::ErrorMessage( "Could not open tube point weights file: "
                          + tubePointWeightsFileValue.asString() );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    tubePointWeightsFile << weightsString;
    }
#endif
  timeCollector.Stop("Compute Model Feature Weights");

  return EXIT_SUCCESS;
}

#endif // End !defined(__PreProcessRegistrationInputs_hxx)
