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
#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// ITK Includes
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class PixelT, unsigned int DimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "ResampleCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Resample code begins
template< class PixelT, unsigned int DimensionI >
int DoIt( int argc, char *argv[] )
{
  PARSE_ARGS;

  typedef PixelT                                              InputPixelType;
  typedef itk::OrientedImage< InputPixelType, DimensionI >    InputImageType;

  typedef InputPixelType                                      OutputPixelType;
  typedef itk::OrientedImage< OutputPixelType, DimensionI >   OutputImageType;

  typedef  itk::ImageFileReader< InputImageType >          InputReaderType;
  typedef  itk::ImageFileWriter< OutputImageType >         OutputWriterType;

  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter reporter( "Resample", CLPProcessInformation );
  reporter.Start( );

  typename InputImageType::Pointer inIm;
    {
    timeCollector.Start( "LoadData" );
    typename InputReaderType::Pointer reader = InputReaderType::New( );
    reader->SetFileName( inputVolume.c_str( ) );
    reader->Update( );
    inIm = reader->GetOutput( );
    timeCollector.Stop( "LoadData" );
    }
  reporter.Report( 0.1 );

  typename InputImageType::SpacingType     inSpacing = inIm->GetSpacing( );
  typename InputImageType::PointType       inOrigin = inIm->GetOrigin( );
  typename InputImageType::SizeType        inSize = 
    inIm->GetLargestPossibleRegion( ).GetSize( );
  typename InputImageType::IndexType       inIndex = 
    inIm->GetLargestPossibleRegion( ).GetIndex( );
  typename InputImageType::DirectionType   inDirection = 
    inIm->GetDirection( );

  typename OutputImageType::SizeType       outSize;
  typename OutputImageType::SpacingType    outSpacing;
  typename OutputImageType::PointType      outOrigin;
  typename OutputImageType::IndexType      outIndex;
  typename OutputImageType::DirectionType  outDirection;
  for( unsigned int i=0; i< DimensionI; i++ )
    {
    outSpacing[i] = inSpacing[i];
    outOrigin[i] = inOrigin[i];
    outIndex[i] = inIndex[i];
    outSize[i] = inSize[i];
    }
  outDirection = inDirection;

  if( matchImage.size( ) > 1 )
    {
    typename InputReaderType::Pointer matchImReader =
      InputReaderType::New( );
    matchImReader->SetFileName( matchImage );
    matchImReader->Update( );
    outSpacing = matchImReader->GetOutput( )->GetSpacing( );
    outOrigin = matchImReader->GetOutput( )->GetOrigin( );
    outDirection = matchImReader->GetOutput( )->GetDirection( );
    outSize = matchImReader->GetOutput( )->GetLargestPossibleRegion()
                           .GetSize( );
    outIndex = matchImReader->GetOutput( )->GetLargestPossibleRegion()
                           .GetIndex( );

    reporter.Report( 0.2 );
    }

  if( spacing.size( ) > 0 )
    {
    for( unsigned int i=0; i< DimensionI; i++ )
      {
      outSpacing[i] = spacing[i];
      }
    }

  if( origin.size( ) > 0 )
    {
    for( unsigned int i=0; i< DimensionI; i++ )
      {
      outOrigin[i] = origin[i];
      }
    }

  if( index.size( ) > 0 )
    {
    for( unsigned int i=0; i< DimensionI; i++ )
      {
      outIndex[i] = index[i];
      }
    }

  if( resampleFactor.size( ) > 0 )
    {
    for( unsigned int i=0; i< DimensionI; i++ )
      {
      outSpacing[i] = outSpacing[i] / resampleFactor[i];
      }
    }

  if( makeIsotropic )
    {
    double iso = outSpacing[0];
    for( unsigned int i=1; i<DimensionI-1; i++ )
      {
      iso += outSpacing[i];
      }
    iso /= ( DimensionI-1 );
    iso += outSpacing[ DimensionI-1 ];
    iso /= 2;
    for( unsigned int i=0; i<DimensionI; i++ )
      {
      outSpacing[i] = iso;
      }
    }

  for( unsigned int i=0; i<DimensionI; i++ )
    {
    if( outSpacing[i]<=0 )
      {
      std::cerr << "ERROR: Illegal or missing output spacing specified." 
                << std::endl;
      return 1;
      }
    }

  if( matchImage.size( ) == 0 )
    {
    std::vector< double > outResampleFactor;
    outResampleFactor.resize( DimensionI );
    for( unsigned int i=0; i< DimensionI; i++ )
      {
      outResampleFactor[i] = inSpacing[i]/outSpacing[i];
      outSize[i] = static_cast<unsigned long>( inSize[i] 
                                              * outResampleFactor[i] );
      }
    }
    
  typename OutputImageType::Pointer outIm;
  {
  timeCollector.Start( "Resample" );
  reporter.Report( 0.25 );
  typedef typename itk::ResampleImageFilter< InputImageType,
          OutputImageType >             ResampleFilterType;

  typename ResampleFilterType::Pointer filter = ResampleFilterType::New( );

  filter->SetInput( inIm );

  typedef typename itk::InterpolateImageFunction< InputImageType,
          double >                      InterpType;
  typename InterpType::Pointer interp;
  if( interpolator == "BSpline" )
    {
    typedef typename itk::BSplineInterpolateImageFunction< InputImageType,
            double >                    MyInterpType;
    interp = MyInterpType::New( );
    }
  else if( interpolator == "NearestNeighbor" )
    {
    typedef typename itk::NearestNeighborInterpolateImageFunction<
            InputImageType, double >    MyInterpType;
    interp = MyInterpType::New( );
    }
  else // default = if( interpolator == "Linear" )
    {
    typedef typename itk::LinearInterpolateImageFunction<
            InputImageType, double >    MyInterpType;
    interp = MyInterpType::New( );
    }
  filter->SetInterpolator( interp );

  if( loadTransform.size() > 0 )
    {
    itk::TransformFileReader::Pointer treader = 
      itk::TransformFileReader::New();
    treader->SetFileName(loadTransform);
    treader->Update();  

    typedef itk::Transform<double, DimensionI, DimensionI> TransformType;
    typename TransformType::Pointer tfm = static_cast< TransformType * >(
      treader->GetTransformList()->front().GetPointer() );

    filter->SetTransform( tfm );
    }

  filter->SetSize( outSize );
  filter->SetOutputStartIndex( outIndex );
  filter->SetOutputOrigin( outOrigin );
  filter->SetOutputSpacing( outSpacing );
  filter->SetOutputDirection( outDirection );
  filter->SetDefaultPixelValue( 0 );
  tube::CLIFilterWatcher  watcher( filter,
                                   "Resample Filter",
                                   CLPProcessInformation,
                                   0.7,
                                   0.25,
                                   true );
  filter->Update( );

  outIm = filter->GetOutput( );
  }
  reporter.Report( 0.95 );
  timeCollector.Stop( "Resample" );

  timeCollector.Start( "Write" );
  typename OutputWriterType::Pointer  writer = OutputWriterType::New( );
  writer->SetFileName( outputVolume.c_str( ) );
  writer->SetInput( outIm );
  writer->SetUseCompression( true );
  writer->Update( );
  timeCollector.Stop( "Write" );

  /*
  itk::PluginFilterWatcher watcher1( smoothing, "Smoothing",
                                    CLPProcessInformation, 0.5, 0.0 );

  itk::PluginFilterWatcher watcher2( confidenceConnected, "Segmenting",
                                    CLPProcessInformation, 0.5, 0.5 );
  */

  timeCollector.Report( );

  reporter.End( );

  return 0;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
