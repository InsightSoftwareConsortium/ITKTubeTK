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


#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following four should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "tubeStringUtilities.h"
#include "itkTubeLDAGenerator.h"
#include "itkTubeMetaLDA.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "tubeNJetLDAGeneratorCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

template < class imageT >
void WriteLDA( typename imageT::Pointer img,
  std::string base, std::string ext, int num )
{
  typedef itk::ImageFileWriter< imageT >     LDAImageWriterType;

  typename LDAImageWriterType::Pointer ldaImageWriter =
    LDAImageWriterType::New();
  std::string fname = base;
  char c[80];
  sprintf( c, ext.c_str(), num );
  fname += std::string( c );
  ldaImageWriter->SetFileName( fname.c_str() );
  ldaImageWriter->SetInput( img );
  ldaImageWriter->Update();
}

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  typedef pixelT                                   InputPixelType;
  typedef itk::OrientedImage< InputPixelType, dimensionT >
                                                   InputImageType;
  typedef itk::OrientedImage< unsigned short, dimensionT >
                                                   MaskImageType;
  typedef itk::OrientedImage< float, dimensionT >  LDAImageType;

  typedef itk::ImageFileReader< LDAImageType >     ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< LDAImageType >     LDAImageWriterType;

  typedef itk::tube::LDAGenerator< LDAImageType, MaskImageType >
    LDAGeneratorType;
  typename LDAGeneratorType::Pointer ldaGenerator = LDAGeneratorType::New();

  timeCollector.Start( "LoadData" );

  unsigned int fCount = 0;
  std::vector< unsigned int > featureSymmetry;
  std::vector< std::string > featureName;
  std::vector< unsigned int > zeroOrderFeatures;
  typename ImageReaderType::Pointer reader;
  for( unsigned int vNum=0; vNum<inputVolumesList.size(); vNum++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[vNum].c_str() );
    std::string featureBaseName = outputBase;
    char s[80];
    sprintf(s, ".img%02d", vNum );
    featureBaseName += std::string( s );
    reader->Update();
    if( vNum == 0 )
      {
      ldaGenerator->SetFeatureImage( reader->GetOutput() );
      sprintf(s, ".f%02d.org", fCount);
      featureName.push_back( featureBaseName + std::string( s ) );
      featureSymmetry.push_back( 0 );
      zeroOrderFeatures.push_back( fCount );
      if( saveFeatureImages.size() > 0 )
        {
        WriteLDA< LDAImageType >( reader->GetOutput(), saveFeatureImages,
          ".f%02d.mha", fCount );
        }
      ++fCount;
      }
    else
      {
      ldaGenerator->AddFeatureImage( reader->GetOutput() );
      sprintf(s, ".f%02d.org", fCount);
      featureName.push_back( featureBaseName + std::string( s ) );
      featureSymmetry.push_back( 0 );
      zeroOrderFeatures.push_back( fCount );
      if( saveFeatureImages.size() > 0 )
        {
        WriteLDA< LDAImageType >( reader->GetOutput(), saveFeatureImages,
          ".f%02d.mha", fCount );
        }
      ++fCount;
      }
    for( unsigned int i=0; i<zeroScales.size(); i++ )
      {
      typedef itk::DiscreteGaussianImageFilter< LDAImageType,
        LDAImageType > FType;
      typename FType::Pointer filter = FType::New();
      filter->SetInput( reader->GetOutput() );
      filter->SetVariance( zeroScales[i]*zeroScales[i] );
      filter->Update();
      ldaGenerator->AddFeatureImage( filter->GetOutput() );
      sprintf(s, ".f%02d.b-%02d", fCount, (int)zeroScales[i]);
      featureName.push_back( featureBaseName + std::string( s ) );
      featureSymmetry.push_back( 0 );
      zeroOrderFeatures.push_back( fCount );
      if( saveFeatureImages.size() > 0 )
        {
        WriteLDA< LDAImageType >( filter->GetOutput(), saveFeatureImages,
          ".f%02d-0.mha", fCount );
        }
      ++fCount;
      }
    for( unsigned int i=0; i<firstScales.size(); i++ )
      {
      int symmetryBase = fCount;
      for( unsigned int d=0; d<dimensionT; d++ )
        {
        typename LDAImageType::Pointer curImage = reader->GetOutput();
        for( unsigned int d2=0; d2<dimensionT; d2++ )
          {
          typedef itk::RecursiveGaussianImageFilter< LDAImageType,
            LDAImageType > FType;
          typename FType::Pointer filter = FType::New();
          filter->SetInput( curImage );
          filter->SetDirection( d2 );
          if( d == d2 )
            {
            filter->SetSigma( firstScales[i] );
            filter->SetOrder( FType::FirstOrder );
            }
          else
            {
            filter->SetSigma( firstScales[i]/2.0 );
            filter->SetOrder( FType::ZeroOrder );
            }
          filter->Update();
          curImage = filter->GetOutput();
          }
        ldaGenerator->AddFeatureImage( curImage );
        sprintf(s, ".f%02d.d1-%02d-%01d", fCount, (int)firstScales[i],
          d );
        featureName.push_back( featureBaseName + std::string( s ) );
        featureSymmetry.push_back( symmetryBase );
        if( saveFeatureImages.size() > 0 )
          {
          WriteLDA< LDAImageType >( curImage, saveFeatureImages,
            ".f%02d-1.mha", fCount );
          }
        ++fCount;
        }
      typename LDAImageType::Pointer curImage = reader->GetOutput();
      for( unsigned int d=0; d<dimensionT; d++ )
        {
        typedef itk::RecursiveGaussianImageFilter< LDAImageType,
          LDAImageType > FType;
        typename FType::Pointer filter = FType::New();
        filter->SetInput( curImage );
        filter->SetDirection( d );
        filter->SetSigma( firstScales[i] );
        filter->SetOrder( FType::FirstOrder );
        filter->Update();
        curImage = filter->GetOutput();
        }
      ldaGenerator->AddFeatureImage( curImage );
      sprintf(s, ".f%02d.d1-%02d-a", fCount, (int)firstScales[i] );
      featureName.push_back( featureBaseName + std::string( s ) );
      featureSymmetry.push_back( 0 );
      if( saveFeatureImages.size() > 0 )
        {
        WriteLDA< LDAImageType >( curImage, saveFeatureImages,
          ".f%02d-1a.mha", fCount );
        }
      ++fCount;
      }
    for( unsigned int i=0; i<secondScales.size(); i++ )
      {
      int symmetryBase = fCount;
      for( unsigned int d=0; d<dimensionT; d++ )
        {
        typename LDAImageType::Pointer curImage = reader->GetOutput();
        for( unsigned int d2=0; d2<dimensionT; d2++ )
          {
          typedef itk::RecursiveGaussianImageFilter< LDAImageType,
            LDAImageType > FType;
          typename FType::Pointer filter = FType::New();
          filter->SetInput( curImage );
          filter->SetDirection( d2 );
          if( d == d2 )
            {
            filter->SetSigma( secondScales[i] );
            filter->SetOrder( FType::SecondOrder );
            }
          else
            {
            filter->SetSigma( secondScales[i]/2.0 );
            filter->SetOrder( FType::ZeroOrder );
            }
          filter->Update();
          curImage = filter->GetOutput();
          }
        ldaGenerator->AddFeatureImage( curImage );
        sprintf(s, ".f%02d.d2-%02d-%01d", fCount, (int)secondScales[i], d );
        featureName.push_back( featureBaseName + std::string( s ) );
        featureSymmetry.push_back( symmetryBase );
        if( saveFeatureImages.size() > 0 )
          {
          WriteLDA< LDAImageType >( curImage, saveFeatureImages,
            ".f%02d-2.mha", fCount );
          }
        ++fCount;
        }
      typename LDAImageType::Pointer curImage = reader->GetOutput();
      for( unsigned int d=0; d<dimensionT; d++ )
        {
        typedef itk::RecursiveGaussianImageFilter< LDAImageType,
          LDAImageType > FType;
        typename FType::Pointer filter = FType::New();
        filter->SetInput( curImage );
        filter->SetDirection( d );
        filter->SetSigma( secondScales[i] );
        filter->SetOrder( FType::SecondOrder );
        filter->Update();
        curImage = filter->GetOutput();
        }
      ldaGenerator->AddFeatureImage( curImage );
      sprintf(s, ".f%02d.d2-%02d-a", fCount, (int)secondScales[i] );
      featureName.push_back( featureBaseName + std::string( s ) );
      featureSymmetry.push_back( 0 );
      if( saveFeatureImages.size() > 0 )
        {
        WriteLDA< LDAImageType >( curImage, saveFeatureImages,
          ".f%02d-2a.mha", fCount );
        }
      ++fCount;
      }
    }

  typename MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
  inMaskReader->SetFileName( labelmap.c_str() );
  inMaskReader->Update();
  ldaGenerator->SetLabelmap( inMaskReader->GetOutput() );

  timeCollector.Stop( "LoadData" );

  if( objectId.size() > 0 )
    {
    ldaGenerator->SetObjectId( objectId[0] );
    if( objectId.size() > 1 )
      {
      for( unsigned int o=1; o<objectId.size(); o++ )
        {
        ldaGenerator->AddObjectId( objectId[o] );
        }
      }
    }

  if( usePCA )
    {
    ldaGenerator->SetPerformPCA( true );
    }

  if( loadLDAInfo.size() > 0 )
    {
    timeCollector.Start( "LoadLDA" );

    itk::tube::MetaLDA ldaReader( loadLDAInfo.c_str() );
    ldaReader.Read();

    ldaGenerator->SetLDAValues( ldaReader.GetLDAValues() );
    ldaGenerator->SetLDAMatrix( ldaReader.GetLDAMatrix() );

    timeCollector.Stop( "LoadLDA" );
    }
  else
    {
    timeCollector.Start( "Update" );

    ldaGenerator->Update();

    timeCollector.Stop( "Update" );
    }

  unsigned int numLDA = ldaGenerator->GetNumberOfLDA();
  if( useNumberOfLDA>0 && useNumberOfLDA < (int)numLDA )
    {
    numLDA = useNumberOfLDA;
    }

  vnl_vector< double > fVal( fCount );
  fVal.fill( 0 );
  for( unsigned int i=0; i<numLDA; i++ )
    {
    for( unsigned int f=0; f<fCount; f++ )
      {
      fVal[f] += vnl_math_abs( ldaGenerator->GetLDAValue( f ) );
      }
    }

  if( forceSymmetry )
    {
    for( unsigned int i=0; i<numLDA; i++ )
      {
      typename LDAGeneratorType::LDAVectorType v;
      v = ldaGenerator->GetLDAVector( i );
      for( unsigned int f=0; f<fCount; f++ )
        {
        if( featureSymmetry[f] != 0 )
          {
          double fSumS = 0;
          unsigned int fStart = f;
          while( featureSymmetry[f] == fStart )
            {
            fSumS += v[f]*v[f];
            ++f;
            }
          double symVal = vcl_sqrt( fSumS / (f-fStart) );
          if( fSumS > 0 )
            {
            f = fStart;
            while( featureSymmetry[f] == fStart )
              {
              std::cout << "sym: f = " << f << " : "
                << v[f] << " -> " << symVal << std::endl;
              v[f] = symVal;
              ++f;
              }
            }
          --f;
          }
        }
      ldaGenerator->SetLDAVector( i, v );
      }
    }

  if( forceSign )
    {
    for( unsigned int i=0; i<numLDA; i++ )
      {
      typename LDAGeneratorType::LDAVectorType v;
      v = ldaGenerator->GetLDAVector( i );
      unsigned int numFeatures = zeroOrderFeatures.size();
      double sum = 0;
      for( unsigned int j=0; j<numFeatures; j++ )
        {
        int f = zeroOrderFeatures[j];
        sum += v[f];
        }
      if( sum < 0 )
        {
        std::cout << "Sign: eVect = " << i << std::endl;
        for( unsigned int f=0; f<fCount; f++ )
          {
          v[f] *= -1;
          }
        ldaGenerator->SetLDAVector( i, v );
        }
      }
    }

  for( unsigned int f=0; f<fCount; f++ )
    {
    std::cout << featureName[f] << " : " << fVal[f] << std::endl;
    }

  if( outputBase.size() > 0 )
    {
    timeCollector.Start( "SaveLDAImages" );

    ldaGenerator->UpdateLDAImages();

    for( unsigned int i=0; i<numLDA; i++ )
      {
      typename LDAImageWriterType::Pointer ldaImageWriter =
        LDAImageWriterType::New();
      std::string fname = outputBase;
      char c[80];
      sprintf( c, ".lda%02d.mha", i );
      fname += std::string( c );
      ldaImageWriter->SetFileName( fname.c_str() );
      ldaImageWriter->SetInput( ldaGenerator->GetLDAImage( i ) );
      ldaImageWriter->Update();
      }
    timeCollector.Stop( "SaveLDAImages" );
    }

  if( saveLDAInfo.size() > 0 )
    {
    timeCollector.Start( "SaveLDA" );
    itk::tube::MetaLDA ldaWriter( ldaGenerator->GetLDAValues(),
      ldaGenerator->GetLDAMatrix() );
    ldaWriter.Write( saveLDAInfo.c_str() );
    timeCollector.Stop( "SaveLDA" );
    }

  timeCollector.Report();

  return 0;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolumesList[0], argc, argv );
}
