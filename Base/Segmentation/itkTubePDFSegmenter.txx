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

#ifndef __itkTubePDFSegmenter_txx
#define __itkTubePDFSegmenter_txx

#include <limits>

#include <vnl/vnl_matrix.h>

#include "itkTubePDFSegmenter.h"

#include "itkVectorImageToListGenerator.h"

#include <itkTimeProbesCollectorBase.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkHistogram.h>
#include <itkJoinImageFilter.h>
#include <itkHistogramToProbabilityImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>


//#include "itkPluginFilterWatcher.h"


namespace itk
{

namespace tube
{

template< class ImageT, unsigned int N, class LabelmapT >
PDFSegmenter< ImageT, N, LabelmapT >
::PDFSegmenter( void )
{
  m_SampleUpToDate = false;
  m_PDFsUpToDate = false;
  m_ImagesUpToDate = false;

  m_InClassList.clear();
  m_OutClassList  = NULL;

  m_InClassHistogram.clear();
  m_OutHistogram = NULL;
  m_HistogramBinMin.Fill( 0 );
  m_HistogramBinMax.Fill( 0 );
  m_HistogramBinScale.Fill( 0 );
  m_HistogramNumBinsND = 100;
  m_HistogramNumBins1D = 200;

  m_InputVolumeList.clear();
  m_InputVolumeList.resize( N, NULL );

  m_Labelmap = NULL;

  m_ObjectIdList.clear();
  m_VoidId = std::numeric_limits< MaskPixelType >::max();
  m_PDFWeightList.clear();

  m_ErodeRadius = 1;
  m_HoleFillIterations = 1;
  m_Draft = false;
  m_ProbabilityImageSmoothingStandardDeviation = 1;
  m_HistogramSmoothingStandardDeviation = 1;
  m_OutlierRejectPortion = 0.02;
  m_ReclassifyObjectMask = false;
  m_ReclassifyNotObjectMask = false;
  m_ForceClassification = false;

  m_ProbabilityImageVector.resize( 0 );

  m_ProgressProcessInfo = NULL;
  m_ProgressFraction = 1.0;
  m_ProgressStart = 0;
}

template< class ImageT, unsigned int N, class LabelmapT >
PDFSegmenter< ImageT, N, LabelmapT >
::~PDFSegmenter( void )
{
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetProgressProcessInformation( void * processInfo, double fraction,
  double start )
{
  m_ProgressProcessInfo = processInfo;
  m_ProgressFraction = fraction;
  m_ProgressStart = start;
}

template< class ImageT, unsigned int N, class LabelmapT >
int
PDFSegmenter< ImageT, N, LabelmapT >
::GetNumberOfClasses( void )
{
  return m_ProbabilityImageVector.size();
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::ClearObjectIds( void )
{
  m_ObjectIdList.clear();
  m_PDFWeightList.clear();
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
  m_PDFWeightList.push_back( 1.0 );
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::AddObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.push_back( objectId );
  m_PDFWeightList.push_back( 1.0 );
}

template< class ImageT, unsigned int N, class LabelmapT >
int
PDFSegmenter< ImageT, N, LabelmapT >
::GetObjectId( unsigned int classNum ) const
{
  if( classNum < m_ObjectIdList.size() )
    {
    return m_ObjectIdList[classNum];
    }

  if( classNum == m_ObjectIdList.size() )
    {
    return m_VoidId;
    }

  return -1;
   
}

template< class ImageT, unsigned int N, class LabelmapT >
unsigned int
PDFSegmenter< ImageT, N, LabelmapT >
::GetObjectNumberFromId( ObjectIdType id ) const
{
  for( unsigned int i=0; i<m_ObjectIdList.size(); i++ )
    {
    if( m_ObjectIdList[i] == id )
      {
      return i;
      }
    }
  throw;
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetObjectPDFWeight( unsigned int num, double weight )
{
  m_PDFWeightList[num] = weight;
}

template< class ImageT, unsigned int N, class LabelmapT >
double
PDFSegmenter< ImageT, N, LabelmapT >
::GetObjectPDFWeight( unsigned int num ) const
{
  return m_PDFWeightList[num];
}

template< class ImageT, unsigned int N, class LabelmapT >
const typename Image< float,
  ImageT::ImageDimension >::Pointer
PDFSegmenter< ImageT, N, LabelmapT >
::GetClassProbabilityVolume( unsigned int classNum ) const
{
  if( classNum < m_ProbabilityImageVector.size() )
    {
    return m_ProbabilityImageVector[classNum];
    }
  return NULL;
}

template< class ImageT, unsigned int N, class LabelmapT >
const typename PDFSegmenter< ImageT, N, LabelmapT >::PDFImageType
::Pointer
PDFSegmenter< ImageT, N, LabelmapT >
::GetClassPDFImage( unsigned int classNum )
{
  if( classNum < m_InClassHistogram.size() )
    {
    return m_InClassHistogram[classNum];
    }
  else if( classNum == m_InClassHistogram.size() )
    {
    return m_OutHistogram;
    }
  else
    {
    return NULL;
    }
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetClassPDFImage( unsigned int classNum,
  typename PDFImageType::Pointer classPDF )
{
  if( classNum < m_InClassHistogram.size() )
    {
    m_InClassHistogram[classNum] = classPDF;
    }
  else if( classNum == m_InClassHistogram.size() )
    {
    m_OutHistogram = classPDF;
    }
  m_SampleUpToDate = false;
  m_PDFsUpToDate = true;
  m_ImagesUpToDate = false;
}

template< class ImageT, unsigned int N, class LabelmapT >
double
PDFSegmenter< ImageT, N, LabelmapT >
::GetPDFBinMin( unsigned int featureNum )
{
  if( featureNum < N )
    {
    return m_HistogramBinMin[featureNum];
    }
  else
    {
    return 0;
    }
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetPDFBinMin( unsigned int featureNum, double val )
{
  if( featureNum < N )
    {
    m_HistogramBinMin[featureNum] = val;
    }
}

template< class ImageT, unsigned int N, class LabelmapT >
double
PDFSegmenter< ImageT, N, LabelmapT >
::GetPDFBinScale( unsigned int featureNum )
{
  if( featureNum < N )
    {
    return m_HistogramBinScale[featureNum];
    }
  else
    {
    return 0;
    }
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetInputVolume( unsigned int featureNum,
  typename ImageType::Pointer vol )
{
  if( featureNum < N )
    {
    m_InputVolumeList[featureNum] = vol;
    }
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetPDFBinScale( unsigned int featureNum, double val )
{
  if( featureNum < N )
    {
    m_HistogramBinScale[featureNum] = val;
    }
}


template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::GenerateSample( void )
{
  m_SampleUpToDate = true;

  unsigned int numClasses = m_ObjectIdList.size();

  //
  //  Convert in/out images to statistical list using masks
  //
  m_InClassList.resize( numClasses );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    m_InClassList[c] = ListSampleType::New();
    }
  m_OutClassList  = ListSampleType::New();

  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateSample" );

  typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType >
    ConstMaskImageIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex< ImageType >
    ConstImageIteratorType;

  ConstMaskImageIteratorType itInMask( m_Labelmap,
    m_Labelmap->GetLargestPossibleRegion() );
  itInMask.GoToBegin();

  ConstImageIteratorType * itInIm[N];
  for( unsigned int i=0; i<N; i++ )
    {
    itInIm[i] = new ConstImageIteratorType( m_InputVolumeList[i],
      m_InputVolumeList[i]->GetLargestPossibleRegion() );
    itInIm[i]->GoToBegin();
    m_HistogramBinMin[i] = itInIm[i]->Get();
    m_HistogramBinMax[i] = itInIm[i]->Get();
    }
  ListVectorType v;
  typename MaskImageType::IndexType indx;
  while( !itInMask.IsAtEnd() )
    {
    int val = itInMask.Get();
    indx = itInMask.GetIndex();
    for( unsigned int i=0; i<N; i++ )
      {
      v[i] = static_cast< PixelType >( itInIm[i]->Get() );
      }
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      v[N+i] = indx[i];
      }
    bool found = false;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      if( val == m_ObjectIdList[c] )
        {
        found = true;
        m_InClassList[c]->PushBack( v );
        break;
        }
      }
    if( !found && val != m_VoidId )
      {
      m_OutClassList->PushBack( v );
      }
    for( unsigned int i=0; i<N; i++ )
      {
      if( v[i]<m_HistogramBinMin[i] )
        {
        m_HistogramBinMin[i] = v[i];
        }
      else if( v[i]>m_HistogramBinMax[i] )
        {
        m_HistogramBinMax[i] = v[i];
        }
      }
    ++itInMask;
    for( unsigned int i=0; i<N; i++ )
      {
      ++( *( itInIm[i] ) );
      }
    if( m_Draft )
      {
      for( unsigned int count=1; count<4; count++ )
        {
        ++itInMask;
        for( unsigned int i=0; i<N; i++ )
          {
          ++( *( itInIm[i] ) );
          }
        }
      }
    }

  for( unsigned int i=0; i<N; i++ )
    {
    m_HistogramBinScale[i] = ( double )m_HistogramNumBinsND /
      ( m_HistogramBinMax[i] - m_HistogramBinMin[i] );
    }

  for( unsigned int i=0; i<N; i++ )
    {
    std::cout << "  Image [" << i << "] : Intensity = "
      << m_HistogramBinMin[i] << " - " << m_HistogramBinMax[i]
      << std::endl;
    }

  for( unsigned int i=0; i<m_ObjectIdList.size(); i++ )
    {
    std::cout << "ObjectId[" << i << "] = " << m_ObjectIdList[i]
      << std::endl;
    std::cout << "  size = " << m_InClassList[i]->Size() << std::endl;
    }
  std::cout << "VoidId = " << m_VoidId << std::endl;
  std::cout << "NotObjectId size = " << m_OutClassList->Size()
    << std::endl;
  for( unsigned int i=0; i<N; i++ )
    {
    delete itInIm[i];
    }

  timeCollector.Stop( "GenerateSample" );

  timeCollector.Report();
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::GeneratePDFs( void )
{
  if( !m_SampleUpToDate )
    {
    this->GenerateSample();
    }
  m_PDFsUpToDate = true;

  itk::TimeProbesCollectorBase timeCollector;

  unsigned int numClasses = m_ObjectIdList.size();

  //
  // Convert lists to histograms that have the same span ( using range
  //   defined above )
  //

  timeCollector.Start( "ListsToHistograms" );
  // Inside

  std::vector< std::vector< float > > clipMin;
  std::vector< std::vector< float > > clipMax;
  if( true ) // creating a local context to limit memory footprint
    {
    clipMin.resize( numClasses );
    clipMax.resize( numClasses );

    std::vector< vnl_matrix< float > > inImHistogram;
    inImHistogram.resize( numClasses );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      clipMin[c].resize( N );
      clipMax[c].resize( N );
      inImHistogram[c].set_size( N, m_HistogramNumBinsND );
      inImHistogram[c].fill( 0 );
      }

    unsigned int totalIn = 0;
    itk::Array<unsigned int> totalInClass( numClasses );
    totalInClass.Fill( 0 );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      typename ListSampleType::ConstIterator
        inClassListIt( m_InClassList[c]->Begin() );
      typename ListSampleType::ConstIterator
        inClassListItEnd( m_InClassList[c]->End() );
      double binV;
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i=0; i<N; i++ )
          {
          binV = inClassListIt.GetMeasurementVector()[i];
          binV = ( int )( ( binV - m_HistogramBinMin[i] )
            * m_HistogramBinScale[i] + 0.5 );
          if( binV>m_HistogramNumBinsND-1 )
            {
            binV = m_HistogramNumBinsND-1;
            }
          else if( binV<0 )
            {
            binV = 0;
            }
          ++( ( inImHistogram[c][i] )[( int )binV] );
          }
        ++totalIn;
        ++totalInClass[c];
        ++inClassListIt;
        }
      }

    unsigned int count;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      double tailReject = totalInClass[c] * ( m_OutlierRejectPortion/2 );
      for( unsigned int i=0; i<N; i++ )
        {
        count = 0;
        for( unsigned int b=0; b<m_HistogramNumBinsND; b++ )
          {
          count += static_cast<unsigned int>( inImHistogram[c][i][b] );
          if( count>=tailReject )
            {
            clipMin[c][i] =  b / m_HistogramBinScale[i]
              + m_HistogramBinMin[i];
            break;
            }
          }
        count = 0;
        for( int b=( int )m_HistogramNumBinsND-1; b>=0; b-- )
          {
          count += static_cast<unsigned int>( inImHistogram[c][i][b] );
          if( count>=tailReject )
            {
            clipMax[c][i] =  b / m_HistogramBinScale[i]
              + m_HistogramBinMin[i];
            break;
            }
          }
        }
      std::cout << "clipMin[" << c << "]=";
      for( unsigned int i=0; i<N; i++ )
        {
        std::cout << clipMin[c][i] << " ";
        }
      std::cout << std::endl;
      std::cout << "clipMax[" << c << "]=";
      for( unsigned int i=0; i<N; i++ )
        {
        std::cout << clipMax[c][i] << " ";
        }
      std::cout << std::endl;
      }
    }
  timeCollector.Stop( "ListsToHistograms" );

  //
  //  Create joint histograms
  //
  timeCollector.Start( "JointClassHistogram" );

  typedef itk::ImageRegionIteratorWithIndex< HistogramImageType >
    HistogramIteratorType;
  typename HistogramImageType::SizeType size;
  size.Fill( m_HistogramNumBinsND );

  m_InClassHistogram.resize( numClasses );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    m_InClassHistogram[c] = HistogramImageType::New();
    m_InClassHistogram[c]->SetRegions( size );
    m_InClassHistogram[c]->Allocate();
    m_InClassHistogram[c]->FillBuffer( 0 );

    unsigned int count = 0;
    typename ListSampleType::ConstIterator
      inClassListIt( m_InClassList[c]->Begin() );
    typename ListSampleType::ConstIterator
      inClassListItEnd( m_InClassList[c]->End() );
    typename HistogramImageType::IndexType indxHistogram;
    while( inClassListIt != inClassListItEnd )
      {
      bool valid = true;
      for( unsigned int i=0; i<N; i++ )
        {
        double binV = inClassListIt.GetMeasurementVector()[i];
        if( binV >= clipMin[c][i] && binV <= clipMax[c][i] )
          {
          binV = ( int )( ( binV - m_HistogramBinMin[i] )
            * m_HistogramBinScale[i] + 0.5 );
          if( binV<0 || binV>=m_HistogramNumBinsND )
            {
            valid = false;
            break;
            }
          indxHistogram[i] = ( int )binV;
          }
        else
          {
          valid = false;
          break;
          }
        }
      if( valid )
        {
        ++count;
        ++( m_InClassHistogram[c]->GetPixel( indxHistogram ) );
        }
      ++inClassListIt;
      }
    std::cout << "m_InClassHistogram size = " << count << std::endl;
    }
  timeCollector.Stop( "JointClassHistogram" );

  //
  //  Convert histograms to images so that we can perform blurring to
  //    generate parzen window density estimates
  //
  timeCollector.Start( "HistogramToPDF" );
  if( true ) // creating a local context to limit memory footprint
    {
    std::cout << "Inside histogram image smoothing..." << std::endl;
    typedef itk::DiscreteGaussianImageFilter< HistogramImageType,
      HistogramImageType > HistogramBlurGenType;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      double inPTotal = 0;

      typename HistogramBlurGenType::Pointer m_InClassHistogramBlurGen =
        HistogramBlurGenType::New();
      m_InClassHistogramBlurGen->SetInput( m_InClassHistogram[c] );
      m_InClassHistogramBlurGen->SetVariance(
        m_HistogramSmoothingStandardDeviation *
        m_HistogramSmoothingStandardDeviation );
      m_InClassHistogramBlurGen->Update();
      m_InClassHistogram[c] = m_InClassHistogramBlurGen->GetOutput();

      itk::ImageRegionIterator<HistogramImageType> m_InClassHistogramIt(
        m_InClassHistogram[c],
        m_InClassHistogram[c]->GetLargestPossibleRegion() );
      m_InClassHistogramIt.GoToBegin();
      while( !m_InClassHistogramIt.IsAtEnd() )
        {
        double tf = m_InClassHistogramIt.Get();
        inPTotal += tf;
        ++m_InClassHistogramIt;
        }
      std::cout << "m_InClassHistogramTotalP = " << inPTotal << std::endl;

      if( inPTotal > 0 )
        {
        m_InClassHistogramIt.GoToBegin();
        while( !m_InClassHistogramIt.IsAtEnd() )
          {
          double tf = m_InClassHistogramIt.Get();
          m_InClassHistogramIt.Set( tf / inPTotal );
          ++m_InClassHistogramIt;
          }
        }
      }
    }
  timeCollector.Stop( "HistogramToPDF" );

  timeCollector.Report();
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::ApplyPDFs( void )
{
  if( !m_SampleUpToDate )
    {
    this->GenerateSample();
    }
  if( !m_PDFsUpToDate )
    {
    this->GeneratePDFs();
    }
  m_ImagesUpToDate = true;

  int erodeRadius = m_ErodeRadius;
  int holeFillIterations = m_HoleFillIterations;
  if( m_Draft )
    {
    erodeRadius /= 2;
    holeFillIterations /= 2;
    }

  unsigned int numClasses = m_ObjectIdList.size();

  //
  //  Compute the probability at each pixel for input images
  //
  std::cout << "Compute probability image..." << std::endl;

  m_ProbabilityImageVector.resize( numClasses + 1 );

  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "ProbabilityImage" );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    m_ProbabilityImageVector[c] = ProbabilityImageType::New();
    m_ProbabilityImageVector[c]->SetRegions(
      m_InputVolumeList[0]->GetLargestPossibleRegion() );
    m_ProbabilityImageVector[c]->CopyInformation( m_InputVolumeList[0] );
    m_ProbabilityImageVector[c]->Allocate();

    itk::ImageRegionIterator<ProbabilityImageType> probIt(
      m_ProbabilityImageVector[c],
      m_ProbabilityImageVector[c]->GetLargestPossibleRegion() );
    probIt.GoToBegin();

    typedef itk::ImageRegionConstIteratorWithIndex< ImageType >
      ConstImageIteratorType;
    std::vector< ConstImageIteratorType * > itInIm(N);
    for( unsigned int i=0; i<N; i++ )
      {
      itInIm[i] = new ConstImageIteratorType( m_InputVolumeList[i],
        m_InputVolumeList[i]->GetLargestPossibleRegion() );
      itInIm[i]->GoToBegin();
      }

    typename HistogramImageType::IndexType binIndex;
    while( !probIt.IsAtEnd() )
      {
      bool valid = true;
      for( unsigned int i=0; i<N; i++ )
        {
        double binV = itInIm[i]->Get();
        binV = ( int )( ( binV - m_HistogramBinMin[i] )
          * m_HistogramBinScale[i] + 0.5 );
        if( binV<0 || binV>m_HistogramNumBinsND-1 )
          {
          valid = false;
          break;
          }
        binIndex[i] = static_cast<long>( binV );
        }
      double prob = 0;
      if( valid )
        {
        if( c < numClasses )
          {
          prob = m_InClassHistogram[c]->GetPixel( binIndex );
          }
        else
          {
          prob = m_OutHistogram->GetPixel( binIndex );
          }
        }
      prob *= m_PDFWeightList[c];
      probIt.Set( prob );
      for( unsigned int i=0; i<N; i++ )
        {
        ++( *( itInIm[i] ) );
        }
      ++probIt;
      }

    for( unsigned int i=0; i<N; i++ )
      {
      delete itInIm[i];
      }
    }
  timeCollector.Stop( "ProbabilityImage" );

  std::cout << "Probability image smoothing..." << std::endl;
  typedef itk::DiscreteGaussianImageFilter< ProbabilityImageType,
          ProbabilityImageType> ProbImageFilterType;
  typename ProbImageFilterType::Pointer probImageFilter;

  timeCollector.Start( "ProbabilityImageDiffusion" );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    probImageFilter = ProbImageFilterType::New();
    probImageFilter->SetInput( m_ProbabilityImageVector[c] );
    probImageFilter->SetVariance(
      m_ProbabilityImageSmoothingStandardDeviation *
      m_ProbabilityImageSmoothingStandardDeviation );
    probImageFilter->Update();
    m_ProbabilityImageVector[c] = probImageFilter->GetOutput();
    }
  timeCollector.Stop( "ProbabilityImageDiffusion" );

  std::cout << "Inside connectivity..." << std::endl;

  //
  //  Create label image
  //

  typename MaskImageType::IndexType labelImageIndex;

  typename MaskImageType::Pointer tmpLabelImage = MaskImageType::New();
  tmpLabelImage->SetRegions( m_InputVolumeList[0]
    ->GetLargestPossibleRegion() );
  tmpLabelImage->CopyInformation( m_InputVolumeList[0] );
  tmpLabelImage->Allocate();

  if( !m_ForceClassification )
    {
    for( unsigned int c=0; c<numClasses; c++ )
      {
      timeCollector.Start( "Connectivity" );

      // For this class, label all pixels for which it is the most
      // likely class.
      itk::ImageRegionIteratorWithIndex<MaskImageType> labelIt(
        tmpLabelImage, tmpLabelImage->GetLargestPossibleRegion() );
      labelIt.GoToBegin();
      while( !labelIt.IsAtEnd() )
        {
        labelImageIndex = labelIt.GetIndex();
        bool maxPC = true;
        double maxP = m_ProbabilityImageVector[c]->GetPixel(
          labelImageIndex );
        for( unsigned int oc=0; oc<numClasses; oc++ )
          {
          if( oc != c &&
              m_ProbabilityImageVector[oc]->GetPixel( labelImageIndex )
              > maxP )
            {
            maxPC = false;
            break;
            }
          }
        if( maxPC )
          {
          labelIt.Set( 128 );
          }
        else
          {
          labelIt.Set( 0 );
          }
        ++labelIt;
        }

      typedef itk::ConnectedThresholdImageFilter<MaskImageType,
        MaskImageType> ConnectedFilterType;
      typename ConnectedFilterType::Pointer insideConnecter =
        ConnectedFilterType::New();

      typename ListSampleType::ConstIterator
        inClassListIt( m_InClassList[c]->Begin() );
      typename ListSampleType::ConstIterator
        inClassListItEnd( m_InClassList[c]->End() );
      typename ImageType::IndexType indx;
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          indx[i] = static_cast<long int>(
            inClassListIt.GetMeasurementVector()[N+i] );
          }
        insideConnecter->AddSeed( indx );
        ++inClassListIt;
        }

      if( !m_ReclassifyObjectMask )
        {
        // The pixels with maximum probability for the current
        // class are all set to 128 before the update of the
        // ConnectedThresholdImageFilter, so if the input label
        // map is not to be reclassified, set the pixels belonging
        // to this class to 128 before updating the filter regardless
        // of the probability.  Setting input labels to 255 before
        // the update
        // of the ConnectedThresholdFilter will cause the filter to
        // return only the values at 255 (the input label map).
        inClassListIt = m_InClassList[c]->Begin();
        inClassListItEnd = m_InClassList[c]->End();
        while( inClassListIt != inClassListItEnd )
          {
          for( unsigned int i=0; i<ImageDimension; i++ )
            {
            indx[i] = static_cast<long int>(
              inClassListIt.GetMeasurementVector()[N+i] );
            }
          tmpLabelImage->SetPixel( indx, 128 );
          ++inClassListIt;
          }
        for( unsigned int oc=0; oc<numClasses; oc++ )
          {
          if( oc != c )
            {
            // Use inside mask to set seed points.  Also draw inside mask
            // in
            // label image to ensure those points are considered object
            // points.
            // Erase other mask from label image
            inClassListIt = m_InClassList[oc]->Begin();
            inClassListItEnd = m_InClassList[oc]->End();
            while( inClassListIt != inClassListItEnd )
              {
              for( unsigned int i=0; i<ImageDimension; i++ )
                {
                indx[i] = static_cast<long int>(
                  inClassListIt.GetMeasurementVector()[N+i] );
                }
              tmpLabelImage->SetPixel( indx, 0 );
              ++inClassListIt;
              }
            }
          }

        // Erase outside mask from label image
        typename ListSampleType::ConstIterator
          outListIt( m_OutClassList->Begin() );
        typename ListSampleType::ConstIterator
          outListItEnd( m_OutClassList->End() );
        while( outListIt != outListItEnd )
          {
          for( unsigned int i=0; i<ImageDimension; i++ )
            {
            indx[i] = static_cast<long int>(
              outListIt.GetMeasurementVector()[N+i] );
            }
          tmpLabelImage->SetPixel( indx, 0 );
          ++outListIt;
          }
        }

      insideConnecter->SetInput( tmpLabelImage );
      insideConnecter->SetLower( 64 );
      insideConnecter->SetUpper( 194 );
      insideConnecter->SetReplaceValue( 255 );
      insideConnecter->Update();
      tmpLabelImage = insideConnecter->GetOutput();

      timeCollector.Stop( "Connectivity" );

      //
      // Fill holes
      //
      if( holeFillIterations > 0 )
        {
        typedef itk::VotingBinaryIterativeHoleFillingImageFilter<
          MaskImageType > HoleFillingFilterType;

        std::cout << "Fill holes..." << std::endl;

        timeCollector.Start( "HoleFiller" );

        typename HoleFillingFilterType::Pointer holeFiller =
          HoleFillingFilterType::New();
        typename MaskImageType::SizeType holeRadius;
        holeRadius.Fill( 1 );
        holeFiller->SetInput( tmpLabelImage );
        holeFiller->SetRadius( holeRadius );
        holeFiller->SetBackgroundValue( 0 );
        holeFiller->SetForegroundValue( 255 );
        holeFiller->SetMajorityThreshold( 2 );
        holeFiller->SetMaximumNumberOfIterations( holeFillIterations );
        holeFiller->Update();
        tmpLabelImage = holeFiller->GetOutput();

        timeCollector.Stop( "HoleFiller" );
        }

      //
      // Erode
      //
      typedef BinaryBallStructuringElement< MaskPixelType,
        ImageType::ImageDimension >
          StructuringElementType;
      typedef BinaryErodeImageFilter< MaskImageType, MaskImageType,
        StructuringElementType >
          ErodeFilterType;

      StructuringElementType sphereOp;
      if( erodeRadius > 0 )
        {
        std::cout << "Inside erode..." << std::endl;

        timeCollector.Start( "Erode" );

        sphereOp.SetRadius( erodeRadius );
        sphereOp.CreateStructuringElement();

        typename ErodeFilterType::Pointer insideMaskErodeFilter =
          ErodeFilterType::New();
        insideMaskErodeFilter->SetKernel( sphereOp );
        insideMaskErodeFilter->SetErodeValue( 255 );
        insideMaskErodeFilter->SetInput( tmpLabelImage );
        insideMaskErodeFilter->Update();
        tmpLabelImage = insideMaskErodeFilter->GetOutput();

        timeCollector.Stop( "Erode" );
        }

      //
      // Re-do connectivity
      //
      if( true ) // creating a local context to limit memory footprint
        {
        typedef itk::ConnectedThresholdImageFilter<MaskImageType,
          MaskImageType> ConnectedMaskFilterType;

        std::cout << "Inside connectivity pass 2..." << std::endl;

        timeCollector.Start( "Connectivity2" );

        typename ConnectedMaskFilterType::Pointer
          insideConnectedMaskFilter = ConnectedMaskFilterType::New();
        insideConnectedMaskFilter->SetInput( tmpLabelImage );
        insideConnectedMaskFilter->SetLower( 194 );
        insideConnectedMaskFilter->SetUpper( 255 );
        insideConnectedMaskFilter->SetReplaceValue( 255 );

        // Use inside mask to set seed points.  Also draw inside mask in
        // label image to ensure those points are considered object points
        inClassListIt = m_InClassList[c]->Begin();
        inClassListItEnd = m_InClassList[c]->End();
        if( !m_ReclassifyObjectMask )
          {
          while( inClassListIt != inClassListItEnd )
            {
            for( unsigned int i=0; i<ImageDimension; i++ )
              {
              indx[i] = static_cast<long int>(
                inClassListIt.GetMeasurementVector()[N+i] );
              }
            insideConnectedMaskFilter->AddSeed( indx );
            tmpLabelImage->SetPixel( indx, 255 );  // Redraw objects
            ++inClassListIt;
            }
          }
        else
          {
          while( inClassListIt != inClassListItEnd )
            {
            for( unsigned int i=0; i<ImageDimension; i++ )
              {
              indx[i] = static_cast<long int>(
                inClassListIt.GetMeasurementVector()[N+i] );
              }

            insideConnectedMaskFilter->AddSeed( indx );
            // Don't redraw objects
            ++inClassListIt;
            }
          }

        insideConnectedMaskFilter->Update();
        tmpLabelImage = insideConnectedMaskFilter->GetOutput();

        timeCollector.Stop( "Connectivity2" );
        }

      //
      // Dilate back to original size
      //
      typedef itk::BinaryDilateImageFilter< MaskImageType,
        MaskImageType, StructuringElementType >           DilateFilterType;

      if( erodeRadius > 0 )
        {
        std::cout << "Inside eroded-connected dilate..." << std::endl;

        timeCollector.Start( "Dilate" );

        typename DilateFilterType::Pointer insideMaskDilateFilter =
          DilateFilterType::New();
        insideMaskDilateFilter->SetKernel( sphereOp );
        insideMaskDilateFilter->SetDilateValue( 255 );
        insideMaskDilateFilter->SetInput( tmpLabelImage );
        insideMaskDilateFilter->Update();
        tmpLabelImage = insideMaskDilateFilter->GetOutput();

        timeCollector.Stop( "Dilate" );
        }

      // Merge with input mask
      typedef itk::ImageRegionIterator< MaskImageType >
        MaskImageIteratorType;
      MaskImageIteratorType itInMask( m_Labelmap,
        m_Labelmap->GetLargestPossibleRegion() );
      itInMask.GoToBegin();

      MaskImageIteratorType itLabel( tmpLabelImage,
        tmpLabelImage->GetLargestPossibleRegion() );
      itLabel.GoToBegin();

      while( !itInMask.IsAtEnd() )
        {
        if( itLabel.Get() == 255 )
          {
          if( itInMask.Get() == m_VoidId
              || ( m_ReclassifyObjectMask && m_ReclassifyNotObjectMask ) )
            {
            itInMask.Set( m_ObjectIdList[c] );
            }
          else
            {
            if( m_ReclassifyObjectMask || m_ReclassifyNotObjectMask )
              {
              bool isObjectId = false;
              for( unsigned int oc=0; oc<numClasses; oc++ )
                {
                if( itInMask.Get() == m_ObjectIdList[oc] )
                  {
                  isObjectId = true;
                  break;
                  }
                }
              if( ( isObjectId && m_ReclassifyObjectMask ) ||
                  ( !isObjectId && m_ReclassifyNotObjectMask ) )
                {
                itInMask.Set( m_ObjectIdList[c] );
                }
              }
            }
          }
        else
          {
          if( itInMask.Get() == m_ObjectIdList[c] )
            {
            if( m_ReclassifyObjectMask )
              {
              itInMask.Set( m_VoidId );
              }
            }
          }
        ++itInMask;
        ++itLabel;
        }
      }
    }
  else
    {
    timeCollector.Start( "ForceClassification" );

    // Merge with input mask
    typedef itk::ImageRegionIteratorWithIndex< MaskImageType >
      MaskImageIteratorType;
    MaskImageIteratorType itInMask( m_Labelmap,
      m_Labelmap->GetLargestPossibleRegion() );
    itInMask.GoToBegin();

    while( !itInMask.IsAtEnd() )
      {
      labelImageIndex = itInMask.GetIndex();
      unsigned int maxPC = 0;
      double maxP = m_ProbabilityImageVector[0]->GetPixel(
        labelImageIndex );
      for( unsigned int c=1; c<numClasses; c++ )
        {
        double p = m_ProbabilityImageVector[c]->GetPixel(
          labelImageIndex );
        if( p > maxP )
          {
          maxP = p;
          maxPC = c;
          }
        }
      if( itInMask.Get() == m_VoidId
        || ( m_ReclassifyObjectMask && m_ReclassifyNotObjectMask ) )
        {
        itInMask.Set( m_ObjectIdList[maxPC] );
        }
      else
        {
        if( m_ReclassifyObjectMask || m_ReclassifyNotObjectMask )
          {
          bool isObjectId = false;
          for( unsigned int oc=0; oc<numClasses; oc++ )
            {
            if( itInMask.Get() == m_ObjectIdList[oc] )
              {
              isObjectId = true;
              break;
              }
            }
          if( ( isObjectId && m_ReclassifyObjectMask ) ||
              ( !isObjectId && m_ReclassifyNotObjectMask ) )
            {
            itInMask.Set( m_ObjectIdList[maxPC] );
            }
          }
        }
      ++itInMask;
      }

    timeCollector.Stop( "ForceClassification" );
    }

  timeCollector.Report();
  m_ImagesUpToDate = true;
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::Update( void )
{
  this->GenerateSample();
  this->GeneratePDFs();
}

template< class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::ClassifyImages( void )
{
  this->ApplyPDFs();
}

template<class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_SampleUpToDate )
    {
    os << indent << "SampleUpToDate = true" << std::endl;
    }
  else
    {
    os << indent << "SampleUpToDate = false" << std::endl;
    }

  if( m_PDFsUpToDate )
    {
    os << indent << "PDFsUpToDate = true" << std::endl;
    }
  else
    {
    os << indent << "PDFsUpToDate = false" << std::endl;
    }

  if( m_ImagesUpToDate )
    {
    os << indent << "ImagesUpToDate = true" << std::endl;
    }
  else
    {
    os << indent << "ImagesUpToDate = false" << std::endl;
    }

  os << indent << "Input volume list size = " << m_InputVolumeList.size()
    << std::endl;

  if( m_Labelmap.IsNotNull() )
    {
    os << indent << "Labelmap = " << m_Labelmap << std::endl;
    }
  else
    {
    os << indent << "Labelmap = NULL" << std::endl;
    }
  os << indent << "Erode radius = " << m_ErodeRadius << std::endl;
  os << indent << "Hole fill iterations = " << m_HoleFillIterations
    << std::endl;
  os << indent << "PDF weight size = " << m_PDFWeightList.size()
    << std::endl;
  os << indent << "Probability Smoothing Standard Deviation = "
    << m_ProbabilityImageSmoothingStandardDeviation << std::endl;
  os << indent << "Histogram Smoothing Standard Deviation = "
    << m_HistogramSmoothingStandardDeviation << std::endl;
  os << indent << "Outlier reject portion = "
    << m_OutlierRejectPortion << std::endl;
  os << indent << "Draft = " << m_Draft << std::endl;
  os << indent << "ReclassifyObjectMask = " << m_ReclassifyObjectMask
    << std::endl;
  os << indent << "ReclassifyNotObjectMask = " << m_ReclassifyNotObjectMask
    << std::endl;
  os << indent << "Number of probability images = "
    << m_ProbabilityImageVector.size() << std::endl;
  os << indent << "InClassList size = "
    << m_InClassList.size() << std::endl;
  if( m_OutClassList.IsNotNull() )
    {
    os << indent << "OutClassList = " << m_OutClassList << std::endl;
    }
  else
    {
    os << indent << "OutClassList = NULL" << std::endl;
    }
  os << indent << "InClassHistogram size = "
    << m_InClassHistogram.size() << std::endl;
  os << indent << "HistogramBinMin = " << m_HistogramBinMin << std::endl;
  os << indent << "HistogramBinMax = " << m_HistogramBinMax << std::endl;
  os << indent << "HistogramBinScale = " << m_HistogramBinScale
    << std::endl;
  os << indent << "HistogramNumBinsND = "
    << m_HistogramNumBinsND << std::endl;
  os << indent << "HistogramNumBins1D = "
    << m_HistogramNumBins1D << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkTubePDFSegmenter_txx)
