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
#ifndef __itkTubePDFSegmenter_txx
#define __itkTubePDFSegmenter_txx

#include <limits>

#include <vnl/vnl_matrix.h>

#include "itkTubePDFSegmenter.h"

#include "itkVectorImageToListGenerator.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHistogram.h"
#include "itkJoinImageFilter.h"
#include "itkHistogramToProbabilityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"


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
  m_OutList = NULL;

  m_InClassHisto.clear();
  m_OutHisto = NULL;
  m_HistoBinMin.Fill( 0 );
  m_HistoBinMax.Fill( 0 );
  m_HistoBinScale.Fill( 0 );
  m_HistoNumBinsND = 100;
  m_HistoNumBins1D = 200;

  m_InputVolumeList.clear();
  m_InputVolumeList.resize( N, NULL );

  m_Labelmap = NULL;

  m_ObjectIdList.clear();
  m_VoidId = std::numeric_limits< MaskPixelType >::max();

  m_ErodeRadius = 1;
  m_HoleFillIterations = 1;
  m_FprWeight = 1.5;
  m_Draft = false;
  m_ProbabilitySmoothingStandardDeviation = 0.5;
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

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetProgressProcessInformation( void * processInfo, double fraction,
  double start )
{
  m_ProgressProcessInfo = processInfo;
  m_ProgressFraction = fraction;
  m_ProgressStart = start;
}

template < class ImageT, unsigned int N, class LabelmapT >
int
PDFSegmenter< ImageT, N, LabelmapT >
::GetNumberOfClasses( void )
{
  return m_ProbabilityImageVector.size();
}

template < class ImageT, unsigned int N, class LabelmapT >
int
PDFSegmenter< ImageT, N, LabelmapT >
::GetClassId( unsigned int classNum )
{
  if( classNum < m_ObjectIdList.size() )
    {
    return m_ObjectIdList[classNum];
    }
  else if( classNum == m_ObjectIdList.size() )
    {
    return m_VoidId;
    }
  else
    {
    return -1;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
const typename Image< float,
  ImageT::ImageDimension >::Pointer
PDFSegmenter< ImageT, N, LabelmapT >
::GetClassProbabilityVolume( unsigned int classNum )
{
  if( classNum < m_ProbabilityImageVector.size() )
    {
    return m_ProbabilityImageVector[classNum];
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
const typename PDFSegmenter< ImageT, N, LabelmapT >::PDFImageType
::Pointer
PDFSegmenter< ImageT, N, LabelmapT >
::GetClassPDFImage( unsigned int classNum )
{
  if( classNum < m_InClassHisto.size() )
    {
    return m_InClassHisto[classNum];
    }
  else if( classNum == m_InClassHisto.size() )
    {
    return m_OutHisto;
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetClassPDFImage( unsigned int classNum,
  typename PDFImageType::Pointer classPDF )
{
  if( classNum < m_InClassHisto.size() )
    {
    m_InClassHisto[classNum] = classPDF;
    }
  else if( classNum == m_InClassHisto.size() )
    {
    m_OutHisto = classPDF;
    }
  m_SampleUpToDate = false;
  m_PDFsUpToDate = true;
  m_ImagesUpToDate = false;
}

template < class ImageT, unsigned int N, class LabelmapT >
double
PDFSegmenter< ImageT, N, LabelmapT >
::GetPDFBinMin( unsigned int featureNum )
{
  if( featureNum < N )
    {
    return m_HistoBinMin[featureNum];
    }
  else
    {
    return 0;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetPDFBinMin( unsigned int featureNum, double val )
{
  if( featureNum < N )
    {
    m_HistoBinMin[featureNum] = val;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
double
PDFSegmenter< ImageT, N, LabelmapT >
::GetPDFBinScale( unsigned int featureNum )
{
  if( featureNum < N )
    {
    return m_HistoBinScale[featureNum];
    }
  else
    {
    return 0;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
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

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::SetPDFBinScale( unsigned int featureNum, double val )
{
  if( featureNum < N )
    {
    m_HistoBinScale[featureNum] = val;
    }
}


template < class ImageT, unsigned int N, class LabelmapT >
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
  m_OutList = ListSampleType::New();

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
    }
  ListVectorType v;
  typename MaskImageType::IndexType indx;
  while( !itInMask.IsAtEnd() )
    {
    int val = itInMask.Get();
    indx = itInMask.GetIndex();
    bool found = false;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      if( val == m_ObjectIdList[c] )
        {
        found = true;
        for( unsigned int i=0; i<N; i++ )
          {
          v[i] = static_cast< PixelType >( itInIm[i]->Get() );
          }
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          v[N+i] = indx[i];
          }
        m_InClassList[c]->PushBack( v );
        break;
        }
      }
    if( !found && itInMask.Get() != m_VoidId )
      {
      for( unsigned int i=0; i<N; i++ )
        {
        v[i] = static_cast< PixelType >( itInIm[i]->Get() );
        }
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        v[N+i] = indx[i];
        }
      m_OutList->PushBack( v );
      }
    for( unsigned int count=0; ( m_Draft && count<4 ) || ( count<1 );
      count++ )
      {
      ++itInMask;
      for( unsigned int i=0; i<N; i++ )
        {
        ++( *( itInIm[i] ) );
        }
      }
    }

  for( unsigned int i=0; i<m_ObjectIdList.size(); i++ )
    {
    std::cout << "ObjectId[" << i << "] = " << m_ObjectIdList[i]
      << std::endl;
    std::cout << "  size = " << m_InClassList[i]->Size() << std::endl;
    }
  std::cout << "VoidId = " << m_VoidId << std::endl;
  std::cout << "  size = " << m_OutList->Size() << std::endl;
  for( unsigned int i=0; i<N; i++ )
    {
    delete itInIm[i];
    }

  timeCollector.Stop( "GenerateSample" );

  timeCollector.Report();
}

template < class ImageT, unsigned int N, class LabelmapT >
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

  timeCollector.Start( "HistogramMinMax" );

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
    m_HistoBinMin[i] = itInIm[i]->Get();
    m_HistoBinMax[i] = itInIm[i]->Get();
    }
  ListVectorType v;
  while( !itInMask.IsAtEnd() )
    {
    for( unsigned int i=0; i<N; i++ )
      {
      v[i] = static_cast< PixelType >( itInIm[i]->Get() );
      }
    for( unsigned int i=0; i<N; i++ )
      {
      if( v[i]<m_HistoBinMin[i] )
        {
        m_HistoBinMin[i] = v[i];
        }
      else if( v[i]>m_HistoBinMax[i] )
        {
        m_HistoBinMax[i] = v[i];
        }
      }
    for( unsigned int count=0; ( m_Draft && count<4 ) || ( count<1 );
      count++ )
      {
      ++itInMask;
      for( unsigned int i=0; i<N; i++ )
        {
        ++( *( itInIm[i] ) );
        }
      }
    }

  for( unsigned int i=0; i<N; i++ )
    {
    delete itInIm[i];
    }

  for( unsigned int i=0; i<N; i++ )
    {
    m_HistoBinScale[i] = ( double )m_HistoNumBinsND /
      ( m_HistoBinMax[i] - m_HistoBinMin[i] );
    }

  for( unsigned int i=0; i<N; i++ )
    {
    std::cout << "  Image" << i << ": " << m_HistoBinMin[i] << " - "
      << m_HistoBinMax[i] << std::endl;
    }
  timeCollector.Stop( "HistogramMinMax" );

  //
  // Convert lists to histograms that have the same span ( using range
  //   defined above )
  //

  timeCollector.Start( "ListsToHistograms" );
  // Inside

  std::vector< vnl_matrix< float > > inImHisto;
  if( true ) // creating a local context to limit memory footprint
    {
    inImHisto.resize( numClasses );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      inImHisto[c].set_size( N, m_HistoNumBinsND );
      inImHisto[c].fill( 0 );
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
          binV = ( int )( ( binV - m_HistoBinMin[i] )
            * m_HistoBinScale[i] + 0.5 );
          if( binV>m_HistoNumBinsND-1 )
            {
            binV = m_HistoNumBinsND-1;
            }
          else if( binV<0 )
            {
            binV = 0;
            }
          ++( ( inImHisto[c][i] )[( int )binV] );
          }
        ++totalIn;
        ++totalInClass[c];
        ++inClassListIt;
        }
      }

    double outlierRejectPercent = 0.01;
    unsigned int count;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      double totalReject = totalInClass[c] * outlierRejectPercent;
      for( unsigned int i=0; i<N; i++ )
        {
        count = 0;
        for( unsigned int b=0; b<m_HistoNumBinsND; b++ )
          {
          count += static_cast<unsigned int>( inImHisto[c][i][b] );
          if( count>=totalReject )
            {
            //m_HistoBinMin[i] = b/m_HistoBinScale[i] + m_HistoBinMin[i];
            inImHisto[c][i][b] = count - totalReject;
            break;
            }
          else
            {
            inImHisto[c][i][b] = 0;
            }
          }
        count = 0;
        for( int b=( int )m_HistoNumBinsND-1; b>=0; b-- )
          {
          count += static_cast<unsigned int>( inImHisto[c][i][b] );
          if( count>=totalReject )
            {
            //m_HistoBinMax[i] = b/m_HistoBinScale[i] + m_HistoBinMin[i];
            inImHisto[c][i][b] = count - totalReject;
            break;
            }
          else
            {
            inImHisto[c][i][b] = 0;
            }
          }
        }
      }
    }

  // Outside
  std::vector< std::vector< float > > outImHisto(N);
  if( true ) // creating a local context to limit memory footprint
    {
    for( unsigned int i=0; i<N; i++ )
      {
      outImHisto[i].resize( m_HistoNumBinsND, 0 );
      }

    // Create Histogram
    unsigned int totalOut = 0;
    typename ListSampleType::ConstIterator
      outListIt( m_OutList->Begin() );
    typename ListSampleType::ConstIterator
      outListItEnd( m_OutList->End() );
    double binV;
    while( outListIt != outListItEnd )
      {
      for( unsigned int i=0; i<N; i++ )
        {
        binV = outListIt.GetMeasurementVector()[i];
        binV = ( int )( ( binV - m_HistoBinMin[i] )
          * m_HistoBinScale[i] + 0.5 );
        if( binV>m_HistoNumBinsND-1 )
          {
          binV = m_HistoNumBinsND-1;
          }
        else if( binV<0 )
          {
          binV = 0;
          }
        ++( ( outImHisto[i] )[( int )binV] );
        }
      ++totalOut;
      ++outListIt;
      }

    // Clip Histograms
    double outlierRejectPercent = 0.01;
    double totalReject = totalOut * outlierRejectPercent;
    for( unsigned int i=0; i<N; i++ )
      {
      unsigned int count = 0;
      for( unsigned int b=0; b<m_HistoNumBinsND; b++ )
        {
        count += static_cast<unsigned int>( outImHisto[i][b] );
        if( count>=totalReject )
          {
          //outImMin[i] = b/outImScale[i] + outImMin[i];
          outImHisto[i][b] = count-totalReject;
          break;
          }
        else
          {
          outImHisto[i][b] = 0;
          }
        }
      count = 0;
      for( int b=m_HistoNumBinsND-1; b>=0; b-- )
        {
        count += static_cast<unsigned int>( outImHisto[i][b] );
        if( count>=totalReject )
          {
          //outImMax[i] = b/outImScale[i] + outImMin[i];
          outImHisto[i][b] = count-totalReject;
          break;
          }
        else
          {
          outImHisto[i][b] = 0;
          }
        }
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
  size.Fill( m_HistoNumBinsND );

  m_InClassHisto.resize( numClasses );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    m_InClassHisto[c] = HistogramImageType::New();
    m_InClassHisto[c]->SetRegions( size );
    m_InClassHisto[c]->Allocate();
    m_InClassHisto[c]->FillBuffer( 0 );

    unsigned int count = 0;
    typename ListSampleType::ConstIterator
      inClassListIt( m_InClassList[c]->Begin() );
    typename ListSampleType::ConstIterator
      inClassListItEnd( m_InClassList[c]->End() );
    typename HistogramImageType::IndexType indxHisto;
    while( inClassListIt != inClassListItEnd )
      {
      bool valid = true;
      for( unsigned int i=0; i<N; i++ )
        {
        double binV = inClassListIt.GetMeasurementVector()[i];
        binV = ( int )( ( binV - m_HistoBinMin[i] )
          * m_HistoBinScale[i] + 0.5 );
        indxHisto[i] = ( int )binV;
        if( binV<0 || binV>=m_HistoNumBinsND
          || inImHisto[c][i][( int )binV] == 0 )
          {
          valid = false;
          break;
          }
        }
      if( valid )
        {
        ++count;
        ++( m_InClassHisto[c]->GetPixel( indxHisto ) );
        }
      ++inClassListIt;
      }
    std::cout << "m_InClassHisto size = " << count << std::endl;
    }

  m_OutHisto = HistogramImageType::New();
  m_OutHisto->SetRegions( size );
  m_OutHisto->Allocate();
  m_OutHisto->FillBuffer( 0 );
  if( m_OutList->Size() > 0 )
    {
    typename ListSampleType::ConstIterator
      outListIt( m_OutList->Begin() );
    typename ListSampleType::ConstIterator
      outListItEnd( m_OutList->End() );
    typename HistogramImageType::IndexType indxHisto;
    while( outListIt != outListItEnd )
      {
      bool valid = true;
      for( unsigned int i=0; i<N; i++ )
        {
        double binV = outListIt.GetMeasurementVector()[i];
        binV = ( int )( ( binV - m_HistoBinMin[i] )
          * m_HistoBinScale[i] + 0.5 );
        if( binV<0 || binV>=m_HistoNumBinsND
          || outImHisto[i][( int )binV] == 0 )
          {
          valid = false;
          break;
          }
        indxHisto[i] = ( int )binV;
        }
      if( valid )
        {
        ++( m_OutHisto->GetPixel( indxHisto ) );
        }
      ++outListIt;
      }
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
      HistogramImageType > HistoBlurGenType;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      double inPTotal = 0;

      typename HistoBlurGenType::Pointer m_InClassHistoBlurGen =
        HistoBlurGenType::New();
      m_InClassHistoBlurGen->SetInput( m_InClassHisto[c] );
      m_InClassHistoBlurGen->SetVariance( 100 );
      m_InClassHistoBlurGen->Update();
      m_InClassHisto[c] = m_InClassHistoBlurGen->GetOutput();

      itk::ImageRegionIterator<HistogramImageType> m_InClassHistoIt(
        m_InClassHisto[c],
        m_InClassHisto[c]->GetLargestPossibleRegion() );
      m_InClassHistoIt.GoToBegin();
      while( !m_InClassHistoIt.IsAtEnd() )
        {
        double tf = m_InClassHistoIt.Get();
        inPTotal += tf;
        ++m_InClassHistoIt;
        }
      std::cout << "m_InClassHistoTotalP = " << inPTotal << std::endl;

      if( inPTotal > 0 )
        {
        m_InClassHistoIt.GoToBegin();
        while( !m_InClassHistoIt.IsAtEnd() )
          {
          double tf = m_InClassHistoIt.Get();
          m_InClassHistoIt.Set( tf / inPTotal );
          ++m_InClassHistoIt;
          }
        }
      }

    if( m_OutList->Size() > 0 )
      {
      std::cout << "Outside histogram image smoothing..." << std::endl;

      typename HistoBlurGenType::Pointer m_OutHistoBlurGen =
        HistoBlurGenType::New();
      m_OutHistoBlurGen->SetInput( m_OutHisto );
      m_OutHistoBlurGen->SetVariance( 100 );
      m_OutHistoBlurGen->Update();
      m_OutHisto = m_OutHistoBlurGen->GetOutput();

      double outPTotal = 0;
      itk::ImageRegionIterator<HistogramImageType> m_OutHistoIt(
        m_OutHisto, m_OutHisto->GetLargestPossibleRegion() );
      m_OutHistoIt.GoToBegin();
      while( !m_OutHistoIt.IsAtEnd() )
        {
        double tf = m_OutHistoIt.Get();
        outPTotal += tf;
        ++m_OutHistoIt;
        }
      std::cout << "m_OutHistoTotalP = " << outPTotal << std::endl;

      if( outPTotal > 0 )
        {
        m_OutHistoIt.GoToBegin();
        while( !m_OutHistoIt.IsAtEnd() )
          {
          double tf = m_OutHistoIt.Get();
          m_OutHistoIt.Set( tf / outPTotal );
          ++m_OutHistoIt;
          }
        }
      }
    }
  timeCollector.Stop( "HistogramToPDF" );

  timeCollector.Report();
}

template < class ImageT, unsigned int N, class LabelmapT >
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
        binV = ( int )( ( binV - m_HistoBinMin[i] )
          * m_HistoBinScale[i] + 0.5 );
        if( binV<0 || binV>m_HistoNumBinsND-1 )
          {
          valid = false;
          break;
          }
        binIndex[i] = static_cast<long>( binV );
        }
      if( valid )
        {
        double prob = m_InClassHisto[c]->GetPixel( binIndex );
        probIt.Set( prob );
        }
      else
        {
        probIt.Set( 0 );
        }
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

  if( true ) // creating a local context to limit memory footprint
    {
    m_ProbabilityImageVector[numClasses] = ProbabilityImageType::New();
    m_ProbabilityImageVector[numClasses]->SetRegions(
      m_InputVolumeList[0]->GetLargestPossibleRegion() );
    m_ProbabilityImageVector[numClasses]->CopyInformation(
      m_InputVolumeList[0] );
    m_ProbabilityImageVector[numClasses]->Allocate();
    itk::ImageRegionIterator<ProbabilityImageType> probIt(
      m_ProbabilityImageVector[numClasses],
      m_ProbabilityImageVector[numClasses]->GetLargestPossibleRegion() );
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
        binV = ( int )( ( binV - m_HistoBinMin[i] )
          * m_HistoBinScale[i] + 0.5 );
        if( binV<0 || binV>m_HistoNumBinsND-1 )
          {
          valid = false;
          break;
          }
        binIndex[i] = static_cast<long>( binV );
        }
      if( valid )
        {
        double prob = m_OutHisto->GetPixel( binIndex );
        probIt.Set( prob );
        }
      else
        {
        probIt.Set( 0 );
        }
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
      m_ProbabilitySmoothingStandardDeviation );
    probImageFilter->Update();
    m_ProbabilityImageVector[c] = probImageFilter->GetOutput();
    }

  if( true ) // creating a local context to limit memory footprint
    {
    probImageFilter = ProbImageFilterType::New();
    probImageFilter->SetInput( m_ProbabilityImageVector[numClasses] );
    probImageFilter->SetVariance(
      m_ProbabilitySmoothingStandardDeviation );
    probImageFilter->Update();
    m_ProbabilityImageVector[numClasses] = probImageFilter->GetOutput();
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

      itk::ImageRegionIteratorWithIndex<MaskImageType> labelIt(
        tmpLabelImage, tmpLabelImage->GetLargestPossibleRegion() );
      labelIt.GoToBegin();
      while( !labelIt.IsAtEnd() )
        {
        labelImageIndex = labelIt.GetIndex();
        bool maxPC = true;
        double maxP = m_ProbabilityImageVector[c]->GetPixel(
          labelImageIndex );
        for( unsigned int oc=0; oc<numClasses+1; oc++ )
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
      insideConnecter->SetInput( tmpLabelImage );
      insideConnecter->SetLower( 64 );
      insideConnecter->SetUpper( 194 );
      insideConnecter->SetReplaceValue( 255 );
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
        // The pixels with maximum probability for the current
        // class are all set to 128 before the update of the
        // ConnectedThresholdImageFilter, so if the input label
        // map is not to be reclassified, set the pixels belonging
        // to this class to 128 before updating the filter regardless
        // of the probability.  Setting input labels to 255 before the update
        // of the ConnectedThresholdFilter will cause the filter to
        // return only the values at 255 (the input label map).
        if( !m_ReclassifyObjectMask )
          {
          tmpLabelImage->SetPixel( indx, 128 );
          }
        ++inClassListIt;
        }
      insideConnecter->Update();
      tmpLabelImage = insideConnecter->GetOutput();

      timeCollector.Stop( "Connectivity" );

      if( !m_ReclassifyObjectMask )
        {
        // Use inside mask to set seed points.  Also draw inside mask in
        // label image to ensure those points are considered object points.
        // Erase other mask from label image
        for( unsigned int oc=0; oc<numClasses; oc++ )
          {
          if( oc != c )
            {
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
        }

      // Erase outside mask from label image
      if( !m_ReclassifyNotObjectMask )
        {
        typename ListSampleType::ConstIterator
          outListIt( m_OutList->Begin() );
        typename ListSampleType::ConstIterator
          outListItEnd( m_OutList->End() );
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

        typename ConnectedMaskFilterType::Pointer insideConnectedMaskFilter =
          ConnectedMaskFilterType::New();
        insideConnectedMaskFilter->SetInput( tmpLabelImage );
        insideConnectedMaskFilter->SetLower( 194 );
        insideConnectedMaskFilter->SetUpper( 255 );
        insideConnectedMaskFilter->SetReplaceValue( 255 );

        // Use inside mask to set seed points.  Also draw inside mask in
        // label image to ensure those points are considered object points
        inClassListIt = m_InClassList[c]->Begin();
        inClassListItEnd = m_InClassList[c]->End();
        while( inClassListIt != inClassListItEnd )
          {
          for( unsigned int i=0; i<ImageDimension; i++ )
            {
            indx[i] = static_cast<long int>(
              inClassListIt.GetMeasurementVector()[N+i] );
            }

          insideConnectedMaskFilter->AddSeed( indx );
          if( !m_ReclassifyObjectMask )
            {
            tmpLabelImage->SetPixel( indx, 255 );
            }
          ++inClassListIt;
          }

        insideConnectedMaskFilter->Update();
        tmpLabelImage = insideConnectedMaskFilter->GetOutput();

        timeCollector.Stop( "Connectivity2" );
        }

      //
      // Dilate back to original size
      //
      typedef itk::BinaryDilateImageFilter< MaskImageType,
        MaskImageType, StructuringElementType >            DilateFilterType;

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
                if( itInMask.Get() == m_ObjectIdList[ oc ] )
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
          if( itInMask.Get() == m_ObjectIdList[ c ] )
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
            if( itInMask.Get() == m_ObjectIdList[ oc ] )
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

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::Update( void )
{
  this->GenerateSample();
  this->GeneratePDFs();
}

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::ClassifyImages( void )
{
  this->ApplyPDFs();
}

template <class ImageT, unsigned int N, class LabelmapT >
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
  os << indent << "Fpr weight = " << m_FprWeight << std::endl;
  os << indent << "Probability Smoothing Standard Deviation = "
    << m_ProbabilitySmoothingStandardDeviation << std::endl;
  os << indent << "Draft = " << m_Draft << std::endl;
  os << indent << "ReclassifyObjectMask = " << m_ReclassifyObjectMask
    << std::endl;
  os << indent << "ReclassifyNotObjectMask = " << m_ReclassifyNotObjectMask
    << std::endl;
  os << indent << "Number of probability images = "
    << m_ProbabilityImageVector.size() << std::endl;
  os << indent << "InClassList size = "
    << m_InClassList.size() << std::endl;
  if( m_OutList.IsNotNull() )
    {
    os << indent << "OutList = " << m_OutList << std::endl;
    }
  else
    {
    os << indent << "OutList = NULL" << std::endl;
    }
  os << indent << "InClassHisto size = "
    << m_InClassHisto.size() << std::endl;
  if( m_OutHisto.IsNotNull() )
    {
    os << indent << "OutHisto = " << m_OutHisto << std::endl;
    }
  else
    {
    os << indent << "OutHisto = NULL" << std::endl;
    }
  os << indent << "HistoBinMin = " << m_HistoBinMin << std::endl;
  os << indent << "HistoBinMax = " << m_HistoBinMax << std::endl;
  os << indent << "HistoBinScale = " << m_HistoBinScale << std::endl;
  os << indent << "HistoNumBinsND = "
    << m_HistoNumBinsND << std::endl;
  os << indent << "HistoNumBins1D = "
    << m_HistoNumBins1D << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkTubePDFSegmenter_txx)
