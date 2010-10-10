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
#ifndef __itkPDFSegmenter_txx
#define __itkPDFSegmenter_txx

#include "itkPDFSegmenter.h"

#include "itkVectorImageToListGenerator.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHistogram.h"
#include "itkListSample.h"
#include "itkJoinImageFilter.h"
#include "itkListSampleToHistogramGenerator.h"
#include "itkHistogramToProbabilityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include <limits>

//#include "itkPluginFilterWatcher.h"

#define nBinsND 200
#define nBins1D 200

namespace itk
{

template< class ImageT, unsigned int N, class LabelmapT >
PDFSegmenter< ImageT, N, LabelmapT >
::PDFSegmenter()
{
  m_InputVolume1 = NULL;
  m_InputVolume2 = NULL;
  m_InputVolume3 = NULL;

  m_Labelmap = NULL;

  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( 1 );
  m_VoidId = std::numeric_limits< MaskPixelType >::max();

  m_UseTexture = false;
  m_ErodeRadius = 1;
  m_HoleFillIterations = 1;
  m_FprWeight = 1.5;
  m_Draft = false;
  m_ProbabilitySmoothingStandardDeviation = 9;
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
::~PDFSegmenter()
{
}

template< class ImageT, unsigned int N, class LabelmapT >
typename ImageT::Pointer 
PDFSegmenter< ImageT, N, LabelmapT >
::GenerateTextureImage( const ImageT * im )
{
  typename ImageT::Pointer text = ImageT::New();
  text->SetRegions( im->GetLargestPossibleRegion() );
  text->CopyInformation( im );
  text->Allocate();

  //Neighborhood iterator here
  typedef itk::ConstNeighborhoodIterator< ImageT > NeighIterType;
  typedef itk::ImageRegionIterator< ImageT >       IterType;

  typename NeighIterType::RadiusType radius;
  radius.Fill( 1 );
  NeighIterType inIter( radius, im, im->GetLargestPossibleRegion() );
  unsigned int size = ( int )( inIter.Size() );

  IterType outIter( text, text->GetLargestPossibleRegion() );

  for( inIter.GoToBegin(), outIter.GoToBegin();
      !inIter.IsAtEnd();
      ++inIter, ++outIter )
    {
    PixelType v = inIter.GetPixel( 0 );
    double minV = v;
    double maxV = v;
    double meanV = v;
    for( unsigned int i=1; i<size; i++ )
      {
      v = inIter.GetPixel( i );
      if( v > maxV )
        {
        maxV = v;
        }
      else if( v < minV )
        {
        minV = v;
        }
      meanV += v;
      }
    // robust mean
    meanV -= ( maxV + minV );
    meanV /= ( ( int )size-2 );
    // robust variance
    v = ( PixelType )( maxV - meanV );
    if( meanV - minV > v )
      {
      v = ( PixelType )( meanV - minV );
      }
    outIter.Set( v );
    }

  return text;
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
const typename itk::OrientedImage< float, 
  ::itk::GetImageDimension< ImageT >::ImageDimension >::Pointer *
PDFSegmenter< ImageT, N, LabelmapT >
::GetProbabilityImage( unsigned int classNum )
{
  if( classNum < m_ProbabilityImageVector.size() )
    {
    return &( m_ProbabilityImageVector[classNum] );
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, unsigned int N, class LabelmapT >
void
PDFSegmenter< ImageT, N, LabelmapT >
::Update()
{
  typename ImageType::Pointer     inIm[N];

  inIm[0] = m_InputVolume1;
  if( N > 1 )
    {
    inIm[1] = m_InputVolume2;
    }
  if( N > 2 )
    {
    inIm[2] = m_InputVolume3;
    }

  if( m_UseTexture )
    {
    inIm[N-1] = this->GenerateTextureImage( inIm[0] );
    }

  int erodeRadius = m_ErodeRadius;
  int holeFillIterations = m_HoleFillIterations;
  if( m_Draft )
    {
    erodeRadius /= 2;
    holeFillIterations /= 2;
    }

  unsigned int numClasses = m_ObjectIdList.size();

  //
  //  Convert in/out images to statistical list using masks
  //
  typedef itk::Vector< PixelType, N+ImageDimension >     
    ListVectorType;
  typedef itk::Statistics::ListSample< ListVectorType >
    ListSampleType;
  typedef std::vector< typename ListSampleType::Pointer > 
    ClassListSampleType;
  typedef itk::Vector< double, N > 
    ListDoubleType;

  ClassListSampleType inList;
  inList.resize( numClasses );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    inList[c] = ListSampleType::New();
    }
  ListDoubleType imMin;
  ListDoubleType imMax;

  typename ListSampleType::Pointer outList = ListSampleType::New();

  itk::TimeProbesCollectorBase timeCollector;

  if( true ) // creating a local context to limit memory footprint
    {
    timeCollector.Start( "CreateLists" );
  
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
      itInIm[i] = new ConstImageIteratorType( inIm[i], 
        inIm[i]->GetLargestPossibleRegion() );
      itInIm[i]->GoToBegin();
      imMin[i] = itInIm[i]->Get();
      imMax[i] = itInIm[i]->Get();
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
      for( unsigned int i=0; i<N; i++ )
        {
        if( v[i]<imMin[i] )
          {
          imMin[i] = v[i];
          }
        else if( v[i]>imMax[i] )
          {
          imMax[i] = v[i];
          }
        }
      bool found = false;
      for( unsigned int c=0; c<numClasses; c++ )
        {
        if( val == m_ObjectIdList[c] )
          {
          found = true;
          for( unsigned int i=0; i<ImageDimension; i++ )
            {
            v[N+i] = indx[i];
            }
          inList[c]->PushBack( v );
          }
        }
      if( !found && itInMask.Get() != m_VoidId )
        {
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          v[N+i] = indx[i];
          }
        outList->PushBack( v );
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
  
    timeCollector.Stop( "CreateLists" );
    }

  ListDoubleType imScale;
  for( unsigned int i=0; i<N; i++ )
    {
    imScale[i] = ( double )nBinsND/( imMax[i]-imMin[i] );
    }

  std::cout << "ObjectId[0] = " << m_ObjectIdList[0] << std::endl;
  std::cout << "VoidId = " << m_VoidId << std::endl;
  std::cout << "Inside list ObjectId[0] size = " << inList[0]->Size() 
    << std::endl;
  for( unsigned int i=0; i<N; i++ )
    {
    std::cout << "  Image" << i << ": " << imMin[i] << " - " << imMax[i] 
      << std::endl;
    }
  std::cout << "Outside list size = " << outList->Size() << std::endl;

  //
  // Convert lists to histograms that have the same span ( using range 
  //   defined above )
  //

  timeCollector.Start( "ListsToHistograms" );
  // Inside
  typedef itk::Vector< float, nBinsND > 
    Histo1DType;
  typedef std::vector< Histo1DType > 
    HistoNDType;
  typedef std::vector< HistoNDType > 
    ClassHistoNDType;

  ClassHistoNDType inImHisto;
  if( true ) // creating a local context to limit memory footprint
    {
    HistoNDType hND;
    inImHisto.resize( numClasses, hND );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      Histo1DType h1D;
      h1D.Fill( 0 );
      inImHisto[c].resize( N, h1D );
      for( unsigned int i=0; i<N; i++ )
        {
        inImHisto[c][i].Fill( 0 );
        }
      }
  
    unsigned int totalIn = 0;
    typename ListSampleType::ConstIterator inListIt;
    typename ListSampleType::ConstIterator inListItEnd;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      inListIt = inList[c]->Begin();
      inListItEnd = inList[c]->End();
      double v;
      while( inListIt != inListItEnd )
        {
        for( unsigned int i=0; i<N; i++ )
          {
          v = inListIt.GetMeasurementVector()[i];
          v = ( int )( ( v - imMin[i] )*imScale[i]+0.5 );
          if( v>nBinsND-1 )
            {
            v = nBinsND-1;
            }
          else if( v<0 )
            {
            v = 0;
            }
          ++( ( inImHisto[c][i] )[( int )v] );
          }
        ++totalIn;
        ++inListIt;
        }
      }
    
    double outlierRejectPercent = 0.01;
    double totalReject = totalIn * outlierRejectPercent;
    unsigned int count;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      for( unsigned int i=0; i<N; i++ )
        {
        count = 0;
        for( unsigned int b=0; b<nBinsND; b++ )
          {
          count += static_cast<unsigned int>( inImHisto[c][i][b] );
          if( count>=totalReject )
            {
            //imMin[i] = b/imScale[i] + imMin[i];
            inImHisto[c][i][b] = count - totalReject;
            break;
            }
          else
            {
            inImHisto[c][i][b] = 0;
            }
          }
        count = 0;
        for( int b=( int )nBinsND-1; b>=0; b-- )
          {
          count += static_cast<unsigned int>( inImHisto[c][i][b] );
          if( count>=totalReject )
            {
            //imMax[i] = b/imScale[i] + imMin[i];
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
  itk::Vector<float, nBinsND> outImHisto[N];
  if( true ) // creating a local context to limit memory footprint
    {
    for( unsigned int i=0; i<N; i++ )
      {
      outImHisto[i].Fill( 0 );
      }
  
    // Create Histogram
    unsigned int totalOut = 0;
    typename ListSampleType::ConstIterator outListIt;
    typename ListSampleType::ConstIterator outListItEnd;
    outListIt = outList->Begin();
    outListItEnd = outList->End();
    double v;
    while( outListIt != outListItEnd )
      {
      for( unsigned int i=0; i<N; i++ )
        {
        v = outListIt.GetMeasurementVector()[i];
        v = ( int )( ( v - imMin[i] )*imScale[i]+0.5 );
        if( v>nBinsND-1 )
          {
          v = nBinsND-1;
          }
        else if( v<0 )
          {
          v = 0;
          }
        ++( ( outImHisto[i] )[( int )v] );
        }
      ++totalOut;
      ++outListIt;
      }
  
    // Clip Histograms
    double outlierRejectPercent = 0.01;
    double totalReject = totalOut * outlierRejectPercent;
    unsigned int count;
    for( unsigned int i=0; i<N; i++ )
      {
      count = 0;
      for( unsigned int b=0; b<nBinsND; b++ )
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
      for( int b=nBinsND-1; b>=0; b-- )
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
  timeCollector.Start( "JoinHistogram" );
  typedef itk::Image<float, N> HistogramType;
  typedef std::vector< typename HistogramType::Pointer > 
    ClassHistogramType;

  typedef itk::ImageRegionIteratorWithIndex< HistogramType >  
    HistogramIteratorType;
  typename HistogramType::SizeType size;
  size.Fill( nBinsND );

  ClassHistogramType inHisto;
  inHisto.resize( numClasses );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    inHisto[c] = HistogramType::New();
    inHisto[c]->SetRegions( size );
    inHisto[c]->Allocate();
    inHisto[c]->FillBuffer( 0 );

    unsigned int count = 0;
    typename ListSampleType::ConstIterator inListIt;
    typename ListSampleType::ConstIterator inListItEnd;
    inListIt = inList[c]->Begin();
    inListItEnd = inList[c]->End();
    typename HistogramType::IndexType indxHisto;
    while( inListIt != inListItEnd )
      {
      bool valid = true;
      double v;
      for( unsigned int i=0; i<N; i++ )
        {
        v = inListIt.GetMeasurementVector()[i];
        v = ( int )( ( v-imMin[i] )*imScale[i]+0.5 );
        indxHisto[i] = ( int )v;
        if( v<0 || v>=nBinsND || inImHisto[c][i][( int )v] == 0 )
          {
          valid = false;
          break;
          }
        }
      if( valid )
        {
        ++count;
        ++( inHisto[c]->GetPixel( indxHisto ) );
        }
      ++inListIt;
      }
    std::cout << "inHisto size = " << count << std::endl;
    }

  typename HistogramType::Pointer outHisto = HistogramType::New();
  outHisto->SetRegions( size );
  outHisto->Allocate();
  outHisto->FillBuffer( 0 );
  if( true ) // creating a local context to limit memory footprint
    {
    typename ListSampleType::ConstIterator outListIt;
    typename ListSampleType::ConstIterator outListItEnd;
    outListIt = outList->Begin();
    outListItEnd = outList->End();
    typename HistogramType::IndexType indxHisto;
    while( outListIt != outListItEnd )
      {
      bool valid = true;
      double v;
      for( unsigned int i=0; i<N; i++ )
        {
        v = outListIt.GetMeasurementVector()[i];
        v = ( int )( ( v-imMin[i] )*imScale[i]+0.5 );
        if( v<0 || v>=nBinsND || outImHisto[i][( int )v] == 0 )
          {
          valid = false;
          break;
          }
        indxHisto[i] = ( int )v;
        }
      if( valid )
        {
        ++( outHisto->GetPixel( indxHisto ) );
        }
      ++outListIt;
      }
    }
  timeCollector.Stop( "JoinHistogram" );

  //
  //  Convert histograms to images so that we can perform blurring to
  //    generate parzen window density estimates
  //
  timeCollector.Start( "HistogramToPDF" );
  if( true ) // creating a local context to limit memory footprint
    {
    std::cout << "Inside histogram image smoothing..." << std::endl;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      double inPTotal = 0;
  
      typedef itk::DiscreteGaussianImageFilter< HistogramType, 
        HistogramType > HistoBlurGenType;
      typename HistoBlurGenType::Pointer inHistoBlurGen = 
        HistoBlurGenType::New();
      inHistoBlurGen->SetInput( inHisto[c] );
      inHistoBlurGen->SetVariance( 100 );
      inHistoBlurGen->Update();
      inHisto[c] = inHistoBlurGen->GetOutput();
  
      itk::ImageRegionIterator<HistogramType> inHistoIt( inHisto[c], 
        inHisto[c]->GetLargestPossibleRegion() );
      inHistoIt.GoToBegin();
      while( !inHistoIt.IsAtEnd() )
        {
        double tf = inHistoIt.Get();
        inPTotal += tf;
        ++inHistoIt;
        }
      std::cout << "inHistoTotalP = " << inPTotal << std::endl;
    
      inHistoIt.GoToBegin();
      while( !inHistoIt.IsAtEnd() )
        {
        double tf = inHistoIt.Get();
        inHistoIt.Set( tf / inPTotal );
        ++inHistoIt;
        }
      }
  
    std::cout << "Outside histogram image smoothing..." << std::endl;

    typedef itk::DiscreteGaussianImageFilter< HistogramType,
      HistogramType > HistoBlurGenType;

    typename HistoBlurGenType::Pointer outHistoBlurGen = 
      HistoBlurGenType::New();
    outHistoBlurGen->SetInput( outHisto );
    outHistoBlurGen->SetVariance( 100 );
    outHistoBlurGen->Update();
    outHisto = outHistoBlurGen->GetOutput();
  
    double outPTotal = 0;
    itk::ImageRegionIterator<HistogramType> outHistoIt( outHisto, 
      outHisto->GetLargestPossibleRegion() );
    outHistoIt.GoToBegin();
    while( !outHistoIt.IsAtEnd() )
      {
      double tf = outHistoIt.Get();
      outPTotal += tf;
      ++outHistoIt;
      }
    std::cout << "outHistoTotalP = " << outPTotal << std::endl;

    outHistoIt.GoToBegin();
    while( !outHistoIt.IsAtEnd() )
      {
      double tf = outHistoIt.Get();
      outHistoIt.Set( tf / outPTotal );
      ++outHistoIt;
      }
    }
  timeCollector.Stop( "HistogramToPDF" );
  
  //
  // Save PDFs
  //
  /*
  std::string fileName;
  typedef itk::ImageFileWriter<HistogramType>  HistoImageWriterType;

  std::cout << "Inside PDF saving..." << std::endl;
  typename HistoImageWriterType::Pointer insideHistoImageWriter =
                                               HistoImageWriterType::New();
  fileName = "insideHistoImage.mha";
  insideHistoImageWriter->SetFileName( fileName );
  insideHistoImageWriter->SetInput( inHisto );
  insideHistoImageWriter->Update();

  std::cout << "Outside PDF saving..." << std::endl;
  typename HistoImageWriterType::Pointer outsideHistoImageWriter =
                                               HistoImageWriterType::New();
  fileName = "outsideHistoImage.mha";
  outsideHistoImageWriter->SetFileName( fileName );
  outsideHistoImageWriter->SetInput( outHisto );
  outsideHistoImageWriter->Update();
  */

  // 
  //  Compute the probability at each pixel for input images
  //
  std::cout << "Compute probability image..." << std::endl;

  m_ProbabilityImageVector.resize( numClasses + 1 );

  timeCollector.Start( "ProbabilityImage" );
  for( unsigned int c=0; c<numClasses; c++ )
    {
    m_ProbabilityImageVector[c] = ProbabilityImageType::New();
    m_ProbabilityImageVector[c]->SetRegions( 
      inIm[0]->GetLargestPossibleRegion() );
    m_ProbabilityImageVector[c]->CopyInformation( inIm[0] );
    m_ProbabilityImageVector[c]->Allocate();
   
    itk::ImageRegionIterator<ProbabilityImageType> probIt( 
      m_ProbabilityImageVector[c], 
      m_ProbabilityImageVector[c]->GetLargestPossibleRegion() );
    probIt.GoToBegin();

    typedef itk::ImageRegionConstIteratorWithIndex< ImageType >  
      ConstImageIteratorType;
    ConstImageIteratorType * itInIm[N];
    for( unsigned int i=0; i<N; i++ )
      {
      itInIm[i] = new ConstImageIteratorType( inIm[i], 
        inIm[i]->GetLargestPossibleRegion() );
      itInIm[i]->GoToBegin();
      }
 
    typename HistogramType::IndexType binIndex;
    while( !probIt.IsAtEnd() )
      {
      bool valid = true;
      for( unsigned int i=0; i<N; i++ )
        {
        double v = itInIm[i]->Get();
        v = ( int )( ( v-imMin[i] )*imScale[i]+0.5 );
        if( v<0 || v>nBinsND-1 )
          {
          valid = false;
          break;
          }
        binIndex[i] = static_cast<long>( v );
        }
      if( valid )
        {
        double prob = inHisto[c]->GetPixel( binIndex );
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
      inIm[0]->GetLargestPossibleRegion() );
    m_ProbabilityImageVector[numClasses]->CopyInformation( inIm[0] );
    m_ProbabilityImageVector[numClasses]->Allocate();
    itk::ImageRegionIterator<ProbabilityImageType> probIt( 
      m_ProbabilityImageVector[ numClasses ], 
      m_ProbabilityImageVector[ numClasses ]->GetLargestPossibleRegion() );
    probIt.GoToBegin();
  
    typedef itk::ImageRegionConstIteratorWithIndex< ImageType >  
                                                   ConstImageIteratorType;
    ConstImageIteratorType * itInIm[N];
    for( unsigned int i=0; i<N; i++ )
      {
      itInIm[i] = new ConstImageIteratorType( inIm[i], 
        inIm[i]->GetLargestPossibleRegion() );
      itInIm[i]->GoToBegin();
      }
   
    typename HistogramType::IndexType binIndex;
    while( !probIt.IsAtEnd() )
      {
      bool valid = true;
      for( unsigned int i=0; i<N; i++ )
        {
        double v = itInIm[i]->Get();
        v = ( int )( ( v-imMin[i] )*imScale[i]+0.5 );
        if( v<0 || v>nBinsND-1 )
          {
          valid = false;
          break;
          }
        binIndex[i] = static_cast<long>( v );
        }
      if( valid )
        {
        double prob = outHisto->GetPixel( binIndex );
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
  tmpLabelImage->SetRegions( inIm[0]->GetLargestPossibleRegion() );
  tmpLabelImage->CopyInformation( inIm[0] );
  tmpLabelImage->Allocate();

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
    typename ListSampleType::ConstIterator inListIt;
    typename ListSampleType::ConstIterator inListItEnd;
    inListIt = inList[c]->Begin();
    inListItEnd = inList[c]->End();
    typename ImageType::IndexType indx;
    while( inListIt != inListItEnd )
      {
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        indx[i] = static_cast<long int>( 
          inListIt.GetMeasurementVector()[N+i] );
        }
      insideConnecter->AddSeed( indx );
      if( !m_ReclassifyObjectMask )
        {
        tmpLabelImage->SetPixel( indx, 255 );
        }
      ++inListIt;
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
          inListIt = inList[oc]->Begin();
          inListItEnd = inList[oc]->End();
          while( inListIt != inListItEnd )
            {
            for( unsigned int i=0; i<ImageDimension; i++ )
              {
              indx[i] = static_cast<long int>( 
                inListIt.GetMeasurementVector()[N+i] );
              }
            tmpLabelImage->SetPixel( indx, 0 );
            ++inListIt;
            }
          }
        }
      }

    // Erase outside mask from label image
    if( !m_ReclassifyNotObjectMask )
      {
      typename ListSampleType::ConstIterator outListIt;
      typename ListSampleType::ConstIterator outListItEnd;
      outListIt = outList->Begin();
      outListItEnd = outList->End();
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
    typedef itk::BinaryBallStructuringElement< MaskPixelType, 
      ::itk::GetImageDimension<ImageType>::ImageDimension >
        StructuringElementType;
    typedef itk::BinaryErodeImageFilter< MaskImageType, MaskImageType, 
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
      inListIt = inList[c]->Begin();
      inListItEnd = inList[c]->End();
      while( inListIt != inListItEnd )
        {
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          indx[i] = static_cast<long int>( 
            inListIt.GetMeasurementVector()[N+i] );
          }
  
        insideConnectedMaskFilter->AddSeed( indx );
        if( !m_ReclassifyObjectMask )
          {
          tmpLabelImage->SetPixel( indx, 255 );
          }
        ++inListIt;
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
          itInMask.Set( m_VoidId );
          }
        }
      ++itInMask;
      ++itLabel;
      }
    }

  if( m_ForceClassification )
    {
    // Fix VoidId's here
    }

  timeCollector.Report();
}

template <class ImageT, unsigned int N, class LabelmapT >
void 
PDFSegmenter< ImageT, N, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_InputVolume1.IsNotNull() )
    {
    os << indent << "Input volume 1 = " << m_InputVolume1 << std::endl;
    }
  else
    {
    os << indent << "Input volume 1 = NULL" << std::endl;
    }
  if( m_InputVolume2.IsNotNull() )
    {
    os << indent << "Input volume 2 = " << m_InputVolume2 << std::endl;
    }
  else
    {
    os << indent << "Input volume 2 = NULL" << std::endl;
    }
  if( m_InputVolume3.IsNotNull() )
    {
    os << indent << "Input volume 3 = " << m_InputVolume3 << std::endl;
    }
  else
    {
    os << indent << "Input volume 3 = NULL" << std::endl;
    }
  if( m_Labelmap.IsNotNull() )
    {
    os << indent << "Input volume 3 = " << m_Labelmap << std::endl;
    }
  else
    {
    os << indent << "Input volume 3 = NULL" << std::endl;
    }
  os << indent << "Use texture = " << m_UseTexture << std::endl;
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
}

}

#endif //PDFSegmenter_h
