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

#ifndef __itktubePDFSegmenter_hxx
#define __itktubePDFSegmenter_hxx

#include "itktubePDFSegmenter.h"
#include "itktubeVectorImageToListGenerator.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkHistogram.h>
#include <itkHistogramToProbabilityImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkJoinImageFilter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>

#include <vnl/vnl_matrix.h>

#include <limits>

namespace itk
{

namespace tube
{

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenter< TImage, N, TLabelMap >
::PDFSegmenter( void )
{
  m_SampleUpToDate = false;
  m_PDFsUpToDate = false;
  m_ImagesUpToDate = false;

  m_InClassList.clear();
  m_OutClassList = NULL;
  m_InClassHistogram.clear();

  m_HistogramBinMin.resize( N, 0 );
  m_HistogramBinSize.resize( N, 0 );
  m_HistogramNumberOfBin.resize( N, 100 );

  m_InputImageList.clear();
  m_InputImageList.resize( N, NULL );

  m_LabelMap = NULL;

  m_ObjectIdList.clear();
  m_VoidId = std::numeric_limits< LabelMapPixelType >::max();
  m_PDFWeightList.clear();

  m_ErodeRadius = 1;
  m_HoleFillIterations = 1;
  m_Draft = false;
  m_ProbabilityImageSmoothingStandardDeviation = 1;
  m_HistogramSmoothingStandardDeviation = 1;
  m_OutlierRejectPortion = 0.02;
  m_ReclassifyObjectLabels = false;
  m_ReclassifyNotObjectLabels = false;
  m_ForceClassification = false;

  m_ProbabilityImageVector.resize( 0 );

  m_LabeledFeatureSpace = NULL;

  m_ProgressProcessInfo = NULL;
  m_ProgressFraction = 1.0;
  m_ProgressStart = 0;
}

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenter< TImage, N, TLabelMap >
::~PDFSegmenter( void )
{
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetProgressProcessInformation( void * processInfo, double fraction,
  double start )
{
  m_ProgressProcessInfo = processInfo;
  m_ProgressFraction = fraction;
  m_ProgressStart = start;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::ClearObjectIds( void )
{
  m_ObjectIdList.clear();
  m_PDFWeightList.clear();
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
  m_PDFWeightList.clear();
  m_PDFWeightList.push_back( 1.0 );
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::AddObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.push_back( objectId );
  m_PDFWeightList.push_back( 1.0 );
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetObjectId( const ObjectIdListType objectId )
{
  m_ObjectIdList = objectId;
  if( m_PDFWeightList.size() != m_ObjectIdList.size() )
    {
    m_PDFWeightList.resize( m_ObjectIdList.size(), 1.0 );
    }
}

template< class TImage, unsigned int N, class TLabelMap >
const std::vector< int > &
PDFSegmenter< TImage, N, TLabelMap >
::GetObjectId( void ) const
{
  return m_ObjectIdList;
}

template< class TImage, unsigned int N, class TLabelMap >
unsigned int
PDFSegmenter< TImage, N, TLabelMap >
::GetNumberOfObjectIds( void ) const
{
  return m_ObjectIdList.size();
}

template< class TImage, unsigned int N, class TLabelMap >
unsigned int
PDFSegmenter< TImage, N, TLabelMap >
::GetNumberOfClasses( void ) const
{
  return m_ObjectIdList.size();
}

template< class TImage, unsigned int N, class TLabelMap >
unsigned int
PDFSegmenter< TImage, N, TLabelMap >
::GetObjectNumberFromId( ObjectIdType id ) const
{
  for( unsigned int i = 0; i < m_ObjectIdList.size(); i++ )
    {
    if( m_ObjectIdList[i] == id )
      {
      return i;
      }
    }
  throw;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetObjectPDFWeight( unsigned int num, double weight )
{
  m_PDFWeightList[num] = weight;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetObjectPDFWeight( const VectorDoubleType & weight )
{
  m_PDFWeightList = weight;
}

template< class TImage, unsigned int N, class TLabelMap >
const std::vector< double > &
PDFSegmenter< TImage, N, TLabelMap >
::GetObjectPDFWeight( void ) const
{
  return m_PDFWeightList;
}

template< class TImage, unsigned int N, class TLabelMap >
typename Image< float, TImage::ImageDimension >::Pointer
PDFSegmenter< TImage, N, TLabelMap >
::GetClassProbabilityForInput( unsigned int classNum ) const
{
  if( classNum < m_ProbabilityImageVector.size() )
    {
    return m_ProbabilityImageVector[classNum];
    }
  return NULL;
}

template< class TImage, unsigned int N, class TLabelMap >
typename PDFSegmenter< TImage, N, TLabelMap >::PDFImageType::Pointer
PDFSegmenter< TImage, N, TLabelMap >
::GetClassPDFImage( unsigned int classNum ) const
{
  if( classNum < m_InClassHistogram.size() )
    {
    return m_InClassHistogram[classNum];
    }
  return NULL;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetClassPDFImage( unsigned int classNum,
  typename PDFImageType::Pointer classPDF )
{
  if( m_ObjectIdList.size() != m_InClassHistogram.size() )
    {
    m_InClassHistogram.resize( m_ObjectIdList.size() );
    }
  m_InClassHistogram[classNum] = classPDF;
  m_SampleUpToDate = false;
  m_PDFsUpToDate = true;
  m_ImagesUpToDate = false;
}

template< class TImage, unsigned int N, class TLabelMap >
const std::vector< int > &
PDFSegmenter< TImage, N, TLabelMap >
::GetNumberOfBinsPerFeature( void ) const
{
  return m_HistogramNumberOfBin;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetNumberOfBinsPerFeature( const VectorIntType & nBins )
{
  m_HistogramNumberOfBin = nBins;
}

template< class TImage, unsigned int N, class TLabelMap >
const std::vector< double > &
PDFSegmenter< TImage, N, TLabelMap >
::GetBinMin( void ) const
{
  return m_HistogramBinMin;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetBinMin( const VectorDoubleType & binMin )
{
  m_HistogramBinMin = binMin;
}

template< class TImage, unsigned int N, class TLabelMap >
const std::vector< double > &
PDFSegmenter< TImage, N, TLabelMap >
::GetBinSize( void ) const
{
  return m_HistogramBinSize;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetBinSize( const VectorDoubleType & scale )
{
  m_HistogramBinSize = scale;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetInput( unsigned int featureNum, typename ImageType::Pointer vol )
{
  if( featureNum < N )
    {
    m_InputImageList[featureNum] = vol;
    }
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::GenerateSample( void )
{
  m_SampleUpToDate = true;

  unsigned int numClasses = m_ObjectIdList.size();

  //
  //  Convert in/out images to statistical list using masks
  //
  m_InClassList.resize( numClasses );
  for( unsigned int c = 0; c < numClasses; c++ )
    {
    m_InClassList[c] = ListSampleType::New();
    }
  m_OutClassList  = ListSampleType::New();

  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateSample" );

  typedef itk::ImageRegionConstIteratorWithIndex< LabelMapType >
    ConstLabelMapIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex< ImageType >
    ConstImageIteratorType;

  ConstLabelMapIteratorType itInLabelMap( m_LabelMap,
    m_LabelMap->GetLargestPossibleRegion() );
  itInLabelMap.GoToBegin();

  ConstImageIteratorType * itInIm[N];
  std::vector< double > histogramBinMax;
  histogramBinMax.resize( N, 0 );
  for( unsigned int i = 0; i < N; i++ )
    {
    itInIm[i] = new ConstImageIteratorType( m_InputImageList[i],
      m_InputImageList[i]->GetLargestPossibleRegion() );
    itInIm[i]->GoToBegin();
    m_HistogramBinMin[i] = itInIm[i]->Get();
    histogramBinMax[i] = itInIm[i]->Get();
    }
  ListVectorType v;
  typename LabelMapType::IndexType indx;
  while( !itInLabelMap.IsAtEnd() )
    {
    int val = itInLabelMap.Get();
    indx = itInLabelMap.GetIndex();
    for( unsigned int i = 0; i < N; i++ )
      {
      v[i] = static_cast< PixelType >( itInIm[i]->Get() );
      }
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      v[N+i] = indx[i];
      }
    bool found = false;
    for( unsigned int c = 0; c < numClasses; c++ )
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
    for( unsigned int i = 0; i < N; i++ )
      {
      if( v[i]<m_HistogramBinMin[i] )
        {
        m_HistogramBinMin[i] = v[i];
        }
      else if( v[i] > histogramBinMax[i] )
        {
        histogramBinMax[i] = v[i];
        }
      }
    ++itInLabelMap;
    for( unsigned int i = 0; i < N; i++ )
      {
      ++( *( itInIm[i] ) );
      }
    if( m_Draft )
      {
      for( unsigned int count = 1; count < 4; count++ )
        {
        ++itInLabelMap;
        for( unsigned int i = 0; i < N; i++ )
          {
          ++( *( itInIm[i] ) );
          }
        }
      }
    }

  for( unsigned int i = 0; i < N; i++ )
    {
    m_HistogramBinSize[i] = ( double )m_HistogramNumberOfBin[i] /
      ( histogramBinMax[i] - m_HistogramBinMin[i] );
    }

  for( unsigned int i = 0; i < N; i++ )
    {
    delete itInIm[i];
    }

  timeCollector.Stop( "GenerateSample" );

  timeCollector.Report();
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
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

  std::vector< std::vector< double > > clipMin;
  std::vector< std::vector< double > > clipMax;
  if( true ) // creating a local context to limit memory footprint
    {
    clipMin.resize( numClasses );
    clipMax.resize( numClasses );

    std::vector< vnl_matrix< double > > inImHistogram;
    inImHistogram.resize( numClasses );
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      clipMin[c].resize( N );
      clipMax[c].resize( N );
      unsigned int maxN = m_HistogramNumberOfBin[0];
      for( unsigned int i = 1; i < N; ++i )
        {
        if( m_HistogramNumberOfBin[i] > maxN )
          {
          maxN = m_HistogramNumberOfBin[i];
          }
        }
      inImHistogram[c].set_size( N, maxN );
      inImHistogram[c].fill( 0 );
      }

    unsigned int totalIn = 0;
    itk::Array<unsigned int> totalInClass( numClasses );
    totalInClass.Fill( 0 );
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      typename ListSampleType::ConstIterator
        inClassListIt( m_InClassList[c]->Begin() );
      typename ListSampleType::ConstIterator
        inClassListItEnd( m_InClassList[c]->End() );
      double binV;
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i = 0; i < N; i++ )
          {
          binV = inClassListIt.GetMeasurementVector()[i];
          binV = ( int )( ( binV - m_HistogramBinMin[i] )
            * m_HistogramBinSize[i] + 0.5 );
          if( binV>m_HistogramNumberOfBin[i]-1 )
            {
            binV = m_HistogramNumberOfBin[i]-1;
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
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      double tailReject = totalInClass[c] * ( m_OutlierRejectPortion/2 );
      for( unsigned int i = 0; i < N; i++ )
        {
        count = 0;
        for( unsigned int b = 0; b < m_HistogramNumberOfBin[i]; b++ )
          {
          count += static_cast<unsigned int>( inImHistogram[c][i][b] );
          if( count>=tailReject )
            {
            clipMin[c][i] =  b / m_HistogramBinSize[i]
              + m_HistogramBinMin[i];
            break;
            }
          }
        count = 0;
        for( int b = ( int )m_HistogramNumberOfBin[i]-1; b >= 0; b-- )
          {
          count += static_cast<unsigned int>( inImHistogram[c][i][b] );
          if( count>=tailReject )
            {
            clipMax[c][i] =  b / m_HistogramBinSize[i]
              + m_HistogramBinMin[i];
            break;
            }
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
  for( unsigned int i = 0; i < N; ++i )
    {
    size[i] = m_HistogramNumberOfBin[i];
    }

  m_InClassHistogram.resize( numClasses );
  for( unsigned int c = 0; c < numClasses; c++ )
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
      for( unsigned int i = 0; i < N; i++ )
        {
        double binV = inClassListIt.GetMeasurementVector()[i];
        if( binV >= clipMin[c][i] && binV <= clipMax[c][i] )
          {
          binV = ( int )( ( binV - m_HistogramBinMin[i] )
            * m_HistogramBinSize[i] + 0.5 );
          if( binV<0 || binV>=m_HistogramNumberOfBin[i] )
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
    }
  timeCollector.Stop( "JointClassHistogram" );

  //
  //  Convert histograms to images so that we can perform blurring to
  //    generate parzen window density estimates
  //
  timeCollector.Start( "HistogramToPDF" );
  if( true ) // creating a local context to limit memory footprint
    {
    typedef itk::DiscreteGaussianImageFilter< HistogramImageType,
      HistogramImageType > HistogramBlurGenType;
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      double inPTotal = 0;

      typename HistogramBlurGenType::Pointer inClassHistogramBlurGen =
        HistogramBlurGenType::New();
      inClassHistogramBlurGen->SetInput( m_InClassHistogram[c] );
      inClassHistogramBlurGen->SetVariance(
        m_HistogramSmoothingStandardDeviation *
        m_HistogramSmoothingStandardDeviation );
      inClassHistogramBlurGen->Update();
      m_InClassHistogram[c] = inClassHistogramBlurGen->GetOutput();

      itk::ImageRegionIterator<HistogramImageType> m_InClassHistogramIt(
        m_InClassHistogram[c],
        m_InClassHistogram[c]->GetLargestPossibleRegion() );
      while( !m_InClassHistogramIt.IsAtEnd() )
        {
        double tf = m_InClassHistogramIt.Get();
        inPTotal += tf;
        ++m_InClassHistogramIt;
        }

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

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::GenerateLabeledFeatureSpace( void )
{
  m_LabeledFeatureSpace = LabeledFeatureSpaceType::New();
  typename LabeledFeatureSpaceType::RegionType region;
  typename LabeledFeatureSpaceType::SpacingType spacing;
  typename LabeledFeatureSpaceType::PointType origin;
  typename LabeledFeatureSpaceType::SizeType size;
  for( unsigned int i = 0; i < N; i++ )
    {
    spacing[i] = m_HistogramBinSize[i];
    origin[i] = m_HistogramBinMin[i];
    size[i] = 100;
    }
  region.SetSize( size );
  m_LabeledFeatureSpace->CopyInformation( m_InClassHistogram[0] );
  m_LabeledFeatureSpace->SetOrigin( origin );
  m_LabeledFeatureSpace->SetRegions( region );
  m_LabeledFeatureSpace->SetSpacing( spacing );
  m_LabeledFeatureSpace->Allocate();

  itk::ImageRegionIterator< LabeledFeatureSpaceType > fsIter(
    m_LabeledFeatureSpace, region );

  unsigned int numClasses = m_ObjectIdList.size();

  typedef itk::ImageRegionIterator< HistogramImageType >
    PDFIteratorType;
  std::vector< PDFIteratorType * > pdfIter( numClasses );
  for( unsigned int i = 0; i < numClasses; i++ )
    {
    pdfIter[i] = new PDFIteratorType( m_InClassHistogram[i],
      m_InClassHistogram[i]->GetLargestPossibleRegion() );
    pdfIter[i]->GoToBegin();
    }

  while( !fsIter.IsAtEnd() )
    {
    double maxP = 0;
    ObjectIdType maxPC = m_VoidId;
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      if( pdfIter[c]->Get() > maxP )
        {
        maxP = pdfIter[c]->Get();
        maxPC = m_ObjectIdList[ c ];
        }
      }
    fsIter.Set( maxPC );

    ++fsIter;
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      ++(*(pdfIter[c]));
      }
    }

  for( unsigned int c = 0; c < numClasses; c++ )
    {
    delete pdfIter[c];
    }
}

template< class TImage, unsigned int N, class TLabelMap >
typename PDFSegmenter< TImage, N, TLabelMap >::LabeledFeatureSpaceType::Pointer
PDFSegmenter< TImage, N, TLabelMap >
::GetLabeledFeatureSpace( void ) const
{
  return m_LabeledFeatureSpace;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::SetLabeledFeatureSpace( typename LabeledFeatureSpaceType::Pointer
  labeledFeatureSpace )
{
  m_LabeledFeatureSpace = labeledFeatureSpace;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::ApplyPDFs( void )
{
  if( m_LabelMap.IsNotNull() && !m_SampleUpToDate )
    {
    this->GenerateSample();
    }
  if( m_LabelMap.IsNotNull() && !m_PDFsUpToDate )
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
  m_ProbabilityImageVector.resize( numClasses );

  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "ProbabilityImage" );
  for( unsigned int c = 0; c < numClasses; c++ )
    {
    m_ProbabilityImageVector[c] = ProbabilityImageType::New();
    m_ProbabilityImageVector[c]->SetRegions(
      m_InputImageList[0]->GetLargestPossibleRegion() );
    m_ProbabilityImageVector[c]->CopyInformation( m_InputImageList[0] );
    m_ProbabilityImageVector[c]->Allocate();

    itk::ImageRegionIterator<ProbabilityImageType> probIt(
      m_ProbabilityImageVector[c],
      m_ProbabilityImageVector[c]->GetLargestPossibleRegion() );
    probIt.GoToBegin();

    typedef itk::ImageRegionConstIteratorWithIndex< ImageType >
      ConstImageIteratorType;
    std::vector< ConstImageIteratorType * > itInIm(N);
    for( unsigned int i = 0; i < N; i++ )
      {
      itInIm[i] = new ConstImageIteratorType( m_InputImageList[i],
        m_InputImageList[i]->GetLargestPossibleRegion() );
      itInIm[i]->GoToBegin();
      }

    typename HistogramImageType::IndexType binIndex;
    while( !probIt.IsAtEnd() )
      {
      bool valid = true;
      for( unsigned int i = 0; i < N; i++ )
        {
        double binV = itInIm[i]->Get();
        binV = ( int )( ( binV - m_HistogramBinMin[i] )
          * m_HistogramBinSize[i] + 0.5 );
        if( binV<0 || binV>m_HistogramNumberOfBin[i]-1 )
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
          prob *= m_PDFWeightList[c];
          }
        }
      probIt.Set( prob );
      for( unsigned int i = 0; i < N; i++ )
        {
        ++( *( itInIm[i] ) );
        }
      ++probIt;
      }

    for( unsigned int i = 0; i < N; i++ )
      {
      delete itInIm[i];
      }
    }
  timeCollector.Stop( "ProbabilityImage" );

  typedef itk::DiscreteGaussianImageFilter< ProbabilityImageType,
          ProbabilityImageType> ProbImageFilterType;
  typename ProbImageFilterType::Pointer probImageFilter;

  timeCollector.Start( "ProbabilityImageDiffusion" );
  for( unsigned int c = 0; c < numClasses; c++ )
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

  //
  //  Create label image
  //

  typename LabelMapType::IndexType labelImageIndex;

  typename LabelMapType::Pointer tmpLabelImage = LabelMapType::New();
  tmpLabelImage->SetRegions( m_InputImageList[0]
    ->GetLargestPossibleRegion() );
  tmpLabelImage->CopyInformation( m_InputImageList[0] );
  tmpLabelImage->Allocate();

  if( m_LabelMap.IsNull() )
    {
    m_LabelMap = tmpLabelImage;
    m_LabelMap->FillBuffer( m_VoidId );
    m_ForceClassification = true;
    }

  if( !m_ForceClassification )
    {
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      timeCollector.Start( "Connectivity" );

      // For this class, label all pixels for which it is the most
      // likely class.
      itk::ImageRegionIteratorWithIndex<LabelMapType> labelIt(
        tmpLabelImage, tmpLabelImage->GetLargestPossibleRegion() );
      labelIt.GoToBegin();
      while( !labelIt.IsAtEnd() )
        {
        labelImageIndex = labelIt.GetIndex();
        bool maxPC = true;
        double maxP = m_ProbabilityImageVector[c]->GetPixel(
          labelImageIndex );
        for( unsigned int oc = 0; oc < numClasses; oc++ )
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

      typedef itk::ConnectedThresholdImageFilter<LabelMapType,
        LabelMapType> ConnectedFilterType;
      typename ConnectedFilterType::Pointer insideConnecter =
        ConnectedFilterType::New();

      typename ListSampleType::ConstIterator
        inClassListIt( m_InClassList[c]->Begin() );
      typename ListSampleType::ConstIterator
        inClassListItEnd( m_InClassList[c]->End() );
      typename ImageType::IndexType indx;
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          indx[i] = static_cast<long int>(
            inClassListIt.GetMeasurementVector()[N+i] );
          }
        insideConnecter->AddSeed( indx );
        ++inClassListIt;
        }

      if( !m_ReclassifyObjectLabels )
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
          for( unsigned int i = 0; i < ImageDimension; i++ )
            {
            indx[i] = static_cast<long int>(
              inClassListIt.GetMeasurementVector()[N+i] );
            }
          tmpLabelImage->SetPixel( indx, 128 );
          ++inClassListIt;
          }
        for( unsigned int oc = 0; oc < numClasses; oc++ )
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
              for( unsigned int i = 0; i < ImageDimension; i++ )
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
          for( unsigned int i = 0; i < ImageDimension; i++ )
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
          LabelMapType > HoleFillingFilterType;

        timeCollector.Start( "HoleFiller" );

        typename HoleFillingFilterType::Pointer holeFiller =
          HoleFillingFilterType::New();
        typename LabelMapType::SizeType holeRadius;
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
      typedef BinaryBallStructuringElement< LabelMapPixelType,
        ImageType::ImageDimension >
          StructuringElementType;
      typedef BinaryErodeImageFilter< LabelMapType, LabelMapType,
        StructuringElementType >
          ErodeFilterType;

      StructuringElementType sphereOp;
      if( erodeRadius > 0 )
        {
        timeCollector.Start( "Erode" );

        sphereOp.SetRadius( erodeRadius );
        sphereOp.CreateStructuringElement();

        typename ErodeFilterType::Pointer insideLabelMapErodeFilter =
          ErodeFilterType::New();
        insideLabelMapErodeFilter->SetKernel( sphereOp );
        insideLabelMapErodeFilter->SetErodeValue( 255 );
        insideLabelMapErodeFilter->SetInput( tmpLabelImage );
        insideLabelMapErodeFilter->Update();
        tmpLabelImage = insideLabelMapErodeFilter->GetOutput();

        timeCollector.Stop( "Erode" );
        }

      //
      // Re-do connectivity
      //
      if( true ) // creating a local context to limit memory footprint
        {
        typedef itk::ConnectedThresholdImageFilter<LabelMapType,
          LabelMapType> ConnectedLabelMapFilterType;

        timeCollector.Start( "Connectivity2" );

        typename ConnectedLabelMapFilterType::Pointer
          insideConnectedLabelMapFilter = ConnectedLabelMapFilterType::New();
        insideConnectedLabelMapFilter->SetInput( tmpLabelImage );
        insideConnectedLabelMapFilter->SetLower( 194 );
        insideConnectedLabelMapFilter->SetUpper( 255 );
        insideConnectedLabelMapFilter->SetReplaceValue( 255 );

        // Use inside mask to set seed points.  Also draw inside mask in
        // label image to ensure those points are considered object points
        inClassListIt = m_InClassList[c]->Begin();
        inClassListItEnd = m_InClassList[c]->End();
        if( !m_ReclassifyObjectLabels )
          {
          while( inClassListIt != inClassListItEnd )
            {
            for( unsigned int i = 0; i < ImageDimension; i++ )
              {
              indx[i] = static_cast<long int>(
                inClassListIt.GetMeasurementVector()[N+i] );
              }
            insideConnectedLabelMapFilter->AddSeed( indx );
            tmpLabelImage->SetPixel( indx, 255 );  // Redraw objects
            ++inClassListIt;
            }
          }
        else
          {
          while( inClassListIt != inClassListItEnd )
            {
            for( unsigned int i = 0; i < ImageDimension; i++ )
              {
              indx[i] = static_cast<long int>(
                inClassListIt.GetMeasurementVector()[N+i] );
              }

            insideConnectedLabelMapFilter->AddSeed( indx );
            // Don't redraw objects
            ++inClassListIt;
            }
          }

        insideConnectedLabelMapFilter->Update();
        tmpLabelImage = insideConnectedLabelMapFilter->GetOutput();

        timeCollector.Stop( "Connectivity2" );
        }

      //
      // Dilate back to original size
      //
      typedef itk::BinaryDilateImageFilter< LabelMapType,
        LabelMapType, StructuringElementType >           DilateFilterType;

      if( erodeRadius > 0 )
        {
        timeCollector.Start( "Dilate" );

        typename DilateFilterType::Pointer insideLabelMapDilateFilter =
          DilateFilterType::New();
        insideLabelMapDilateFilter->SetKernel( sphereOp );
        insideLabelMapDilateFilter->SetDilateValue( 255 );
        insideLabelMapDilateFilter->SetInput( tmpLabelImage );
        insideLabelMapDilateFilter->Update();
        tmpLabelImage = insideLabelMapDilateFilter->GetOutput();

        timeCollector.Stop( "Dilate" );
        }

      // Merge with input mask
      typedef itk::ImageRegionIterator< LabelMapType >
        LabelMapIteratorType;
      LabelMapIteratorType itInLabelMap( m_LabelMap,
        m_LabelMap->GetLargestPossibleRegion() );
      itInLabelMap.GoToBegin();

      LabelMapIteratorType itLabel( tmpLabelImage,
        tmpLabelImage->GetLargestPossibleRegion() );
      itLabel.GoToBegin();

      while( !itInLabelMap.IsAtEnd() )
        {
        if( itLabel.Get() == 255 )
          {
          if( itInLabelMap.Get() == m_VoidId
              || ( m_ReclassifyObjectLabels && m_ReclassifyNotObjectLabels ) )
            {
            itInLabelMap.Set( m_ObjectIdList[c] );
            }
          else
            {
            if( m_ReclassifyObjectLabels || m_ReclassifyNotObjectLabels )
              {
              bool isObjectId = false;
              for( unsigned int oc = 0; oc < numClasses; oc++ )
                {
                if( itInLabelMap.Get() == m_ObjectIdList[oc] )
                  {
                  isObjectId = true;
                  break;
                  }
                }
              if( ( isObjectId && m_ReclassifyObjectLabels ) ||
                  ( !isObjectId && m_ReclassifyNotObjectLabels ) )
                {
                itInLabelMap.Set( m_ObjectIdList[c] );
                }
              }
            }
          }
        else
          {
          if( itInLabelMap.Get() == m_ObjectIdList[c] )
            {
            if( m_ReclassifyObjectLabels )
              {
              itInLabelMap.Set( m_VoidId );
              }
            }
          }
        ++itInLabelMap;
        ++itLabel;
        }
      }
    }
  else
    {
    timeCollector.Start( "ForceClassification" );

    // Merge with input mask
    typedef itk::ImageRegionIteratorWithIndex< LabelMapType >
      LabelMapIteratorType;
    LabelMapIteratorType itInLabelMap( m_LabelMap,
      m_LabelMap->GetLargestPossibleRegion() );
    itInLabelMap.GoToBegin();

    while( !itInLabelMap.IsAtEnd() )
      {
      labelImageIndex = itInLabelMap.GetIndex();
      unsigned int maxPC = 0;
      double maxP = m_ProbabilityImageVector[0]->GetPixel(
        labelImageIndex );
      for( unsigned int c = 1; c < numClasses; c++ )
        {
        double p = m_ProbabilityImageVector[c]->GetPixel(
          labelImageIndex );
        if( p > maxP )
          {
          maxP = p;
          maxPC = c;
          }
        }
      if( itInLabelMap.Get() == m_VoidId
        || ( m_ReclassifyObjectLabels && m_ReclassifyNotObjectLabels ) )
        {
        itInLabelMap.Set( m_ObjectIdList[maxPC] );
        }
      else
        {
        if( m_ReclassifyObjectLabels || m_ReclassifyNotObjectLabels )
          {
          bool isObjectId = false;
          for( unsigned int oc = 0; oc < numClasses; oc++ )
            {
            if( itInLabelMap.Get() == m_ObjectIdList[oc] )
              {
              isObjectId = true;
              break;
              }
            }
          if( ( isObjectId && m_ReclassifyObjectLabels ) ||
              ( !isObjectId && m_ReclassifyNotObjectLabels ) )
            {
            itInLabelMap.Set( m_ObjectIdList[maxPC] );
            }
          }
        }
      ++itInLabelMap;
      }

    timeCollector.Stop( "ForceClassification" );
    }

  timeCollector.Report();
  m_ImagesUpToDate = true;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::Update( void )
{
  this->GenerateSample();
  this->GeneratePDFs();
  this->GenerateLabeledFeatureSpace();
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
::ClassifyImages( void )
{
  this->ApplyPDFs();
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenter< TImage, N, TLabelMap >
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

  os << indent << "Input volume list size = " << m_InputImageList.size()
    << std::endl;

  if( m_LabelMap.IsNotNull() )
    {
    os << indent << "LabelMap = " << m_LabelMap << std::endl;
    }
  else
    {
    os << indent << "LabelMap = NULL" << std::endl;
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
  os << indent << "ReclassifyObjectLabels = " << m_ReclassifyObjectLabels
    << std::endl;
  os << indent << "ReclassifyNotObjectLabels = "
    << m_ReclassifyNotObjectLabels << std::endl;
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
  os << indent << "HistogramBinMin = " << m_HistogramBinMin[0]
    << std::endl;
  os << indent << "HistogramBinSize = " << m_HistogramBinSize[0]
    << std::endl;
  os << indent << "HistogramNumberOfBin[0] = "
    << m_HistogramNumberOfBin[0] << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubePDFSegmenter_hxx)
