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

#ifndef __itktubePDFSegmenterParzen_hxx
#define __itktubePDFSegmenterParzen_hxx

#include "itktubePDFSegmenterParzen.h"
#include "itktubeVectorImageToListGenerator.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itktubeSmoothingRecursiveGaussianImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkHistogram.h>
#include <itkHistogramToProbabilityImageFilter.h>
#include <itkNormalizeToConstantImageFilter.h>
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
PDFSegmenterParzen< TImage, N, TLabelMap >
::PDFSegmenterParzen( void )
{
  m_InClassHistogram.clear();

  m_HistogramBinMin.resize( N, 0 );
  m_HistogramBinSize.resize( N, 0 );
  m_HistogramNumberOfBin.resize( N, 100 );

  m_HistogramSmoothingStandardDeviation = 2;

  m_LabeledFeatureSpace = NULL;
}

template< class TImage, unsigned int N, class TLabelMap >
PDFSegmenterParzen< TImage, N, TLabelMap >
::~PDFSegmenterParzen( void )
{
}

template< class TImage, unsigned int N, class TLabelMap >
typename PDFSegmenterParzen< TImage, N, TLabelMap >::PDFImageType::Pointer
PDFSegmenterParzen< TImage, N, TLabelMap >
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
PDFSegmenterParzen< TImage, N, TLabelMap >
::SetClassPDFImage( unsigned int classNum,
  typename PDFImageType::Pointer classPDF )
{
  if( this->m_ObjectIdList.size() != m_InClassHistogram.size() )
    {
    m_InClassHistogram.resize( this->m_ObjectIdList.size() );
    }
  m_InClassHistogram[classNum] = classPDF;
  this->m_SampleUpToDate = false;
  this->m_PDFsUpToDate = true;
  this->m_ClassProbabilityImagesUpToDate = false;
}

template< class TImage, unsigned int N, class TLabelMap >
const typename PDFSegmenterParzen< TImage, N, TLabelMap >::VectorUIntType &
PDFSegmenterParzen< TImage, N, TLabelMap >
::GetNumberOfBinsPerFeature( void ) const
{
  return m_HistogramNumberOfBin;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::SetNumberOfBinsPerFeature( const VectorUIntType & nBins )
{
  m_HistogramNumberOfBin = nBins;
}

template< class TImage, unsigned int N, class TLabelMap >
const typename PDFSegmenterParzen< TImage, N, TLabelMap >::VectorDoubleType &
PDFSegmenterParzen< TImage, N, TLabelMap >
::GetBinMin( void ) const
{
  return m_HistogramBinMin;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::SetBinMin( const VectorDoubleType & binMin )
{
  m_HistogramBinMin = binMin;
}

template< class TImage, unsigned int N, class TLabelMap >
const typename PDFSegmenterParzen< TImage, N, TLabelMap >::VectorDoubleType &
PDFSegmenterParzen< TImage, N, TLabelMap >
::GetBinSize( void ) const
{
  return m_HistogramBinSize;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::SetBinSize( const VectorDoubleType & scale )
{
  m_HistogramBinSize = scale;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::GeneratePDFs( void )
{
  if( !this->m_SampleUpToDate )
    {
    this->GenerateSample();
    }
  this->m_PDFsUpToDate = true;

  unsigned int numClasses = this->m_ObjectIdList.size();

  //
  // Convert lists to histograms that have the same span ( using range
  //   defined above )
  //

  // Inside

  VectorDoubleType histogramBinMax;
  histogramBinMax.resize( N );
  for( unsigned int i = 0; i < N; i++ )
    {
    m_HistogramBinMin[i] = 99999999999;
    histogramBinMax[i] = -99999999999;
    }

  for( unsigned int c = 0; c < numClasses; c++ )
    {
    typename ListSampleType::ConstIterator
      inClassListIt( this->m_InClassList[c]->Begin() );
    typename ListSampleType::ConstIterator
      inClassListItEnd( this->m_InClassList[c]->End() );
    while( inClassListIt != inClassListItEnd )
      {
      for( unsigned int i = 0; i < N; i++ )
        {
        double binV = inClassListIt.GetMeasurementVector()[i];
        if( binV < m_HistogramBinMin[i] )
          {
          m_HistogramBinMin[i] = binV;
          }
        else if( binV > histogramBinMax[i] )
          {
          histogramBinMax[i] = binV;
          }
        }
      ++inClassListIt;
      }
    }

  for( unsigned int i = 0; i < N; i++ )
    {
    double buffer = 0.025 * ( histogramBinMax[i] - m_HistogramBinMin[i] );
    m_HistogramBinMin[i] -= buffer;
    histogramBinMax[i] += buffer;
    m_HistogramBinSize[i] = ( histogramBinMax[i] - m_HistogramBinMin[i] ) /
      ( double )( m_HistogramNumberOfBin[i] );
    }

  std::vector< VectorDoubleType > clipMin;
  std::vector< VectorDoubleType > clipMax;
  if( true ) // creating a local context to limit memory footprint
    {
    clipMin.resize( numClasses );
    clipMax.resize( numClasses );

    std::vector< vnl_matrix< double > > inImHistogram;
    inImHistogram.resize( numClasses );
    unsigned int maxNB = m_HistogramNumberOfBin[0];
    for( unsigned int i = 1; i < N; ++i )
      {
      if( m_HistogramNumberOfBin[i] > maxNB )
        {
        maxNB = m_HistogramNumberOfBin[i];
        }
      }
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      clipMin[c].resize( N );
      clipMax[c].resize( N );
      inImHistogram[c].set_size( N, maxNB );
      inImHistogram[c].fill( 0 );
      }

    unsigned int totalIn = 0;
    itk::Array<unsigned int> totalInClass( numClasses );
    totalInClass.Fill( 0 );
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      typename ListSampleType::ConstIterator
        inClassListIt( this->m_InClassList[c]->Begin() );
      typename ListSampleType::ConstIterator
        inClassListItEnd( this->m_InClassList[c]->End() );
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i = 0; i < N; i++ )
          {
          double binV = inClassListIt.GetMeasurementVector()[i];
          int binN = static_cast< int >( ( binV - m_HistogramBinMin[i] )
            / m_HistogramBinSize[i] );
          if( binN < 0 )
            {
            binN = 0;
            }
          else if( static_cast< unsigned int >( binN )
            >= m_HistogramNumberOfBin[i] )
            {
            binN = m_HistogramNumberOfBin[i] - 1;
            }
          ++( ( inImHistogram[c][i] )[ binN ] );
          }
        ++totalIn;
        ++totalInClass[c];
        ++inClassListIt;
        }
      }

    for( unsigned int c = 0; c < numClasses; c++ )
      {
      for( unsigned int i = 0; i < N; i++ )
        {
        clipMin[c][i] = m_HistogramBinMin[i];
        clipMax[c][i] = m_HistogramNumberOfBin[i] * m_HistogramBinSize[i]
          + m_HistogramBinMin[i];
        }
      }
    if( m_OutlierRejectPortion > 0 )
      {
      unsigned int count;
      unsigned int prevCount;
      for( unsigned int c = 0; c < numClasses; c++ )
        {
        double tailReject = totalInClass[c] * ( m_OutlierRejectPortion /
          2.0 );
        for( unsigned int i = 0; i < N; i++ )
          {
          count = 0;
          for( unsigned int b = 0; b < m_HistogramNumberOfBin[i]; ++b )
            {
            prevCount = count;
            count += static_cast<unsigned int>( inImHistogram[c][i][b] );
            if( count >= tailReject )
              {
              if( b > 0 )
                {
                clipMin[c][i] = (b-1) * m_HistogramBinSize[i]
                  + m_HistogramBinMin[i];
                }
              if( count - prevCount != 0 )
                {
                clipMin[c][i] += m_HistogramBinSize[i]
                  * ((tailReject-prevCount) / (count-prevCount));
                }
              break;
              }
            }
          count = 0;
          for( int b = ( int )m_HistogramNumberOfBin[i]-1; b >= 0; --b )
            {
            prevCount = count;
            count += static_cast<unsigned int>( inImHistogram[c][i][b] );
            if( count >= tailReject )
              {
              clipMax[c][i] = (b + 1) * m_HistogramBinSize[i]
                + m_HistogramBinMin[i];
              if( count - prevCount != 0 )
                {
                clipMax[c][i] -= m_HistogramBinSize[i]
                  * ((tailReject-prevCount) / (count-prevCount));
                }
              break;
              }
            }
          //std::cout << "Class " << c << " : Feature " << i << " : using "
            //<< clipMin[c][i] << "(" << m_HistogramBinMin[i] << ") - "
            //<< clipMax[c][i] << "(" << histogramBinMax[i] << ")"
            //<< std::endl;
          }
        }
      }

    inImHistogram[0].fill( 0 );
    for( unsigned int i = 0; i < N; i++ )
      {
      m_HistogramBinMin[i] = clipMin[0][i];
      histogramBinMax[i] = clipMax[0][i];
      }
    for( unsigned int c = 1; c < numClasses; c++ )
      {
      inImHistogram[c].fill( 0 );
      for( unsigned int i = 0; i < N; i++ )
        {
        if( clipMin[c][i] < m_HistogramBinMin[i] )
          {
          m_HistogramBinMin[i] = clipMin[c][i];
          }
        else if( clipMax[c][i] > histogramBinMax[i] )
          {
          histogramBinMax[i] = clipMax[c][i];
          }
        }
      }

    for( unsigned int i = 0; i < N; i++ )
      {
      double buffer = 0.025 * ( histogramBinMax[i] -
        m_HistogramBinMin[i] );
      m_HistogramBinMin[i] -= buffer;
      histogramBinMax[i] += buffer;
      m_HistogramBinSize[i] = ( histogramBinMax[i] - m_HistogramBinMin[i] )
        / ( double )( m_HistogramNumberOfBin[i] );
      //std::cout << "Feature " << i << " : buffered : "
        //<< m_HistogramBinMin[i] << " - " << histogramBinMax[i]
        //<< " = " << m_HistogramBinSize[i] << " * "
        //<< m_HistogramNumberOfBin[i] << std::endl;
      }

    for( unsigned int c = 0; c < numClasses; c++ )
      {
      typename ListSampleType::ConstIterator
        inClassListIt( this->m_InClassList[c]->Begin() );
      typename ListSampleType::ConstIterator
        inClassListItEnd( this->m_InClassList[c]->End() );
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i = 0; i < N; i++ )
          {
          double binV = inClassListIt.GetMeasurementVector()[i];
          int binN = static_cast< int >( ( binV - m_HistogramBinMin[i] )
            / m_HistogramBinSize[i] );
          if( binN < 0 )
            {
            binN = 0;
            }
          else if( static_cast< unsigned int >( binN )
            >= m_HistogramNumberOfBin[i] )
            {
            binN = m_HistogramNumberOfBin[i] - 1;
            }
          ++( ( inImHistogram[c][i] )[ binN ] );
          }
        ++totalIn;
        ++totalInClass[c];
        ++inClassListIt;
        }
      }
    }

  //
  //  Create joint histograms
  //
  typename HistogramImageType::SizeType size;
  typename HistogramImageType::SpacingType spacing;
  typename HistogramImageType::PointType origin;
  typename HistogramImageType::RegionType region;
  for( unsigned int i = 0; i < N; i++ )
    {
    spacing[i] = m_HistogramBinSize[i];
    origin[i] = m_HistogramBinMin[i];
    size[i] = m_HistogramNumberOfBin[i];
    }
  region.SetSize( size );

  m_InClassHistogram.resize( numClasses );
  for( unsigned int c = 0; c < numClasses; c++ )
    {
    m_InClassHistogram[c] = HistogramImageType::New();
    m_InClassHistogram[c]->SetRegions( size );
    m_InClassHistogram[c]->SetOrigin( origin );
    m_InClassHistogram[c]->SetSpacing( spacing );
    m_InClassHistogram[c]->Allocate();
    m_InClassHistogram[c]->FillBuffer( 0 );

    typename ListSampleType::ConstIterator
      inClassListIt( this->m_InClassList[c]->Begin() );
    typename ListSampleType::ConstIterator
      inClassListItEnd( this->m_InClassList[c]->End() );
    typename HistogramImageType::IndexType indxHistogram;
    while( inClassListIt != inClassListItEnd )
      {
      for( unsigned int i = 0; i < N; i++ )
        {
        double binV = inClassListIt.GetMeasurementVector()[i];
        int binN = static_cast< int >( ( binV - m_HistogramBinMin[i] )
          / m_HistogramBinSize[i] );
        if( binN < 0 )
          {
          binN = 0;
          }
        else if( static_cast< unsigned int >( binN )
          >= m_HistogramNumberOfBin[i] )
          {
          binN = m_HistogramNumberOfBin[i] - 1;
          }
        indxHistogram[i] = binN;
        }
      ++( m_InClassHistogram[c]->GetPixel( indxHistogram ) );
      ++inClassListIt;
      }
    }

  //
  //  Convert histograms to images so that we can perform blurring to
  //    generate parzen window density estimates
  //
  if( true )
    {
    typedef itk::tube::SmoothingRecursiveGaussianImageFilter<
      HistogramImageType, HistogramImageType > HistogramBlurGenType;
    typename HistogramBlurGenType::SigmaArrayType sigmas;
    for( unsigned int d = 0; d < N; ++d )
      {
      sigmas[d] = spacing[d] * m_HistogramSmoothingStandardDeviation;
      }

    for( unsigned int c = 0; c < numClasses; c++ )
      {
      if( m_HistogramSmoothingStandardDeviation > 0 )
        {
        typename HistogramBlurGenType::Pointer blurFilter =
          HistogramBlurGenType::New();
        blurFilter->SetInput( m_InClassHistogram[c] );
        blurFilter->SetSigmaArray( sigmas );
        blurFilter->Update();
        m_InClassHistogram[c] = blurFilter->GetOutput();

        typedef itk::ThresholdImageFilter< HistogramImageType >
          ThresholdFilterType;
        typename ThresholdFilterType::Pointer thresholdFilter =
          ThresholdFilterType::New();
        thresholdFilter->SetInput( m_InClassHistogram[c] );
        thresholdFilter->SetOutsideValue( 0 );
        thresholdFilter->ThresholdBelow( 0 );
        thresholdFilter->Update();
        m_InClassHistogram[c] = thresholdFilter->GetOutput();
        }

      typedef itk::NormalizeToConstantImageFilter< HistogramImageType,
        HistogramImageType > NormalizeImageFilterType;
      typename NormalizeImageFilterType::Pointer normFilter =
        NormalizeImageFilterType::New();
      normFilter->SetInput( m_InClassHistogram[c] );
      normFilter->Update();
      m_InClassHistogram[c] = normFilter->GetOutput();
      }
    }
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
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
    size[i] = m_HistogramNumberOfBin[i];
    }
  region.SetSize( size );
  m_LabeledFeatureSpace->CopyInformation( m_InClassHistogram[0] );
  m_LabeledFeatureSpace->SetOrigin( origin );
  m_LabeledFeatureSpace->SetRegions( region );
  m_LabeledFeatureSpace->SetSpacing( spacing );
  m_LabeledFeatureSpace->Allocate();

  itk::ImageRegionIterator< LabeledFeatureSpaceType > fsIter(
    m_LabeledFeatureSpace, region );

  unsigned int numClasses = this->m_ObjectIdList.size();

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
    ObjectIdType maxPC = this->m_VoidId;
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      if( pdfIter[c]->Get() > maxP )
        {
        maxP = pdfIter[c]->Get();
        maxPC = this->m_ObjectIdList[ c ];
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
typename PDFSegmenterParzen< TImage, N, TLabelMap >::LabeledFeatureSpaceType::Pointer
PDFSegmenterParzen< TImage, N, TLabelMap >
::GetLabeledFeatureSpace( void ) const
{
  return m_LabeledFeatureSpace;
}

template< class TImage, unsigned int N, class TLabelMap >
typename PDFSegmenterParzen< TImage, N, TLabelMap >::ProbabilityPixelType
PDFSegmenterParzen< TImage, N, TLabelMap >
::GetClassProbability( unsigned int classNum, const FeatureVectorType & fv)
  const
{
  typename HistogramImageType::IndexType binIndex;
  for( unsigned int i = 0; i < N; i++ )
    {
    int binN = static_cast< int >( ( fv[i] - m_HistogramBinMin[i] )
      / m_HistogramBinSize[i] );
    if( binN < 0 )
      {
      binN = 0;
      }
    else if( static_cast< unsigned int >( binN )
      >= m_HistogramNumberOfBin[i] )
      {
      binN = m_HistogramNumberOfBin[i] - 1;
      }
    binIndex[i] = binN;
    }
  return m_InClassHistogram[classNum]->GetPixel( binIndex );
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::SetLabeledFeatureSpace( typename LabeledFeatureSpaceType::Pointer
  labeledFeatureSpace )
{
  m_LabeledFeatureSpace = labeledFeatureSpace;
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::Update( void )
{
  this->GenerateSample();
  this->GeneratePDFs();
  this->GenerateLabeledFeatureSpace();
}

template< class TImage, unsigned int N, class TLabelMap >
void
PDFSegmenterParzen< TImage, N, TLabelMap >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Histogram Smoothing Standard Deviation = "
    << m_HistogramSmoothingStandardDeviation << std::endl;
  os << indent << "InClassHistogram size = "
    << m_InClassHistogram.size() << std::endl;
  os << indent << "HistogramBinMin = " << m_HistogramBinMin[0]
    << std::endl;
  os << indent << "HistogramBinSize = " << m_HistogramBinSize[0]
    << std::endl;
  os << indent << "HistogramNumberOfBin[0] = "
    << m_HistogramNumberOfBin[0] << std::endl;

  os << indent << "Outlier reject portion = "
    << m_OutlierRejectPortion << std::endl;

  if( m_LabeledFeatureSpace.IsNotNull() )
    {
    os << indent << "LabeledFeatureSpace = " << m_LabeledFeatureSpace
      << std::endl;
    }
  else
    {
    os << indent << "LabeledFeatureSpace = NULL" << std::endl;
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubePDFSegmenterParzen_hxx)
