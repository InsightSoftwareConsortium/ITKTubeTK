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

#ifndef __itktubePDFSegmenterParzen_hxx
#define __itktubePDFSegmenterParzen_hxx

#include "itktubeVectorImageToListGenerator.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
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

template< class TImage, class TLabelMap >
PDFSegmenterParzen< TImage, TLabelMap >
::PDFSegmenterParzen( void )
{
  m_InClassHistogram.clear();

  m_HistogramBinMin.clear();
  m_HistogramBinSize.clear();
  m_HistogramNumberOfBin.clear();

  m_HistogramSmoothingStandardDeviation = 4;

  m_OutlierRejectPortion = 0.001;

  m_LabeledFeatureSpace = NULL;
}

template< class TImage, class TLabelMap >
PDFSegmenterParzen< TImage, TLabelMap >
::~PDFSegmenterParzen( void )
{
}

template< class TImage, class TLabelMap >
typename PDFSegmenterParzen< TImage, TLabelMap >::PDFImageType::Pointer
PDFSegmenterParzen< TImage, TLabelMap >
::GetClassPDFImage( unsigned int classNum ) const
{
  if( classNum < m_InClassHistogram.size() )
    {
    return m_InClassHistogram[classNum];
    }
  return nullptr;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::SetClassPDFImage( unsigned int classNum,
  PDFImageType * classPDF )
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

template< class TImage, class TLabelMap >
const typename PDFSegmenterParzen< TImage, TLabelMap >::VectorUIntType &
PDFSegmenterParzen< TImage, TLabelMap >
::GetNumberOfBinsPerFeature( void ) const
{
  return m_HistogramNumberOfBin;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::SetNumberOfBinsPerFeature( const VectorUIntType & nBins )
{
  m_HistogramNumberOfBin = nBins;
}

template< class TImage, class TLabelMap >
const typename PDFSegmenterParzen< TImage, TLabelMap >::VectorDoubleType &
PDFSegmenterParzen< TImage, TLabelMap >
::GetBinMin( void ) const
{
  return m_HistogramBinMin;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::SetBinMin( const VectorDoubleType & binMin )
{
  m_HistogramBinMin = binMin;
}

template< class TImage, class TLabelMap >
const typename PDFSegmenterParzen< TImage, TLabelMap >::VectorDoubleType &
PDFSegmenterParzen< TImage, TLabelMap >
::GetBinSize( void ) const
{
  return m_HistogramBinSize;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::SetBinSize( const VectorDoubleType & scale )
{
  m_HistogramBinSize = scale;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::GeneratePDFs( void )
{
  //std::cout << "GeneratePDFs" << std::endl;
  if( !this->m_SampleUpToDate )
    {
    this->GenerateSample();
    }
  this->m_PDFsUpToDate = true;

  unsigned int numClasses = this->m_ObjectIdList.size();
  unsigned int numFeatures = this->m_FeatureVectorGenerator->
    GetNumberOfFeatures();

  m_HistogramBinMin.resize( numFeatures, 0 );
  m_HistogramBinSize.resize( numFeatures, 0 );
  m_HistogramNumberOfBin.resize( numFeatures, 50 );

  //
  // Convert lists to histograms that have the same span ( using range
  //   defined above )
  //

  // Inside

  VectorDoubleType histogramBinMax;

  //std::cout << "Init" << std::endl;
  histogramBinMax.resize( numFeatures );
  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    m_HistogramBinMin[i] = 99999999999;
    histogramBinMax[i] = -99999999999;
    }

  for( unsigned int c = 0; c < numClasses; c++ )
    {
    typename ListSampleType::const_iterator
      inClassListIt( this->m_InClassList[c].begin() );
    typename ListSampleType::const_iterator
      inClassListItEnd( this->m_InClassList[c].end() );
    while( inClassListIt != inClassListItEnd )
      {
      for( unsigned int i = 0; i < numFeatures; i++ )
        {
        double binV = ( *inClassListIt )[i];
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

  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    double buffer = 0.1 * ( histogramBinMax[i] - m_HistogramBinMin[i] );
    m_HistogramBinMin[i] -= buffer;
    histogramBinMax[i] += buffer;
    m_HistogramBinSize[i] = ( histogramBinMax[i] - m_HistogramBinMin[i] ) /
      ( double )( m_HistogramNumberOfBin[i] );
    //std::cout << "BinMin[" << i << "] = " << m_HistogramBinMin[i]
      //<< std::endl;
    //std::cout << "  BinMax[" << i << "] = " << histogramBinMax[i]
      //<< std::endl;
    //std::cout << "  BinSize[" << i << "] = " << m_HistogramBinSize[i]
      //<< std::endl;
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
    for( unsigned int i = 1; i < numFeatures; ++i )
      {
      if( m_HistogramNumberOfBin[i] > maxNB )
        {
        maxNB = m_HistogramNumberOfBin[i];
        }
      }
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      clipMin[c].resize( numFeatures );
      clipMax[c].resize( numFeatures );
      inImHistogram[c].set_size( numFeatures, maxNB );
      inImHistogram[c].fill( 0 );
      }

    unsigned int totalIn = 0;
    itk::Array<unsigned int> totalInClass( numClasses );
    totalInClass.Fill( 0 );
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      typename ListSampleType::const_iterator
        inClassListIt( this->m_InClassList[c].begin() );
      typename ListSampleType::const_iterator
        inClassListItEnd( this->m_InClassList[c].end() );
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i = 0; i < numFeatures; i++ )
          {
          double binV = ( *inClassListIt )[i];
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
          ++( inImHistogram[c][i][binN] );
          }
        ++totalIn;
        ++totalInClass[c];
        ++inClassListIt;
        }
      }

    for( unsigned int c = 0; c < numClasses; c++ )
      {
      for( unsigned int i = 0; i < numFeatures; i++ )
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
        //std::cout << "Class tail = " << tailReject << std::endl;
        for( unsigned int i = 0; i < numFeatures; i++ )
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
                clipMin[c][i] = ( b-1 ) * m_HistogramBinSize[i]
                  + m_HistogramBinMin[i];
                }
              if( count - prevCount != 0 )
                {
                //std::cout << "   binMin = " << b << std::endl;
                //std::cout << "   count = " << count << std::endl;
                //std::cout << "   prevcount = " << prevCount << std::endl;
                //std::cout << "   clipMin = " << clipMin[c][i] << std::endl;
                clipMin[c][i] += m_HistogramBinSize[i]
                  * ( ( tailReject-prevCount ) / ( count-prevCount ) );
                //std::cout << "   clipMin = " << clipMin[c][i] << std::endl;
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
              clipMax[c][i] = ( b + 1 ) * m_HistogramBinSize[i]
                + m_HistogramBinMin[i];
              if( count - prevCount != 0 )
                {
                clipMax[c][i] -= m_HistogramBinSize[i]
                  * ( ( tailReject-prevCount ) / ( count-prevCount ) );
                }
              break;
              }
            }
          //std::cout << "Class " << c << " : Feature " << i << " : using "
            //<< clipMin[c][i] << "( " << m_HistogramBinMin[i] << " ) - "
            //<< clipMax[c][i] << "( " << histogramBinMax[i] << " )"
            //<< std::endl;
          }
        }
      }

    inImHistogram[0].fill( 0 );
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      m_HistogramBinMin[i] = clipMin[0][i];
      histogramBinMax[i] = clipMax[0][i];
      }
    for( unsigned int c = 1; c < numClasses; c++ )
      {
      inImHistogram[c].fill( 0 );
      for( unsigned int i = 0; i < numFeatures; i++ )
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

    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      double buffer = 0.1 * ( histogramBinMax[i] - m_HistogramBinMin[i] );
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
      typename ListSampleType::const_iterator
        inClassListIt( this->m_InClassList[c].begin() );
      typename ListSampleType::const_iterator
        inClassListItEnd( this->m_InClassList[c].end() );
      while( inClassListIt != inClassListItEnd )
        {
        for( unsigned int i = 0; i < numFeatures; i++ )
          {
          double binV = ( *inClassListIt )[i];
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
          ++( inImHistogram[c][i][binN] );
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
  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    spacing[i] = m_HistogramBinSize[i];
    if(spacing[i] == 0)
      {
      spacing[i] = 1;
      }
    origin[i] = m_HistogramBinMin[i];
    size[i] = m_HistogramNumberOfBin[i];
    }
  for( unsigned int i = numFeatures; i < PARZEN_MAX_NUMBER_OF_FEATURES;
    ++i )
    {
    spacing[i] = 1;
    origin[i] = 0;
    size[i] = 1;
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

    typename ListSampleType::const_iterator
      inClassListIt( this->m_InClassList[c].begin() );
    typename ListSampleType::const_iterator
      inClassListItEnd( this->m_InClassList[c].end() );
    typename HistogramImageType::IndexType indxHistogram;
    indxHistogram.Fill( 0 );
    while( inClassListIt != inClassListItEnd )
      {
      for( unsigned int i = 0; i < numFeatures; i++ )
        {
        double binV = ( *inClassListIt )[i];
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
      m_InClassHistogram[c]->SetPixel( indxHistogram,
        m_InClassHistogram[c]->GetPixel( indxHistogram ) + 1 );
      ++inClassListIt;
      }
    }

  //
  //  Convert histograms to images so that we can perform blurring to
  //    generate parzen window density estimates
  //
  if( true )
    {
    typedef itk::RecursiveGaussianImageFilter<
      HistogramImageType, HistogramImageType > HistogramBlurGenType;

    for( unsigned int c = 0; c < numClasses; c++ )
      {
      if( m_HistogramSmoothingStandardDeviation > 0 )
        {
        for( unsigned int dir=0; dir<numFeatures; ++dir )
          {
          typename HistogramBlurGenType::Pointer blurFilter =
            HistogramBlurGenType::New();
          blurFilter->SetInput( m_InClassHistogram[c] );
          blurFilter->SetSigma( m_HistogramSmoothingStandardDeviation *
            m_InClassHistogram[c]->GetSpacing()[dir] );
          blurFilter->SetDirection( dir );
          blurFilter->SetZeroOrder();
          blurFilter->SetNormalizeAcrossScale( false );
          blurFilter->Update();
          m_InClassHistogram[c] = blurFilter->GetOutput();
          }

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

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::GenerateLabeledFeatureSpace( void )
{
  unsigned int numFeatures = this->m_FeatureVectorGenerator->
    GetNumberOfFeatures();
  m_LabeledFeatureSpace = LabeledFeatureSpaceType::New();
  typename LabeledFeatureSpaceType::RegionType region;
  typename LabeledFeatureSpaceType::SpacingType spacing;
  typename LabeledFeatureSpaceType::PointType origin;
  typename LabeledFeatureSpaceType::SizeType size;
  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    spacing[i] = m_HistogramBinSize[i];
    if( spacing[i] == 0 )
      {
      spacing[i] = 1;
      }
    origin[i] = m_HistogramBinMin[i];
    size[i] = m_HistogramNumberOfBin[i];
    }
  for( unsigned int i = numFeatures; i < PARZEN_MAX_NUMBER_OF_FEATURES;
    ++i )
    {
    spacing[i] = 1;
    origin[i] = 0;
    size[i] = 1;
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
      ++( *( pdfIter[c] ) );
      }
    }

  for( unsigned int c = 0; c < numClasses; c++ )
    {
    delete pdfIter[c];
    }
}

template< class TImage, class TLabelMap >
typename PDFSegmenterParzen< TImage, TLabelMap >::LabeledFeatureSpaceType::
Pointer
PDFSegmenterParzen< TImage, TLabelMap >
::GetLabeledFeatureSpace( void ) const
{
  return m_LabeledFeatureSpace;
}

template< class TImage, class TLabelMap >
typename PDFSegmenterParzen< TImage, TLabelMap >::ProbabilityVectorType
PDFSegmenterParzen< TImage, TLabelMap >
::GetProbabilityVector( const FeatureVectorType & fv ) const
{
  unsigned int numFeatures = this->m_FeatureVectorGenerator->
    GetNumberOfFeatures();
  typename HistogramImageType::IndexType binIndex;
  binIndex.Fill( 0 );
  for( unsigned int i = 0; i < numFeatures; i++ )
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

  unsigned int numClasses = this->m_ObjectIdList.size();
  ProbabilityVectorType prob( numClasses );
  for( unsigned int c=0; c<numClasses; ++c )
    {
    prob[c] = m_InClassHistogram[c]->GetPixel( binIndex );
    }
  return prob;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::SetLabeledFeatureSpace( LabeledFeatureSpaceType * labeledFeatureSpace )
{
  m_LabeledFeatureSpace = labeledFeatureSpace;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::Update( void )
{
  this->GenerateSample();
  this->GeneratePDFs();
  this->GenerateLabeledFeatureSpace();
}

template< class TImage, class TLabelMap >
void
PDFSegmenterParzen< TImage, TLabelMap >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Histogram Smoothing Standard Deviation = "
    << m_HistogramSmoothingStandardDeviation << std::endl;
  os << indent << "InClassHistogram size = "
    << m_InClassHistogram.size() << std::endl;
  if( m_HistogramBinMin.size() > 0 )
    {
    os << indent << "HistogramBinMin = " << m_HistogramBinMin[0]
      << std::endl;
    os << indent << "HistogramBinSize = " << m_HistogramBinSize[0]
      << std::endl;
    os << indent << "HistogramNumberOfBin[0] = "
      << m_HistogramNumberOfBin[0] << std::endl;
    }
  else
    {
    os << indent << "HistogramBinMin = NULL" << std::endl;
    os << indent << "HistogramBinSize = NULL" << std::endl;
    os << indent << "HistogramNumberOfBin = NULL" << std::endl;
    }

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

#endif // End !defined( __itktubePDFSegmenterParzen_hxx )
