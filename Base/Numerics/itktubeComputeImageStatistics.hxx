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
#ifndef __itktubeComputeImageStatistics_hxx
#define __itktubeComputeImageStatistics_hxx

//ITK includes
#include "itkImageRegionIterator.h"

// TubeTK includes
#include "itktubeComputeImageStatistics.h"

#include <fstream>

namespace itk
{

namespace tube
{
/**
 * Constructor
 */
template< class TPixel, unsigned int VDimension >
ComputeImageStatistics< TPixel, VDimension >
::ComputeImageStatistics()
{
  m_InputMask = NULL;
  m_NumberOfComponents = 0;

}

template< class TPixel, unsigned int VDimension >
void
ComputeImageStatistics< TPixel, VDimension >
::SetQuantiles (const std::vector<float> _arg)
{
  if( this->m_Quantiles != _arg )
    {
    this->m_Quantiles = _arg;
    this->Modified();
    }
}

template< class TPixel, unsigned int VDimension >
void
ComputeImageStatistics< TPixel, VDimension >
::GenerateData( void )
{
  //Sanity checks
  if( !this->GetInput() )
    {
    itkExceptionMacro("Input Image is not set");
    }
  if( !m_InputMask )
    {
    itkExceptionMacro("Input Mask is not set");
    }

  typename ConnCompType::Pointer curConnComp = ConnCompType::New();
  curConnComp->CopyInformation( this->GetInput() );
  curConnComp->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
  curConnComp->Allocate();

  this->GetOutput()->CopyInformation( this->GetInput() );
  this->GetOutput()->SetRegions(
    this->GetInput()->GetLargestPossibleRegion() );
  this->GetOutput()->Allocate();

  typedef std::map< TPixel, unsigned int > MapType;
  MapType maskMap;

  itk::ImageRegionIterator< MaskType > maskIter( m_InputMask,
    m_InputMask->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ConnCompType > connCompIter( curConnComp,
    curConnComp->GetLargestPossibleRegion() );

  // Calculate number of components
  typename MapType::iterator mapIter;
  unsigned int id = 0;
  TPixel lastMaskV = maskIter.Get() + 1;
  TPixel maskV = 0;
  while( !connCompIter.IsAtEnd() )
    {
    maskV = maskIter.Get();

    if( maskV != lastMaskV )
      {
      lastMaskV = maskV;
      mapIter = maskMap.find( maskV );
      if( mapIter != maskMap.end() )
        {
        id = mapIter->second;
        }
      else
        {
        id = m_NumberOfComponents;
        maskMap[ maskV ] = id;
        ++m_NumberOfComponents;
        }
      }
    ++maskIter;
    ++connCompIter;
    }

  //Compute connected components statistics
  m_CompMean.resize( m_NumberOfComponents, 0 );
  m_CompMin.resize( m_NumberOfComponents, 0 );
  m_CompMax.resize( m_NumberOfComponents, 0 );
  m_CompStdDev.resize( m_NumberOfComponents, 0 );
  m_CompCount.resize( m_NumberOfComponents, 0 );
  m_CompValue.resize( m_NumberOfComponents, 0 );

  maskIter.GoToBegin();
  connCompIter.GoToBegin();
  itk::ImageRegionConstIterator< VolumeType > volumeIter( this->GetInput(),
    this->GetInput()->GetLargestPossibleRegion() );
  id = 0;
  lastMaskV = maskIter.Get() + 1;
  maskV = 0;
  while( !connCompIter.IsAtEnd() )
    {
    maskV = maskIter.Get();

    if( maskV != lastMaskV )
      {
      lastMaskV = maskV;
      mapIter = maskMap.find( maskV );
      if( mapIter != maskMap.end() )
        {
        id = mapIter->second;
        }
      else
        {
        id = m_NumberOfComponents;
        maskMap[ maskV ] = id;
        m_CompValue[ id ] = maskV;
        ++m_NumberOfComponents;
        }
      }

    double volumeV = volumeIter.Get();
    m_CompMean[ id ] += volumeV;
    m_CompStdDev[ id ] += volumeV * volumeV;
    if( m_CompCount[ id ] == 0 )
      {
      m_CompMin[ id ] = volumeV;
      m_CompMax[ id ] = volumeV;
      }
    else
      {
      if( volumeV < m_CompMin[ id ] )
        {
        m_CompMin[ id ] = volumeV;
        }
      else
        {
        if( volumeV > m_CompMax[ id ] )
          {
          m_CompMax[ id ] = volumeV;
          }
        }
      }
    ++m_CompCount[ id ];

    connCompIter.Set( id );

    ++maskIter;
    ++volumeIter;
    ++connCompIter;
    }

  for( unsigned int i=0; i<m_NumberOfComponents; ++i )
    {
    if( m_CompCount[ i ] > 0 )
      {
      m_CompStdDev[ i ] = std::sqrt( ( m_CompStdDev[ i ]
        - ( ( m_CompMean[ i ] * m_CompMean[ i ] ) / m_CompCount[ i ] ) )
        / ( m_CompCount[ i ] - 1 ) );
      m_CompMean[ i ] /= m_CompCount[ i ];
      }
    }

  //Compute components histogram
  unsigned int maxNumBins = 2000;
  vnl_matrix< unsigned int > compHisto( m_NumberOfComponents, maxNumBins );
  compHisto.fill( 0 );

  maskIter.GoToBegin();
  volumeIter.GoToBegin();
  while( !volumeIter.IsAtEnd() )
    {
    maskV = maskIter.Get();
    mapIter = maskMap.find( maskV );
    id = mapIter->second;
    int bin = static_cast< int >( ( ( volumeIter.Get() -
      m_CompMin[id] ) /
      ( m_CompMax[id] - m_CompMin[id] ) ) * maxNumBins + 0.5 );
    if( bin < 0 )
      {
      bin = 0;
      }
    else if( bin >= static_cast< int >( maxNumBins ) )
      {
      bin = maxNumBins - 1;
      }
    ++compHisto[id][bin];
    ++maskIter;
    ++volumeIter;
    }

  //Compute quantile value
  unsigned int numberOfQuantiles = m_Quantiles.size();
  m_QuantileValue.set_size( m_NumberOfComponents,
    numberOfQuantiles );
  for( unsigned int quantileNumber=0; quantileNumber<numberOfQuantiles;
    ++quantileNumber )
    {
    double quant = m_Quantiles[ quantileNumber ];
    for( unsigned int comp=0; comp<m_NumberOfComponents; ++comp )
      {
      unsigned int targetCount = static_cast< unsigned int >( quant *
        m_CompCount[ comp ] );
      unsigned int bin = 0;
      unsigned int binCount = 0;
      while( binCount + compHisto[ comp ][ bin ] < targetCount )
        {
        binCount += compHisto[ comp ][ bin ];
        ++bin;
        }
      double binPortion = static_cast< double >( targetCount - binCount )
        / compHisto[ comp ][ bin ];
      double binSize = ( m_CompMax[ comp ] - m_CompMin[ comp ] ) / maxNumBins;
      m_QuantileValue[ comp ][ quantileNumber ] = ( bin + binPortion )
        * binSize + m_CompMin[ comp ];
      }
    }

  //Set output image
  itk::ImageRegionIterator< VolumeType > outputVolumeIter( this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );
  connCompIter.GoToBegin();
  outputVolumeIter.GoToBegin();
  while( !connCompIter.IsAtEnd() )
    {
    id = connCompIter.Get();
    outputVolumeIter.Set( m_CompMean[ id ] );
    ++connCompIter;
    ++outputVolumeIter;
    }
}

template< class TPixel, unsigned int VDimension >
void
ComputeImageStatistics< TPixel, VDimension >
::WriteCSVStatistics( std::string csvStatisticsFile ) const
{
  std::cout << "Number of components = " << m_NumberOfComponents << std::endl;
  std::ofstream writeStream;
  unsigned int numberOfQuantiles = m_Quantiles.size();
  if( ! csvStatisticsFile.empty() )
    {
    writeStream.open( csvStatisticsFile.c_str(),
      std::ios::binary | std::ios::out );
    if( !writeStream.rdbuf()->is_open() )
      {
      std::cerr << "Cannot write to file " << csvStatisticsFile << std::endl;
      return;
      }
    }

  std::cout << "id, Value, Count, Mean, StdDev, Min, Max";
  for( unsigned int q=0; q<numberOfQuantiles; ++q )
    {
    std::cout << ", " << m_Quantiles[ q ];
    }
  std::cout << std::endl;
  if( ! csvStatisticsFile.empty() )
    {
    writeStream << "id, Value, Count, Mean, StdDev, Min, Max";
    for( unsigned int q=0; q<numberOfQuantiles; ++q )
      {
      writeStream << ", " << m_Quantiles[ q ];
      }
    writeStream << std::endl;
    }
  for( unsigned int i=0; i<m_NumberOfComponents; ++i )
    {
    std::cout << i << ", ";
    std::cout << static_cast< double >( m_CompValue[ i ] ) << ", ";
    std::cout << m_CompCount[ i ] << ", ";
    if( ! csvStatisticsFile.empty() )
      {
      writeStream << i << ", ";
      writeStream << static_cast< double >( m_CompValue[ i ] ) << ", ";
      writeStream << m_CompCount[ i ] << ", ";
      }

    std::cout << m_CompMean[ i ] << ", ";
    std::cout << m_CompStdDev[ i ] << ", ";
    std::cout << m_CompMin[ i ] << ", ";
    std::cout << m_CompMax[ i ];
    for( unsigned int q=0; q<numberOfQuantiles; ++q )
      {
      std::cout << ", " << m_QuantileValue[ i ][ q ];
      }
    std::cout << std::endl;
    if( ! csvStatisticsFile.empty() )
      {
      writeStream << m_CompMean[ i ] << ", ";
      writeStream << m_CompStdDev[ i ] << ", ";
      writeStream << m_CompMin[ i ] << ", ";
      writeStream << m_CompMax[ i ];
      for( unsigned int q=0; q<numberOfQuantiles; ++q )
        {
        writeStream << ", " << m_QuantileValue[ i ][ q ];
        }
      writeStream << std::endl;
      }
    }
  if( ! csvStatisticsFile.empty() )
    {
    writeStream.close();
    }
}

template< class TPixel, unsigned int VDimension >
void
ComputeImageStatistics< TPixel, VDimension >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // End namespace tube

} // End namespace itk

#endif
