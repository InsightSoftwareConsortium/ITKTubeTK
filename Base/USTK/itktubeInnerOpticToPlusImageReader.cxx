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

#include "itktubeInnerOpticToPlusImageReader.h"
#include <itkMath.h>
#include <itksys/SystemTools.hxx>
#include <itkImageScanlineIterator.h>
#include <itkMetaDataObject.h>


namespace itk
{

namespace tube
{

InnerOpticToPlusImageReader
::InnerOpticToPlusImageReader( void ):
  m_StartIndex( 0 ),
  m_EndIndex( NumericTraits< SizeValueType >::max() ),
  m_IncrementIndex( 1 )
{
  m_SyncRecordManager = new SyncRecordManager;
}


InnerOpticToPlusImageReader
::~InnerOpticToPlusImageReader( void )
{
  delete m_SyncRecordManager;
}


void
InnerOpticToPlusImageReader
::GenerateOutputInformation( void )
{
  OutputImageType * output = this->GetOutput();

  typedef OutputImageType::RegionType RegionType;
  RegionType::SizeType regionSize;
  RegionType::IndexType regionIndex;

  if( m_FileName.length() == 0 )
    {
    itkExceptionMacro( << "Must set the FileName." );
    }

  // SyncRecordManager currently needs the current working directory to be the
  // directory containing the image files.
  const std::string cwdPre = itksys::SystemTools::GetCurrentWorkingDirectory();
  const std::string filenamePath =
    itksys::SystemTools::GetFilenamePath( m_FileName.c_str() );
  // push
  itksys::SystemTools::ChangeDirectory( filenamePath.c_str() );
  const bool loaded = m_SyncRecordManager->load( m_FileName.c_str() );
  if( !loaded )
    {
    itkExceptionMacro( << "Could not load InnerOpticMetadataFile" );
    }

  SyncRecord * syncRecord = m_SyncRecordManager->getNextRecord();
  if( syncRecord == NULL )
    {
    itkExceptionMacro( << "Could not get first Sync Record." );
    }

  SizeValueType xMin = NumericTraits< SizeValueType >::max();
  SizeValueType yMin = NumericTraits< SizeValueType >::max();
  SizeValueType xMax = NumericTraits< SizeValueType >::min();
  SizeValueType yMax = NumericTraits< SizeValueType >::min();
  for( int ii = 0; ii < MAX_US_SCAN_CROP_POLYGON_PTS; ++ii )
    {
    double xReal;
    double yReal;
    const bool gotIt = syncRecord->getScanCropVertex_in_ruf( ii, xReal, yReal );
    if( !gotIt )
      {
      itkExceptionMacro( << "Could not get scan crop vertex" );
      }
    SizeValueType xPixel = Math::RoundHalfIntegerUp< double >( xReal );
    SizeValueType yPixel = Math::RoundHalfIntegerUp< double >( yReal );
    xMin = std::min( xMin, xPixel );
    yMin = std::min( yMin, yPixel );
    xMax = std::max( xMax, xPixel );
    yMax = std::max( yMax, yPixel );
    }

  MetaDataDictionary & metaDataDict = output->GetMetaDataDictionary();
  SizeValueType zCount = 0;
  SizeValueType frameIndex = 0;
  for( SizeValueType ii = 0; ii < m_StartIndex; ++ii )
    {
    syncRecord = m_SyncRecordManager->getNextRecord();
    ++frameIndex;
    }
  while( syncRecord != NULL && frameIndex <= m_EndIndex )
    {
    double transformationMatrix[16];
    syncRecord->getTrackerFromRufMatrix( transformationMatrix );
    std::ostringstream keyPrefix;
    keyPrefix << "Seq_Frame";
    keyPrefix.fill( '0' );
    keyPrefix.width( 4 );
    keyPrefix << zCount;
    std::ostringstream value;
    for( unsigned int ii = 0; ii < 4; ++ii )
      {
      for( unsigned int jj = 0; jj < 4; ++jj )
        {
        // TODO: Use Google DoubleConversion ?
        value << transformationMatrix[jj * 4 + ii];
        if( !( ii == 3 && jj == 3 ) )
          {
          value << ' ';
          }
        }
      }
    EncapsulateMetaData< std::string >( metaDataDict,
      keyPrefix.str() + "_ProbeToTrackerTransform",
      value.str() );
    EncapsulateMetaData< std::string >( metaDataDict,
      keyPrefix.str() + "_ProbeToTrackerTransformStatus",
      "OK" );
    value.str( "" );
    value << syncRecord->getTimestamp();
    EncapsulateMetaData< std::string >( metaDataDict,
                                        keyPrefix.str() + "_Timestamp",
                                        value.str() );

    ++zCount;
    for( SizeValueType ii = 0; ii < m_IncrementIndex; ++ii )
      {
      syncRecord = m_SyncRecordManager->getNextRecord();
      ++frameIndex;
      }
    }

  // pop
  itksys::SystemTools::ChangeDirectory( cwdPre.c_str() );
  regionIndex[0] = xMin;
  regionIndex[1] = yMin;
  regionSize[0] = xMax - xMin;
  regionSize[1] = yMax - yMin;

  regionIndex[2] = 0;
  regionSize[2] = zCount;

  RegionType largestRegion;
  largestRegion.SetIndex( regionIndex );
  largestRegion.SetSize( regionSize );

  output->SetLargestPossibleRegion( largestRegion );
}


void
InnerOpticToPlusImageReader
::GenerateData( void )
{
  const std::string cwdPre = itksys::SystemTools::GetCurrentWorkingDirectory();
  const std::string filenamePath =
    itksys::SystemTools::GetFilenamePath( m_FileName.c_str() );
  // push
  itksys::SystemTools::ChangeDirectory( filenamePath.c_str() );
  m_SyncRecordManager->rewind();
  this->AllocateOutputs();

  OutputImageType * output = this->GetOutput();
  typedef OutputImageType::RegionType RegionType;
  const RegionType region = output->GetBufferedRegion();
  const RegionType::IndexType index = region.GetIndex();
  const RegionType::SizeType size = region.GetSize();

  SyncRecord * syncRecord = m_SyncRecordManager->getNextRecord();
  typedef ImageScanlineIterator< OutputImageType > ImageIteratorType;
  ImageIteratorType outputIt( output, region );
  outputIt.GoToBegin();
  const SizeValueType pixelBytes = 3 * sizeof( PixelComponentType );
  const SizeValueType rufXWidth = 960; // currently hard coded
  const SizeValueType rufXWidthBytes = rufXWidth * pixelBytes;
  //const SizeValueType rufYHeight = 768;
  const SizeValueType xPrePaddingBytes = pixelBytes * index[0];
  const SizeValueType xPostPaddingBytes =
    pixelBytes * ( rufXWidth - index[0] - size[0] );
  SizeValueType frameIndex = 0;
  for( SizeValueType ii = 0; ii < m_StartIndex; ++ii )
    {
    syncRecord = m_SyncRecordManager->getNextRecord();
    ++frameIndex;
    }
  while( syncRecord != NULL && !outputIt.IsAtEnd() )
    {
    const unsigned char * rgbRUFPixelsIt = syncRecord->loadRawRgbPixels();
    rgbRUFPixelsIt += rufXWidthBytes * index[1];
    for( SizeValueType yCount = 0; yCount < size[1]; ++yCount )
      {
      rgbRUFPixelsIt += xPrePaddingBytes;
      for( SizeValueType xCount = 0; xCount < size[0]; ++xCount )
        {
        outputIt.Set( rgbRUFPixelsIt );
        ++outputIt;
        rgbRUFPixelsIt += pixelBytes;
        }
      rgbRUFPixelsIt += xPostPaddingBytes;
      outputIt.NextLine();
      }
    // TODO: SyncRecord avoid allocation/deallocation?
    syncRecord->unloadRawRgbPixels();

    for( SizeValueType ii = 0; ii < m_IncrementIndex; ++ii )
      {
      syncRecord = m_SyncRecordManager->getNextRecord();
      ++frameIndex;
      }
    }

  // pop
  itksys::SystemTools::ChangeDirectory( cwdPre.c_str() );
}


void
InnerOpticToPlusImageReader
::EnlargeOutputRequestedRegion( DataObject * output )
{
  typedef OutputImageType::RegionType RegionType;
  OutputImageType * outputImage =
    itkDynamicCastInDebugMode< OutputImageType * >( output );

  const RegionType & largestRegion = outputImage->GetLargestPossibleRegion();
  const RegionType::IndexType & largestIndex = largestRegion.GetIndex();
  const RegionType::SizeType & largestSize = largestRegion.GetSize();

  RegionType requestedRegion = outputImage->GetRequestedRegion();
  RegionType::IndexType requestedIndex = requestedRegion.GetIndex();
  RegionType::SizeType requestedSize = requestedRegion.GetSize();
  for( unsigned int ii = 0; ii < 2; ++ii )
    {
    requestedIndex[ii] = largestIndex[ii];
    requestedSize[ii] = largestSize[ii];
    }
  requestedRegion.SetIndex( requestedIndex );
  requestedRegion.SetSize( requestedSize );
  outputImage->SetRequestedRegion( requestedRegion );
}

} // End namespace tubetk

} // End namespace itk
