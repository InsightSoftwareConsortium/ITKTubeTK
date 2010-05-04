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
#ifndef __tubePdPfaScorer_txx
#define __tubePdPfaScorer_txx

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "tubePdPfaScorer.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <iostream>

namespace tube
{

const unsigned char                                  BG = 0;
const unsigned char                                  WIRE = 255;
const unsigned char                                  NOTWIRE = 128;

template< class pixelT, unsigned int dimensionT >
void
PdPfaScorer<pixelT,dimensionT>
::SegmentChanges( typename ImageType::Pointer changeImage,
                  ChangeMapType& changes,
                  LabelMMapType& changesByLabel )
{
  // The mapping of pixels to labels is one-to-one, while, the
  // mapping of labels to pixels is many-to-one
  ChangeMapType changeMap;
  LabelMMapType inverseChangeMap;

  typename SegmentFilter::Pointer segmentFilter = SegmentFilter::New();
  segmentFilter->SetInput( changeImage);
  segmentFilter->Update();

  segmentFilter->SetBackgroundValue( BG );
  typename ImageType::Pointer segmentImage = segmentFilter->GetOutput();

  typename RelabelFilter::Pointer relabelFilter = RelabelFilter::New();

  relabelFilter->SetInput( segmentImage );
  relabelFilter->Update();

  typename ImageType::Pointer relabelImage = relabelFilter->GetOutput();
  ConstIteratorType it( relabelImage,
                        relabelImage->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    typename ImageType::IndexType ind;
    ind = it.GetIndex();
    Pixel2D tmpPixel;
    if( it.Get() > 0 )
      {
      tmpPixel.m_X = ind[0];
      tmpPixel.m_Y = ind[1];
      int label = it.Get();
      changeMap[tmpPixel] = label;
      inverseChangeMap.insert( LabelMMapType::value_type( label, tmpPixel ) );
      }
    }

  changes = changeMap;
  changesByLabel = inverseChangeMap;
}

// Assume that the input images are in unsigned char format
template< class pixelT, unsigned int dimensionT >
void
PdPfaScorer<pixelT,dimensionT>
::ComputeBoundingBoxes( typename ImageType::Pointer im,
                        typename ImageType::Pointer maskImage,
                        BBoxListType& BBList )
{

  int sizex = im->GetLargestPossibleRegion().GetSize()[0];
  int sizey = im->GetLargestPossibleRegion().GetSize()[1];

  ChangeMapType changes;
  LabelMMapType changesByLabel;
  this->SegmentChanges( im , changes, changesByLabel );

  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while( !allSegmentsProcessed )
    {
    std::pair<LabelMMapType::iterator,LabelMMapType::iterator> ret;
    ret = changesByLabel.equal_range( currentSegment );
    if( !( ret.first == ret.second ) )  // Elements with the given label found
      {
      BoundingBox BB;
      BB.m_NumPixels = 0;
      BB.m_StartPixelX = sizex;
      BB.m_EndPixelX = 0;
      BB.m_StartPixelY = sizey;
      BB.m_EndPixelY = 0;

      for( LabelMMapType::iterator it = ret.first; it != ret.second; it++ )
        {
        Pixel2D tmpPixel = (*it).second;

        if( tmpPixel.m_X < BB.m_StartPixelX )
          {
          BB.m_StartPixelX = tmpPixel.m_X;
          }

        if( tmpPixel.m_X > BB.m_EndPixelX )
          {
          BB.m_EndPixelX = tmpPixel.m_X;
          }

        if( tmpPixel.m_Y < BB.m_StartPixelY )
          {
          BB.m_StartPixelY = tmpPixel.m_Y;
          }

        if( tmpPixel.m_Y > BB.m_EndPixelY )
          {
          BB.m_EndPixelY = tmpPixel.m_Y;
          }

        BB.m_NumPixels++;
        }

      BBList.push_back(BB);

      }
    else
      {
      allSegmentsProcessed = true;
      }
    currentSegment++;
    }

}

template< class pixelT, unsigned int dimensionT >
void
PdPfaScorer<pixelT,dimensionT>
::ComputeChangeStatistics( typename ImageType::Pointer trueChangeImage,
                           typename ImageType::Pointer foundChangeImage,
                           int minPixels,
                           int& totalNumChanges,
                           int& totalNumChangesFound,
                           int& totalNumFalsePositives )
{

  ChangeMapType trueChanges;
  LabelMMapType trueChangesByLabel;
  ChangeMapType foundChanges;
  LabelMMapType foundChangesByLabel;

  this->SegmentChanges( trueChangeImage, trueChanges, trueChangesByLabel );
  this->SegmentChanges( foundChangeImage, foundChanges, foundChangesByLabel );

  totalNumChanges = 0;
  totalNumChangesFound = 0;
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while( !allSegmentsProcessed )
    {
    std::pair<LabelMMapType::iterator, LabelMMapType::iterator> ret;
    ret = trueChangesByLabel.equal_range( currentSegment );
    if( !( ret.first == ret.second ) )  // Elements with the given label found
      {

      int trueRegionSize = 0;
      for( LabelMMapType::iterator it = ret.first; it != ret.second; it++ )
        {
        trueRegionSize++;
        }

      if( trueRegionSize >= minPixels )
        {
        totalNumChanges++;

        // Iterate over all of the pixels within this change
        bool changeFound = false;
        for( LabelMMapType::iterator it = ret.first; it != ret.second; it++ )
          {
          if( foundChanges.find( it->second ) != foundChanges.end() )
            {
            changeFound = true;
            break;
            }
          }
        if( changeFound)
          {
          totalNumChangesFound++;
          }
        }
      }
    else
      {
      allSegmentsProcessed = true;
      }
    currentSegment++;
    }

  // Now find false-positives
  totalNumFalsePositives = 0;
  allSegmentsProcessed = false;
  currentSegment = 1;
  while( !allSegmentsProcessed )
    {
    std::pair<LabelMMapType::iterator, LabelMMapType::iterator> ret;
    ret = foundChangesByLabel.equal_range( currentSegment );
    if( !( ret.first == ret.second))  // Elements with the given label found
      {
      // Iterate over all of the pixels within this change
      bool changeFound = false;
      for( LabelMMapType::iterator it = ret.first; it != ret.second; it++ )
        {
        if( trueChanges.find( it->second ) != trueChanges.end() )
          {
          changeFound = true;
          break;
          }
        }
      if( !( changeFound ) )
        {
        totalNumFalsePositives++;
        }
      }
    else
      {
      allSegmentsProcessed = true;
      }
    currentSegment++;
    }
}

}

#endif
