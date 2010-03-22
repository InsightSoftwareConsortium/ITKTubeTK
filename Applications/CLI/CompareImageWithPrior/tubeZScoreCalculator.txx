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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "tubeZScoreCalculator.h"

#include "tubeSubImageGenerator.h"
#include "tubeJointHistogramGenerator.h"

#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"

namespace tube
{

template< class pixelT, unsigned int dimensionT>
ZScoreCalculator<pixelT,dimensionT>
::ZScoreCalculator()
  :
  m_InputVolume(NULL),
  m_InputPrior(NULL),
  m_OutputVolume(NULL),
  m_NumberOfBins(100)
{
}

template< class pixelT, unsigned int dimensionT>
ZScoreCalculator<pixelT,dimensionT>
::~ZScoreCalculator()
{
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::Update(tube::CLIProgressReporter& progressReporter,
         double start,
         double proportion,
         double samples)
{
  double progress = start;
  double increment = proportion/samples;

  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> FullItrType;
  typedef itk::ImageRegionConstIterator<HistogramType>      HistIteratorType;
  typedef itk::ImageRegionConstIterator<SelectionMaskType>  SelectionMaskItrType;
  
  FullItrType imageItr( m_InputVolume, m_Region );
  SelectionMaskItrType maskItr( m_SelectionMask, m_Region );

  typename HistogramType::Pointer hist;

  unsigned int dimensions = m_InputVolume->GetImageDimension();

  m_OutputVolume = ImageType::New();
  m_OutputVolume->SetRegions(m_InputVolume->GetLargestPossibleRegion());
  m_OutputVolume->Allocate();
  m_OutputVolume->FillBuffer(0);

  imageItr.GoToBegin();
  maskItr.GoToBegin();
  while( !imageItr.IsAtEnd() )
    {

    typename ImageType::IndexType curIndex = imageItr.GetIndex();
    VectorType roiCenter(dimensions);
    for( unsigned int i = 0; i < dimensionT; ++i )
      {
      roiCenter[i] = curIndex[i];
      }

    if( maskItr.Get() )
      {
      tube::SubImageGenerator<PixelType,dimensionT> subGenerator;
      subGenerator.SetRoiCenter( roiCenter );
      subGenerator.SetRoiSize( m_ROISize );
      subGenerator.SetInputVolume( m_InputVolume );
      subGenerator.SetInputMask( m_InputPrior );
      subGenerator.Update();

      tube::JointHistogramGenerator<PixelType,dimensionT> histGenerator;
      histGenerator.SetInputVolume( subGenerator.GetOutputVolume() );
      histGenerator.SetInputMask( subGenerator.GetOutputMask() );
      histGenerator.SetNumberOfBins( m_NumberOfBins );
      histGenerator.SetInputMin( m_InputMin );
      histGenerator.SetInputMax( m_InputMax );
      histGenerator.SetMaskMin( m_MaskMin );
      histGenerator.SetMaskMax( m_MaskMax );
      histGenerator.Update();
      hist = histGenerator.GetOutputVolume();

      HistIteratorType histItr( hist, hist->GetLargestPossibleRegion() );
      HistIteratorType meanItr( m_Mean, m_Mean->GetLargestPossibleRegion() );
      HistIteratorType stdItr( m_Stdev, m_Stdev->GetLargestPossibleRegion() );
      histItr.GoToBegin();
      meanItr.GoToBegin();
      stdItr.GoToBegin();
      typename HistogramType::PixelType val = 0;
      typename HistogramType::PixelType histCount = 0;
      while( !histItr.IsAtEnd() && !meanItr.IsAtEnd() && !stdItr.IsAtEnd() )
        {
        typename HistogramType::PixelType t = histItr.Get();
        typename HistogramType::PixelType m = meanItr.Get();
        typename HistogramType::PixelType s = stdItr.Get();
        s += 0.0001;
        val += vnl_math_abs(t-m)/s;
        ++histItr;
        ++meanItr;
        ++stdItr;
        ++histCount;
        }
      val /= histCount;
      m_OutputVolume->SetPixel(curIndex,val);
      progress += increment;
      progressReporter.Report( progress );
      }
    ++imageItr;
    ++maskItr;
    }
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetNumberOfBins( unsigned int numBins )
{
  m_NumberOfBins = numBins;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetInputMin( pixelT min )
{
  m_InputMin = min;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetInputMax( pixelT max )
{
  m_InputMax = max;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetMaskMin( pixelT min )
{
  m_MaskMin = min;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetMaskMax( pixelT max )
{
  m_MaskMax = max;
}
 
template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetInputVolume( typename ImageType::Pointer inputVolume )
{
  m_InputVolume = inputVolume;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetInputPrior( typename ImageType::Pointer inputPrior )
{
  m_InputPrior = inputPrior;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetSelectionMask( typename SelectionMaskType::Pointer selectionMask )
{
  m_SelectionMask = selectionMask;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetRegion( typename ImageType::RegionType region )
{
  m_Region = region;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetROISize( VectorType roiSize )
{
  m_ROISize = roiSize;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetMeanHistogram( typename HistogramType::Pointer mean )
{
  m_Mean = mean;
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::SetStdevHistogram( typename HistogramType::Pointer stdev )
{
  m_Stdev = stdev;
}

template< class pixelT, unsigned int dimensionT>
typename ZScoreCalculator<pixelT,dimensionT>::ImageType::Pointer
ZScoreCalculator<pixelT,dimensionT>
::GetOutputVolume()
{
  return m_OutputVolume;
}

} // end namespace tube
