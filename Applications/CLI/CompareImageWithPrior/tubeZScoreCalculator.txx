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

#include <algorithm>

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
  m_ScoreVector = NULL;
}

template< class pixelT, unsigned int dimensionT>
ZScoreCalculator<pixelT,dimensionT>
::~ZScoreCalculator()
{
  if( m_ScoreVector )
    {
    delete m_ScoreVector;
    }
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

  FullItrType imageItr( m_InputVolume, m_Region );
  SelectionMaskItrType maskItr( m_SelectionMask, m_Region );

  typename HistogramType::Pointer hist;

  unsigned int dimensions = m_InputVolume->GetImageDimension();

  m_OutputVolume = ImageType::New();
  m_OutputVolume->CopyInformation( m_InputVolume );
  m_OutputVolume->Allocate();
  m_OutputVolume->FillBuffer(0);

  imageItr.GoToBegin();
  maskItr.GoToBegin();
  if( m_ScoreVector )
    {
    delete m_ScoreVector;
    m_ScoreVector = new BigVectorType(samples);
    }
  else
    {
    m_ScoreVector = new BigVectorType(samples);
    }
  double count = 0;
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
        s += 0.001;
        val += vnl_math_abs((t-m)/s);
        ++histItr;
        ++meanItr;
        ++stdItr;
        ++histCount;
        }
      val /= histCount;
      (*m_ScoreVector)[count] = val;
      m_OutputVolume->SetPixel(curIndex,val);
      progress += increment;
      count++;
      progressReporter.Report( progress );
      }
    ++imageItr;
    ++maskItr;
    }
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::CalculateRobustMeanAndStdev(tube::CLIProgressReporter& progressReporter,
                              double start,
                              double proportion,
                              double samples,
                              double percentageToKeep=100,
                              double scoreThreshold=0)
{
  PrecisionType newSampleNumber = samples;
  PrecisionType upperZScore = 0;
  if( percentageToKeep > 0 && percentageToKeep < 100 )
    {
    std::sort( m_ScoreVector->begin(), m_ScoreVector->end() );
    double fraction = (percentageToKeep/100) * samples;
    upperZScore = (*m_ScoreVector)[static_cast<unsigned int>(fraction)];
    }
  upperZScore = scoreThreshold;

  double progress = start;
  double increment = proportion/samples;

  FullItrType imageItr( m_InputVolume, m_Region );
  SelectionMaskItrType maskItr( m_SelectionMask, m_Region );

  typename HistogramType::Pointer hist;

  imageItr.GoToBegin();
  maskItr.GoToBegin();
  if( m_ScoreVector )
    {
    delete m_ScoreVector;
    m_ScoreVector = new BigVectorType(samples);
    }
  double count = 0;
  while( !imageItr.IsAtEnd() )
    {

    typename ImageType::IndexType curIndex = imageItr.GetIndex();

    VectorType roiCenter(dimensionT);
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
      if( val > upperZScore )
        {
        typename SubtracterType::Pointer subtracter = SubtracterType::New();
        typename SquareType::Pointer square = SquareType::New();
        square->SetInput(hist);
        square->Update();
        subtracter->SetInput1( m_SumSqrs );
        subtracter->SetInput2( square->GetOutput() );
        subtracter->Update();
        m_SumSqrs = subtracter->GetOutput();
        subtracter = SubtracterType::New();
        subtracter->SetInput1( m_Sum );
        subtracter->SetInput2( hist );
        subtracter->Update();
        m_Sum = subtracter->GetOutput();
        newSampleNumber--;
        }
      progress += increment;
      count++;
      progressReporter.Report( progress );
      }
    ++imageItr;
    ++maskItr;
    }

  // Calculate the mean
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput( m_Sum );
  divider->SetConstant( newSampleNumber );
  divider->Update();
  m_Mean = divider->GetOutput();
  
  // Calculate the standard deviation
  typename SubtracterType::Pointer subtracter = SubtracterType::New();
  typename SquareType::Pointer square = SquareType::New();
  typename MultiplierType::Pointer multiplier = MultiplierType::New();
  typename SqrtType::Pointer sqrt = SqrtType::New();
  divider = DividerType::New();
  typename HistogramType::Pointer meanSquaredDivided;
  typename HistogramType::Pointer sumSquaresDivided;
  typename HistogramType::PixelType meanCo = 
    newSampleNumber/(newSampleNumber-1);
  square->SetInput(m_Mean);
  square->Update();
  multiplier->SetInput(square->GetOutput());
  multiplier->SetConstant( meanCo );
  multiplier->Update();
  meanSquaredDivided = multiplier->GetOutput();
  divider->SetInput(m_SumSqrs);
  divider->SetConstant(newSampleNumber-1);
  divider->Update();
  sumSquaresDivided = divider->GetOutput();
  subtracter->SetInput1(sumSquaresDivided);
  subtracter->SetInput2(meanSquaredDivided);
  subtracter->Update();
  sqrt->SetInput(subtracter->GetOutput());
  sqrt->Update();
  m_Stdev = sqrt->GetOutput();  
}

template< class pixelT, unsigned int dimensionT>
void
ZScoreCalculator<pixelT,dimensionT>
::CalculateMeanAndStdev(tube::CLIProgressReporter& progressReporter,
                        double start,
                        double proportion,
                        double samples)
{
  double progress = start;
  double increment = proportion/samples;

  typename HistogramType::Pointer hist;
  typename HistogramType::Pointer meanHist;
  typename HistogramType::Pointer stdevHist;
  typename HistogramType::Pointer sumHist;
  typename HistogramType::Pointer sumSqrHist;

  FullItrType imageItr( m_InputVolume, m_Region );
  SelectionMaskItrType maskItr( m_SelectionMask, m_Region );

  imageItr.GoToBegin();
  maskItr.GoToBegin();
  bool firstPass = true;
  while( !imageItr.IsAtEnd() )
    {
    if( maskItr.Get() )
      {
      typename ImageType::IndexType curIndex = imageItr.GetIndex();
      std::vector<int> roiCenter = std::vector<int>(dimensionT);
      for( unsigned int i = 0; i < dimensionT; ++i )
        {
        roiCenter[i] = curIndex[i];
        }

      tube::SubImageGenerator<PixelType,dimensionT> subGenerator;
      subGenerator.SetRoiCenter(roiCenter);
      subGenerator.SetRoiSize(m_ROISize);
      subGenerator.SetInputVolume(m_InputVolume);
      subGenerator.SetInputMask(m_InputPrior);
      subGenerator.Update();

      tube::JointHistogramGenerator<PixelType,dimensionT> histGenerator;
      histGenerator.SetInputVolume(subGenerator.GetOutputVolume());
      histGenerator.SetInputMask(subGenerator.GetOutputMask());
      histGenerator.SetNumberOfBins(m_NumberOfBins);
      histGenerator.SetInputMin(m_InputMin);
      histGenerator.SetInputMax(m_InputMax);
      histGenerator.SetMaskMin(m_MaskMin);
      histGenerator.SetMaskMax(m_MaskMax);
      histGenerator.Update();
      hist = histGenerator.GetOutputVolume();
      
      if( firstPass )
        {
        typename HistogramType::RegionType histRegion =
          hist->GetLargestPossibleRegion();
        sumHist = HistogramType::New();
        sumHist->SetRegions(histRegion);
        sumHist->Allocate();
        sumHist->FillBuffer(0);
        sumSqrHist = HistogramType::New();
        sumSqrHist->SetRegions(histRegion);
        sumSqrHist->Allocate();
        sumSqrHist->FillBuffer(0);
        firstPass = false;
        }
      
      typename AdderType::Pointer adder = AdderType::New();
      typename SquareType::Pointer square = SquareType::New();
      
      // Calculate Running Sum
      adder->SetInput1(sumHist);
      adder->SetInput2(hist);
      adder->Update();
      sumHist = adder->GetOutput();
      
      // Calculate Running Sum of Squares
      adder = AdderType::New();
      square->SetInput(hist);
      square->Update();
      adder->SetInput1(sumSqrHist);
      adder->SetInput2(square->GetOutput());
      adder->Update();
      sumSqrHist = adder->GetOutput();
      
      progress += increment;
      progressReporter.Report( progress );
      }
    
    ++imageItr;
    ++maskItr;
    }
  
  // Calculate the mean
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput( sumHist );
  divider->SetConstant( samples );
  divider->Update();
  meanHist = divider->GetOutput();
  
  // Calculate the standard deviation
  typename SubtracterType::Pointer subtracter = SubtracterType::New();
  typename SquareType::Pointer square = SquareType::New();
  typename MultiplierType::Pointer multiplier = MultiplierType::New();
  typename SqrtType::Pointer sqrt = SqrtType::New();
  divider = DividerType::New();
  typename HistogramType::Pointer meanSquaredDivided;
  typename HistogramType::Pointer sumSquaresDivided;
  typename HistogramType::PixelType meanCo = samples/(samples-1);
  square->SetInput(meanHist);
  square->Update();
  multiplier->SetInput(square->GetOutput());
  multiplier->SetConstant( meanCo );
  multiplier->Update();
  meanSquaredDivided = multiplier->GetOutput();
  divider->SetInput(sumSqrHist);
  divider->SetConstant(samples-1);
  divider->Update();
  sumSquaresDivided = divider->GetOutput();
  subtracter->SetInput1(sumSquaresDivided);
  subtracter->SetInput2(meanSquaredDivided);
  subtracter->Update();
  sqrt->SetInput(subtracter->GetOutput());
  sqrt->Update();
  stdevHist = sqrt->GetOutput();
  
  m_Sum = sumHist;
  m_SumSqrs = sumSqrHist;

  m_Mean = meanHist;
  m_Stdev = stdevHist;
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
ZScoreCalculator<pixelT,dimensionT>::HistogramType::Pointer 
ZScoreCalculator<pixelT,dimensionT>
::GetMeanHistogram()
{
  return m_Mean;
}

template< class pixelT, unsigned int dimensionT>
ZScoreCalculator<pixelT,dimensionT>::HistogramType::Pointer 
ZScoreCalculator<pixelT,dimensionT>
::GetStdevHistogram()
{
  return m_Stdev;
}

template< class pixelT, unsigned int dimensionT>
typename ZScoreCalculator<pixelT,dimensionT>::ImageType::Pointer
ZScoreCalculator<pixelT,dimensionT>
::GetOutputVolume()
{
  return m_OutputVolume;
}

} // end namespace tube
