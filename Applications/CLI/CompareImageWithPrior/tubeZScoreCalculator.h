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
#ifndef __tubeZScoreCalculator_h
#define __tubeZScoreCalculator_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <vector>

// progress reporting
#include "tubeCLIProgressReporter.h"

// It is important to use OrientedImages
#include "itkOrientedImage.h"

namespace tube
{

template< class pixelT, unsigned int dimensionT >
class ZScoreCalculator
{
public:

  // input and storage typedefs
  typedef itk::OrientedImage< pixelT, dimensionT >           ImageType;
  typedef pixelT                                             PixelType;
  typedef double                                             PrecisionType;
  typedef itk::Image<PrecisionType, 2>                       HistogramType;
  typedef itk::OrientedImage<bool, dimensionT>               SelectionMaskType;
  typedef std::vector<int>                                   VectorType;
  typedef std::vector<PrecisionType>                         BigVectorType;

  // typedefs for mathematical filters
  typedef itk::DivideByConstantImageFilter< HistogramType, double,
    HistogramType >                                          DividerType;
  typedef itk::MultiplyByConstantImageFilter< HistogramType, double,
    HistogramType >                                          MultiplierType;
  typedef itk::AddImageFilter< HistogramType, HistogramType, HistogramType>
                                                             AdderType;
  typedef itk::SubtractImageFilter< HistogramType, HistogramType,
    HistogramType>                                           SubtracterType;
  typedef itk::SquareImageFilter< HistogramType, 
    HistogramType >                                          SquareType;
  typedef itk::SqrtImageFilter< HistogramType, 
    HistogramType >                                          SqrtType;
  typedef itk::MinimumMaximumImageCalculator<ImageType>      CalculatorType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType>  FullItrType;
  typedef itk::ImageRegionConstIterator<HistogramType>       HistIteratorType;
  typedef itk::ImageRegionConstIterator<SelectionMaskType>   SelectionMaskItrType;

  /// Default Constructor
  ZScoreCalculator();

  /// Default Destructor
  virtual ~ZScoreCalculator();

  /// Update function for doing the processing and producing the output 
  /// volumes.
  void Update(tube::CLIProgressReporter& progressReporter,
              double start,
              double proportion,
              double samples);

  /// Calculate the mean and standard deviation to be used in the scoring
  void CalculateMeanAndStdev(tube::CLIProgressReporter& progressReporter,
                             double start,
                             double proportion,
                             double samples);

  void CalculateRobustMeanAndStdev(tube::CLIProgressReporter& progressReporter,
                                   double start,
                                   double proportion,
                                   double samples,
                                   double percentageToKeep,
                                   double scoreThreshold );

  void SetInputVolume( typename ImageType::Pointer inputVolume );
  void SetInputPrior( typename ImageType::Pointer inputPrior );
  void SetSelectionMask( typename SelectionMaskType::Pointer selectionMask );

  void SetRegion( typename ImageType::RegionType region );

  void SetNumberOfBins( unsigned int numBins );

  void SetInputMin( pixelT min );
  void SetInputMax( pixelT max );
  void SetMaskMin( pixelT min );
  void SetMaskMax( pixelT max );

  void SetROISize( VectorType roiSize );

  void SetMeanHistogram( typename HistogramType::Pointer mean );
  void SetStdevHistogram( typename HistogramType::Pointer stdev );

  typename HistogramType::Pointer GetMeanHistogram();
  typename HistogramType::Pointer GetStdevHistogram();

  typename ImageType::Pointer GetOutputVolume();

protected:
  
  typename ImageType::Pointer            m_InputVolume;
  typename ImageType::Pointer            m_InputPrior;
  typename SelectionMaskType::Pointer    m_SelectionMask;
  typename ImageType::Pointer            m_OutputVolume;
  unsigned int                           m_NumberOfBins;
  typename ImageType::PixelType          m_InputMin;
  typename ImageType::PixelType          m_InputMax;
  typename ImageType::PixelType          m_MaskMin;
  typename ImageType::PixelType          m_MaskMax;
  typename ImageType::RegionType         m_Region;
  VectorType                             m_ROISize;
  typename HistogramType::Pointer        m_Mean;
  typename HistogramType::Pointer        m_Stdev;
  typename HistogramType::Pointer        m_Sum;
  typename HistogramType::Pointer        m_SumSqrs;

  BigVectorType*                         m_ScoreVector;
};

}

#include "tubeZScoreCalculator.txx"

#endif
