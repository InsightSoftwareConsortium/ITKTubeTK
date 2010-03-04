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

#include "tubeJointHistogramGenerator.h"

#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"

namespace tube
{

template< class pixelT, unsigned int dimensionT>
JointHistogramGenerator<pixelT,dimensionT>
::JointHistogramGenerator()
  :
  m_InputVolume(NULL),
  m_InputMask(NULL),
  m_OutputVolume(NULL),
  m_NumberOfBins(100)
{
}

template< class pixelT, unsigned int dimensionT>
JointHistogramGenerator<pixelT,dimensionT>
::~JointHistogramGenerator()
{
}

template< class pixelT, unsigned int dimensionT>
void
JointHistogramGenerator<pixelT,dimensionT>
::Update()
{
  typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
  typedef itk::Image<unsigned int, 2>                   JointHistogramType;

  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetInput(m_InputVolume);
  calculator->Compute();
  pixelT minInput = calculator->GetMinimum();
  pixelT maxInput = calculator->GetMaximum();
  pixelT rangeInput = maxInput - minInput;
  pixelT stepInput = rangeInput / 100;
  calculator->SetInput(m_InputMask);
  calculator->Compute();
  pixelT minMask = calculator->GetMinimum();
  pixelT maxMask = calculator->GetMaximum();
  pixelT rangeMask = maxMask - minMask;
  pixelT stepMask = rangeInput / 100;

  typename JointHistogramType::IndexType start;
  typename JointHistogramType::SizeType size;
  start[0] = start[1] = 0;
  size[0] = size[1] = m_NumberOfBins;
  typename JointHistogramType::RegionType region;
  region.SetIndex(start);
  region.SetSize(size);
  typename JointHistogramType::Pointer hist = JointHistogramType::New();
  hist->SetRegions(region);
  hist->Allocate();
  hist->FillBuffer(0);
  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  ConstIteratorType inputItr( m_InputVolume, 
                              m_InputVolume->GetLargestPossibleRegion() );
  ConstIteratorType maskItr( m_InputMask, 
                             m_InputMask->GetLargestPossibleRegion() );
  inputItr.GoToBegin();
  maskItr.GoToBegin();
  while( !inputItr.IsAtEnd() && !maskItr.IsAtEnd() )

    {
    typename JointHistogramType::IndexType cur;
    cur[0] = inputItr.Get() / stepInput;
    cur[1] = maskItr.Get() / stepMask;
    hist->SetPixel(cur, hist->GetPixel(cur) + 1);
    ++inputItr;
    ++maskItr;
    }

}

template< class pixelT, unsigned int dimensionT>
void
JointHistogramGenerator<pixelT,dimensionT>
::SetNumberOfBins( unsigned int numBins )
{
  m_NumberOfBins = numBins;
}
 
template< class pixelT, unsigned int dimensionT>
void
JointHistogramGenerator<pixelT,dimensionT>
::SetInputVolume( typename ImageType::Pointer inputVolume )
{
  m_InputVolume = inputVolume;
}

template< class pixelT, unsigned int dimensionT>
void
JointHistogramGenerator<pixelT,dimensionT>
::SetInputMask( typename ImageType::Pointer inputMask )
{
  m_InputMask = inputMask;
}

template< class pixelT, unsigned int dimensionT>
typename JointHistogramGenerator<pixelT,dimensionT>::ImageType::Pointer
JointHistogramGenerator<pixelT,dimensionT>
::GetOutputVolume()
{
  return m_OutputVolume;
}

} // end namespace tube
