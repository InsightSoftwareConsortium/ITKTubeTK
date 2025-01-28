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

#ifndef __itktubeDifferenceImageFilter_hxx
#define __itktubeDifferenceImageFilter_hxx


#include <itkImageRegionIterator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkProgressReporter.h>
#include <itkZeroFluxNeumannBoundaryCondition.h>

#include <cmath>
#include <math.h>

namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
DifferenceImageFilter<TInputImage, TOutputImage>::DifferenceImageFilter(void)
{
  // We require two inputs to execute.
  this->SetNumberOfRequiredInputs(2);

  // Set the default DifferenceThreshold.
  m_DifferenceThreshold = NumericTraits<OutputPixelType>::Zero;

  // Set the default ToleranceRadius.
  m_ToleranceRadius = 0;

  // Initialize statistics about difference image.
  m_MeanDifference = NumericTraits<RealType>::Zero;
  m_TotalDifference = NumericTraits<AccumulateType>::Zero;
  m_NumberOfPixelsWithDifferences = 0;
  m_IgnoreBoundaryPixels = false;
  // Use the ITKv4 Threading Model (call ThreadedGenerateData instead of DynamicThreadedGenerateData)
  this->DynamicMultiThreadingOff();
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "ToleranceRadius: " << m_ToleranceRadius << "\n";
  os << indent << "DifferenceThreshold: " << m_DifferenceThreshold << "\n";
  os << indent << "MeanDifference: " << m_MeanDifference << "\n";
  os << indent << "TotalDifference: " << m_TotalDifference << "\n";
  os << indent << "NumberOfPixelsWithDifferences: " << m_NumberOfPixelsWithDifferences << "\n";
  os << indent << "IgnoreBoundaryPixels: " << m_IgnoreBoundaryPixels << "\n";
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter<TInputImage, TOutputImage>::SetValidInput(const InputImageType * validImage)
{
  // The valid image should be input 0.
  this->SetInput(0, validImage);
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter<TInputImage, TOutputImage>::SetTestInput(const InputImageType * testImage)
{
  // The test image should be input 1.
  this->SetInput(1, testImage);
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData(void)
{
  int numberOfThreads = this->GetNumberOfWorkUnits();

  // Initialize statistics about difference image.
  m_MeanDifference = NumericTraits<RealType>::Zero;
  m_TotalDifference = NumericTraits<AccumulateType>::Zero;
  m_NumberOfPixelsWithDifferences = 0;

  // Resize the thread temporaries
  m_ThreadDifferenceSum.SetSize(numberOfThreads);
  m_ThreadNumberOfPixelsWithDifferences.SetSize(numberOfThreads);

  // Initialize the temporaries
  m_ThreadDifferenceSum.Fill(NumericTraits<AccumulateType>::Zero);
  m_ThreadNumberOfPixelsWithDifferences.Fill(0);
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(const OutputImageRegionType & threadRegion,
                                                                       ThreadIdType                  threadId)
{
  typedef ConstNeighborhoodIterator<InputImageType>                           SmartIterator;
  typedef ImageRegionConstIterator<InputImageType>                            InputIterator;
  typedef ImageRegionIterator<OutputImageType>                                OutputIterator;
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> FacesCalculator;
  typedef typename FacesCalculator::RadiusType                                RadiusType;
  typedef typename FacesCalculator::FaceListType                              FaceListType;
  typedef typename FaceListType::iterator                                     FaceListIterator;
  typedef typename InputImageType::PixelType                                  InputPixelType;

  // Prepare standard boundary condition.
  ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;

  // Get a pointer to each image.
  const InputImageType * validImage = this->GetInput(0);
  const InputImageType * testImage = this->GetInput(1);
  OutputImageType *      outputPtr = this->GetOutput();

  // Create a radius of pixels.
  RadiusType radius;
  if (m_ToleranceRadius > 0)
  {
    radius.Fill(m_ToleranceRadius);
  }
  else
  {
    radius.Fill(0);
  }

  // Find the data-set boundary faces.
  FacesCalculator boundaryCalculator;
  FaceListType    faceList = boundaryCalculator(testImage, threadRegion, radius);

  // Support progress methods/callbacks.
  ProgressReporter progress(this, threadId, threadRegion.GetNumberOfPixels());

  // Process the internal face and each of the boundary faces.
  for (FaceListIterator face = faceList.begin(); face != faceList.end(); ++face)
  {
    SmartIterator test(radius, testImage, *face);
    // Iterate over test image.

    InputIterator valid(validImage, *face);
    // Iterate over valid image.

    OutputIterator out(outputPtr, *face);
    // Iterate over output image.

    if (!test.GetNeedToUseBoundaryCondition() || !m_IgnoreBoundaryPixels)
    {
      test.OverrideBoundaryCondition(&nbc);

      for (valid.GoToBegin(), test.GoToBegin(), out.GoToBegin(); !valid.IsAtEnd(); ++valid, ++test, ++out)
      {
        // Get the current valid pixel.
        InputPixelType t = valid.Get();

        //  Assume a good match - so test center pixel first, for speed
        RealType difference = static_cast<RealType>(t) - test.GetCenterPixel();
        if (std::isnan(difference))
        {
          difference = m_DifferenceThreshold * 2;
        }
        RealType absDifference = difference;
        if (NumericTraits<RealType>::IsNegative(difference))
        {
          absDifference = -absDifference;
        }
        RealType minimumDifference = difference;
        RealType minimumAbsDifference = absDifference;

        unsigned int neighborhoodSize = test.Size();
        for (unsigned int i = 0; i < neighborhoodSize; ++i)
        {
          // Use the RealType for the difference to make sure we get the
          // sign.
          difference = static_cast<RealType>(t) - test.GetPixel(i);
          if (std::isnan(difference))
          {
            difference = m_DifferenceThreshold * 2;
          }
          absDifference = difference;
          if (NumericTraits<RealType>::IsNegative(difference))
          {
            absDifference = -absDifference;
          }
          if (absDifference < minimumAbsDifference)
          {
            minimumDifference = difference;
            minimumAbsDifference = absDifference;
          }
        }

        // Store the minimum difference value in the output image.
        out.Set(static_cast<OutputPixelType>(minimumDifference));

        m_ThreadDifferenceSum[threadId] += minimumAbsDifference;
        if (minimumAbsDifference > m_DifferenceThreshold)
        {
          m_ThreadNumberOfPixelsWithDifferences[threadId]++;
        }

        // Update progress.
        progress.CompletedPixel();
      }
    }
    else
    {
      for (out.GoToBegin(); !out.IsAtEnd(); ++out)
      {
        out.Set(NumericTraits<OutputPixelType>::Zero);
        progress.CompletedPixel();
      }
    }
  }
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter<TInputImage, TOutputImage>::AfterThreadedGenerateData(void)
{
  // Set statistics about difference image.
  int numberOfThreads = this->GetNumberOfWorkUnits();
  for (int i = 0; i < numberOfThreads; ++i)
  {
    m_TotalDifference += m_ThreadDifferenceSum[i];
    m_NumberOfPixelsWithDifferences += m_ThreadNumberOfPixelsWithDifferences[i];
  }

  // Get the total number of pixels processed in the region.
  // This is different from the m_TotalNumberOfPixels which
  // is the number of pixels that actually have differences
  // above the intensity threshold.
  OutputImageRegionType region = this->GetOutput()->GetRequestedRegion();
  AccumulateType        numberOfPixels = region.GetNumberOfPixels();

  // Calculate the mean difference.
  m_MeanDifference = m_TotalDifference / numberOfPixels;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeDifferenceImageFilter_hxx )
