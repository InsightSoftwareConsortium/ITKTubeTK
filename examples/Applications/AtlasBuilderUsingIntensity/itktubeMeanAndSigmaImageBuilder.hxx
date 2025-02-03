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

#ifndef __itktubeMeanAndSigmaImageBuilder_hxx
#define __itktubeMeanAndSigmaImageBuilder_hxx


template <class TInputImageType, class TOutputMeanImageType, class TOutputSigmaImageType>
MeanAndSigmaImageBuilder<TInputImageType, TOutputMeanImageType, TOutputSigmaImageType>::MeanAndSigmaImageBuilder(void)
  : m_ImageCountThreshold(1)
  , m_ThresholdInputImageBelowOn(false)
  , m_ThresholdInputImageBelow(0)
  , m_IsProcessing(false)
  , m_UseStandardDeviation(true)
  , m_DynamicallyAdjustOutputSize(false)
{
  m_OutputSize.Fill(0);
  m_OutputSpacing.Fill(0);
  m_OutputOrigin.Fill(0);
}

template <class TInputImageType, class TOutputMeanImageType, class TOutputSigmaImageType>
void
MeanAndSigmaImageBuilder<TInputImageType, TOutputMeanImageType, TOutputSigmaImageType>::AddImage(InputImagePointer i)
{
  if (!(this->GetIsProcessing()))
  {
    BuildProcessingImages(i);
  }

  // Iterate through current_image and copy values into the Sum Image,
  // Sum Sqared Image, and Valid Count Image
  //
  // Note: Must transform the base image region so that the new added
  // image has the appropriate region value
  ProcessImagePointer sumImage = this->GetSumImage();
  ProcessImagePointer sumSquareImage = this->GetSumSquareImage();
  CountImagePointer   validImage = this->GetValidCountImage();

  typename ProcessImageType::RegionType region;
  region.SetSize(this->GetOutputSize());

  InputConstIteratorType it_image(i, region);
  CountIteratorType      it_valid(validImage, region);
  ProcessIteratorType    it_sum(sumImage, region);
  ProcessIteratorType    it_sum_sqre(sumSquareImage, region);

  it_image.GoToBegin();
  it_valid.GoToBegin();
  it_sum.GoToBegin();
  it_sum_sqre.GoToBegin();

  while (!it_image.IsAtEnd() && !it_sum.IsAtEnd())
  {
    InputPixelType p = it_image.Get();

    // If requested to be threshold, ensure that the value
    // added has a value above the threshold ( and not out of the
    // image area ), else ignore.
    if (this->GetThresholdInputImageBelowOn() == false ||
        (this->GetThresholdInputImageBelowOn() == true && p > this->GetThresholdInputImageBelow()))
    {
      it_valid.Set(it_valid.Get() + 1);
      it_sum.Set(it_sum.Get() + p);
      it_sum_sqre.Set(it_sum_sqre.Get() + (p * p));
    }
    ++it_image;
    ++it_sum;
    ++it_sum_sqre;
    ++it_valid;
  }
  this->SetSumImage(sumImage);
  this->SetSumSquareImage(sumSquareImage);
  this->SetValidCountImage(validImage);
}

template <class TInputImageType, class TOutputMeanImageType, class TOutputSigmaImageType>
void
MeanAndSigmaImageBuilder<TInputImageType, TOutputMeanImageType, TOutputSigmaImageType>::FinalizeOutput(void)
{
  if (!(this->GetIsProcessing()))
  {
    ::tube::FmtErrorMessage("Must call Start() before Finalizing");
    return;
  }

  ProcessImagePointer sumImage = this->GetSumImage();
  ProcessImagePointer sumSquareImage = this->GetSumSquareImage();
  CountImagePointer   validImages = this->GetValidCountImage();

  RegionType  outputRegion = sumImage->GetLargestPossibleRegion();
  SpacingType outputSpacing = sumImage->GetSpacing();
  PointType   outputOrigin = sumImage->GetOrigin();

  // Build Mean and Variance Images
  OutputMeanImagePointer meanImage = OutputMeanImageType::New();
  meanImage->SetRegions(outputRegion);
  meanImage->SetSpacing(outputSpacing);
  meanImage->SetOrigin(outputOrigin);
  meanImage->Allocate();

  OutputSigmaImagePointer sigmaImage = OutputSigmaImageType::New();
  sigmaImage->SetRegions(outputRegion);
  sigmaImage->SetSpacing(outputSpacing);
  sigmaImage->SetOrigin(outputOrigin);
  sigmaImage->Allocate();

  CountConstIteratorType   it_valid(validImages, outputRegion);
  ProcessConstIteratorType it_sum(sumImage, outputRegion);
  ProcessConstIteratorType it_sum_sqre(sumSquareImage, outputRegion);
  OutputMeanIteratorType   it_mean(meanImage, outputRegion);
  OutputSigmaIteratorType  it_dev(sigmaImage, outputRegion);

  it_sum.GoToBegin();
  it_sum_sqre.GoToBegin();
  it_dev.GoToBegin();
  it_mean.GoToBegin();
  it_valid.GoToBegin();

  // The minimum number of contributing images to count the voxel
  // in the output mean and variance
  const unsigned short countThreshold = this->GetImageCountThreshold();

  // Calculate standard deviation or varaince
  const bool isStdDeviation = this->GetUseStandardDeviation();

  while (!it_sum.IsAtEnd())
  {
    // Ensure that the number of valid images at the point is above
    // the set minimum image threshold
    CountPixelType number = it_valid.Get();
    if (number >= countThreshold)
    {
      if (number > 1) // Prevent potential division by zero for variance
      {
        ProcessPixelType sum = it_sum.Get();
        ProcessPixelType sumSqr = it_sum_sqre.Get();
        // Variance Calc. s^2 = ( ( sum( x^2 ) - ( sum( x )^2/n ) )/ ( n-1 ) )
        ProcessPixelType variance = (sumSqr - ((sum * sum) / number)) / (number - 1);

        // If Standard Deviation Calc. s = std::sqrt( s^2 )
        if (isStdDeviation)
        {
          variance = std::sqrt(variance);
        }
        it_dev.Set((OutputSigmaPixelType)variance);
        it_mean.Set((OutputMeanPixelType)(sum / ((ProcessPixelType)number)));
      }
      else
      {
        it_dev.Set(0);
        it_mean.Set((OutputMeanPixelType)it_sum.Get());
      }
    }
    else
    {
      it_dev.Set(0);
      it_mean.Set(0);
    }

    ++it_valid;
    ++it_mean;
    ++it_dev;
    ++it_sum;
    ++it_sum_sqre;
  }

  this->SetOutputMeanImage(meanImage);
  this->SetOutputSigmaImage(sigmaImage);

  // Unset the processing flag and dump the processing image pointers
  this->SetIsProcessing(false);
}

template <class TInputImageType, class TOutputMeanImageType, class TOutputSigmaImageType>
void
MeanAndSigmaImageBuilder<TInputImageType, TOutputMeanImageType, TOutputSigmaImageType>::BuildProcessingImages(
  InputImagePointer i)
{
  ProcessImagePointer sumImage = ProcessImageType::New();
  ProcessImagePointer sumSquareImage = ProcessImageType::New();
  CountImagePointer   validCountImage = CountImageType::New();

  RegionType  region = i->GetLargestPossibleRegion();
  SpacingType spacing = i->GetSpacing();
  PointType   origin = i->GetOrigin();

  sumImage->SetRegions(region);
  sumImage->SetSpacing(spacing);
  sumImage->SetOrigin(origin);
  sumImage->Allocate();
  sumImage->FillBuffer(0);

  sumSquareImage->SetRegions(region);
  sumSquareImage->SetSpacing(spacing);
  sumSquareImage->SetOrigin(origin);
  sumSquareImage->Allocate();
  sumSquareImage->FillBuffer(0);

  validCountImage->SetRegions(region);
  validCountImage->SetSpacing(spacing);
  validCountImage->SetOrigin(origin);
  validCountImage->Allocate();
  validCountImage->FillBuffer(0);

  this->SetSumImage(sumImage);
  this->SetSumSquareImage(sumSquareImage);
  this->SetValidCountImage(validCountImage);

  // Set the new output parameters
  this->SetOutputSize(region.GetSize());
  this->SetOutputSpacing(spacing);
  this->SetOutputOrigin(origin);

  this->SetIsProcessing(true);
}

template <class TInputImageType, class TOutputMeanImageType, class TOutputSigmaImageType>
void
MeanAndSigmaImageBuilder<TInputImageType, TOutputMeanImageType, TOutputSigmaImageType>::UpdateOutputImageSize(
  SizeType inputSize)
{
  if (!(this->GetIsProcessing()))
  {
    ::tube::ErrorMessage("Call AddImage() before updating image size!");
    return;
  }

  ProcessImagePointer sumImage = this->GetSumImage();
  ProcessImagePointer sumSquareImage = this->GetSumSquareImage();
  CountImagePointer   validImage = this->GetValidCountImage();

  typedef ResampleImageFilter<ProcessImageType, ProcessImageType> ResampleProcessImageType;
  typename ResampleProcessImageType::Pointer                      processFilter = ResampleProcessImageType::New();

  processFilter->SetInput(sumImage);
  processFilter->SetSize(inputSize);
  processFilter->SetOutputSpacing(sumImage->GetSpacing());
  processFilter->SetOutputOrigin(sumImage->GetOrigin());
  processFilter->Update();

  this->SetSumImage(processFilter->GetOutput());

  typename ResampleProcessImageType::Pointer processFilter2 = ResampleProcessImageType::New();
  processFilter2->SetInput(sumSquareImage);
  processFilter2->SetSize(inputSize);
  // Keep all the spacing an origins constant for all images
  // ( so using only one to declare )
  processFilter2->SetOutputSpacing(sumImage->GetSpacing());
  // Keep all the spacing an origins constant for all images
  // ( so using only one to declare )
  processFilter2->SetOutputOrigin(sumImage->GetOrigin());
  processFilter2->Update();
  this->SetSumSquareImage(processFilter2->GetOutput());


  typedef ResampleImageFilter<CountImageType, CountImageType> ResampleCountImageType;
  typename ResampleCountImageType::Pointer                    countFilter = ResampleCountImageType::New();

  countFilter->SetInput(validImage);
  countFilter->SetSize(inputSize);

  // Keep all the spacing an origins constant for all images ( so using only
  // one to declare )
  countFilter->SetOutputSpacing(sumImage->GetSpacing());
  // Keep all the spacing an origins constant for all images ( so using only
  // one to declare )
  countFilter->SetOutputOrigin(sumImage->GetOrigin());
  countFilter->Update();
  this->SetValidCountImage(countFilter->GetOutput());

  this->SetOutputSize(inputSize);
}

#endif // End !defined( __itktubeMeanAndSigmaImageBuilder_hxx )
