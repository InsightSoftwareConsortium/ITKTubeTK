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
#ifndef itktubeComputeImageSimilarityMetrics_hxx
#define itktubeComputeImageSimilarityMetrics_hxx

// ITK includes
#include <itkIdentityTransform.h>
#include <itkMutualInformationImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkNormalizeImageFilter.h>

// TubeTK includes

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template <class TInputImage>
ComputeImageSimilarityMetrics<TInputImage>::ComputeImageSimilarityMetrics()
{
  m_Input1 = NULL;
  m_Input2 = NULL;
  m_SamplingRate = 0.05;
  m_UseCorrelation = false;
}

template <class TInputImage>
void
ComputeImageSimilarityMetrics<TInputImage>::Update(void)
{
  // check if both input images are set
  if (m_Input1.IsNull())
  {
    itkExceptionMacro("Input Image 1 is not set");
  }

  if (m_Input1.IsNull())
  {
    itkExceptionMacro("Input Image 2 is not set");
  }

  // Normalize the images
  using NormFilterType = itk::NormalizeImageFilter<ImageType, ImageType>;

  typename NormFilterType::Pointer norm1 = NormFilterType::New();
  norm1->SetInput(m_Input1);
  norm1->Update();

  typename NormFilterType::Pointer norm2 = NormFilterType::New();
  norm2->SetInput(m_Input2);
  norm2->Update();

  // Compute similarity
  using MetricType = itk::ImageToImageMetric<ImageType, ImageType>;
  typename MetricType::Pointer metric;

  using TransformType = itk::IdentityTransform<double, TInputImage::ImageDimension>;
  typename TransformType::Pointer transform = TransformType::New();

  using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(norm2->GetOutput());

  if (!m_UseCorrelation)
  {
    using MIMetricType = itk::MutualInformationImageToImageMetric<ImageType, ImageType>;
    metric = MIMetricType::New();
  }
  else
  {
    using CorMetricType = itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType>;
    metric = CorMetricType::New();
  }

  typename ImageType::SizeType size = norm1->GetOutput()->GetLargestPossibleRegion().GetSize();

  metric->SetFixedImage(norm1->GetOutput());
  metric->SetMovingImage(norm2->GetOutput());
  metric->SetFixedImageRegion(norm1->GetOutput()->GetLargestPossibleRegion());
  metric->SetTransform(transform);
  metric->SetInterpolator(interpolator);
  metric->SetNumberOfSpatialSamples(size[0] * size[1] * m_SamplingRate);
  metric->Initialize();
  metric->MultiThreadingInitialize();

  if (!m_UseCorrelation)
  {
    m_Output = metric->GetValue(transform->GetParameters());
  }
  else
  {
    m_Output = -metric->GetValue(transform->GetParameters());
  }
}

template <class TInputImage>
void
ComputeImageSimilarityMetrics<TInputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << "Use Correlation: " << m_UseCorrelation << std::endl;
  os << "Sampling Rate: " << m_SamplingRate << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif
