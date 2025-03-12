/*=========================================================================
   Library:   TubeTK

   Copyright Kitware Inc.

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

#ifndef itktubeComputeTrainingMaskFilter_h
#define itktubeComputeTrainingMaskFilter_h

#include <itkImageToImageFilter.h>

#include <itkBinaryThinningImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>

namespace itk
{

namespace tube
{

/**
 * This class returns expert vessel and not vessel mask.
 *
 * \sa ComputeTrainingMaskFilter
 */

template <class TInputImage, class TLabelMap = Image<typename TInputImage::PixelType, TInputImage::ImageDimension>>
class ComputeTrainingMaskFilter : public ImageToImageFilter<TInputImage, TLabelMap>
{
public:
  using Self = ComputeTrainingMaskFilter;
  using Superclass = ImageToImageFilter<TInputImage, TLabelMap>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;
  using ImageType = TInputImage;
  using LabelMapType = TLabelMap;

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  itkNewMacro(Self);
  const LabelMapType *
  GetObjectMask();
  const LabelMapType *
  GetNotObjectMask();
  itkSetMacro(Gap, double);
  itkSetMacro(ObjectWidth, double);
  itkSetMacro(NotObjectWidth, double);
  itkGetMacro(Gap, double);
  itkGetMacro(ObjectWidth, double);
  itkGetMacro(NotObjectWidth, double);

protected:
  ComputeTrainingMaskFilter();
  virtual ~ComputeTrainingMaskFilter();

  virtual void
  GenerateData() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  using BallType = itk::BinaryBallStructuringElement<short, ImageType::ImageDimension>;
  using DilateFilterType = itk::DilateObjectMorphologyImageFilter<ImageType, ImageType, BallType>;
  using BinaryThinningFilterType = itk::BinaryThinningImageFilter<ImageType, ImageType>;
  using ThresholdFilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
  using SubtractFilterType = itk::SubtractImageFilter<ImageType, ImageType, ImageType>;
  using MultiplyFilterType = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>;
  using AddFilterType = itk::AddImageFilter<ImageType, ImageType, ImageType>;
  using CastFilterType = itk::CastImageFilter<ImageType, LabelMapType>;

  ComputeTrainingMaskFilter(const Self &);
  void
  operator=(const Self &);
  void
  ApplyDilateMorphologyFilter(typename ImageType::Pointer & input, int size = 1);

  typename AddFilterType::Pointer            m_Add;
  typename MultiplyFilterType::Pointer       m_Multiply;
  typename ThresholdFilterType::Pointer      m_Threshold;
  typename BinaryThinningFilterType::Pointer m_BinaryThinning;
  typename DilateFilterType::Pointer         m_Dilate;
  typename SubtractFilterType::Pointer       m_Subtract;
  typename MultiplyFilterType::Pointer       m_MultiplyCenterLine;
  typename MultiplyFilterType::Pointer       m_MultiplyOutside;
  typename CastFilterType::Pointer           m_Cast;
  typename CastFilterType::Pointer           m_CastObject;
  typename CastFilterType::Pointer           m_CastNotObject;

  BallType m_Ball;
  double   m_Gap;
  double   m_ObjectWidth;
  double   m_NotObjectWidth;
};

} // namespace tube
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeComputeTrainingMaskFilter.hxx"
#endif

#endif
