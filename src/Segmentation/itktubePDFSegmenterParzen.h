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

#ifndef itktubePDFSegmenterParzen_h
#define itktubePDFSegmenterParzen_h

#include "itktubePDFSegmenterBase.h"

#include <itkImage.h>
#include <itkListSample.h>

#include <vector>

namespace itk
{

namespace tube
{

#define PARZEN_MAX_NUMBER_OF_FEATURES 4

template <class TImage, class TLabelMap>
class PDFSegmenterParzen : public PDFSegmenterBase<TImage, TLabelMap>
{
public:
  using Self = PDFSegmenterParzen;
  using Superclass = PDFSegmenterBase<TImage, TLabelMap>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkTypeMacro(PDFSegmenterParzen, PDFSegmenterBase);

  itkNewMacro(Self);

  //
  // Template Args Typedefs
  //
  using InputImageType = TImage;
  using LabelMapType = TLabelMap;

  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  //
  // Superclass Typedefs
  //
  typedef typename Superclass::FeatureVectorGeneratorType FeatureVectorGeneratorType;
  using FeatureValueType = typename Superclass::FeatureValueType;
  using FeatureVectorType = typename Superclass::FeatureVectorType;
  using FeatureImageType = typename Superclass::FeatureImageType;

  using LabelMapPixelType = typename Superclass::LabelMapPixelType;

  using ObjectIdType = typename Superclass::ObjectIdType;
  using ObjectIdListType = typename Superclass::ObjectIdListType;

  typedef typename Superclass::ProbabilityPixelType ProbabilityPixelType;

  typedef typename Superclass::ProbabilityVectorType ProbabilityVectorType;
  typedef typename Superclass::ProbabilityImageType  ProbabilityImageType;

  using VectorDoubleType = typename Superclass::VectorDoubleType;
  using VectorIntType = typename Superclass::VectorIntType;
  using VectorUIntType = typename Superclass::VectorUIntType;

  //
  // Custom Typedefs
  //
  using HistogramPixelType = float;

  using HistogramImageType = Image<HistogramPixelType, PARZEN_MAX_NUMBER_OF_FEATURES>;

  using PDFPixelType = HistogramPixelType;
  using PDFImageType = HistogramImageType;

  using LabeledFeatureSpaceType = Image<LabelMapPixelType, PARZEN_MAX_NUMBER_OF_FEATURES>;

  //
  // Methods
  //
  itkSetMacro(HistogramSmoothingStandardDeviation, double);
  itkGetMacro(HistogramSmoothingStandardDeviation, double);
  itkSetMacro(OutlierRejectPortion, double);
  itkGetMacro(OutlierRejectPortion, double);

  typename PDFImageType::Pointer
  GetClassPDFImage(unsigned int classNum) const;

  void
  SetClassPDFImage(unsigned int classNum, PDFImageType * classPDF);

  const VectorUIntType &
  GetNumberOfBinsPerFeature(void) const;
  void
  SetNumberOfBinsPerFeature(const VectorUIntType & nBin);
  const VectorDoubleType &
  GetBinMin(void) const;
  void
  SetBinMin(const VectorDoubleType & binMin);
  const VectorDoubleType &
  GetBinSize(void) const;
  void
  SetBinSize(const VectorDoubleType & binMin);

  /** Given one PDF per class, generate a labelmap of feature space */
  void
  GenerateLabeledFeatureSpace(void);

  void
  SetLabeledFeatureSpace(LabeledFeatureSpaceType * labeledFeatureSpace);

  typename LabeledFeatureSpaceType::Pointer
  GetLabeledFeatureSpace(void) const;

  virtual void
  Update(void) override;

  //
  // Must overwrite
  //
  virtual ProbabilityVectorType
  GetProbabilityVector(const FeatureVectorType & fv) const override;

protected:
  PDFSegmenterParzen(void);
  virtual ~PDFSegmenterParzen(void);

  virtual void
  GeneratePDFs(void) override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  PDFSegmenterParzen(const Self &); // Purposely not implemented
  void
  operator=(const Self &); // Purposely not implemented

  // Superclass type alias
  using ProbabilityImageVectorType = std::vector<typename ProbabilityImageType::Pointer>;

  using ListVectorType = std::vector<ProbabilityPixelType>;
  using ListSampleType = std::vector<ListVectorType>;
  using ClassListSampleType = std::vector<ListSampleType>;

  // Custom type alias
  using ClassHistogramImageType = std::vector<typename HistogramImageType::Pointer>;

  ClassHistogramImageType m_InClassHistogram;
  VectorDoubleType        m_HistogramBinMin;
  VectorDoubleType        m_HistogramBinSize;
  VectorUIntType          m_HistogramNumberOfBin;

  double m_OutlierRejectPortion;

  double m_HistogramSmoothingStandardDeviation;

  typename LabeledFeatureSpaceType::Pointer m_LabeledFeatureSpace;

}; // End class PDFSegmenterParzen

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubePDFSegmenterParzen.hxx"
#endif

#endif // End !defined( __itktubePDFSegmenterParzen_h )
