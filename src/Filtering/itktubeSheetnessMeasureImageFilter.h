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

#ifndef itktubeSheetnessMeasureImageFilter_h
#define itktubeSheetnessMeasureImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{

namespace tube
{

/** \class SheetnessMeasureImageFilter
 *
 * \brief Computes a measure of sheetness from the Hessian Eigenvalues
 *
 * Based on the "Sheetness" measure proposed by Decouteaux et. al.
 *
 * M.Descoteaux, M.Audette, K.Chinzei, el al.:
 * "Bone enhancement filtering: Application to sinus bone segmentation
 *  and simulation of pituitary surgery."
 *

 * \sa HessianRecursiveGaussianImageFilter
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 *
 * \ingroup IntensityImageFilters TensorObjects
 *
 * \ingroup ITKImageFeature
 */

template <class TPixel>
class SheetnessMeasureImageFilter
  : public ImageToImageFilter<Image<SymmetricSecondRankTensor<double, 3>, 3>, Image<TPixel, 3>>
{
public:
  /** Standard class type alias. */
  using Self = SheetnessMeasureImageFilter;
  using Superclass = ImageToImageFilter<Image<SymmetricSecondRankTensor<double, 3>, 3>, Image<TPixel, 3>>;

  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using InputImageType = typename Superclass::InputImageType;
  using OutputImageType = typename Superclass::OutputImageType;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = TPixel;

  /** Image dimension = 3. */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);
  itkStaticConstMacro(InputPixelDimension, unsigned int, InputPixelType::Dimension);

  using EigenValueArrayType = FixedArray<double, itkGetStaticConstMacro(InputPixelDimension)>;
  using EigenValueImageType = Image<EigenValueArrayType, itkGetStaticConstMacro(ImageDimension)>;
  using EigenAnalysisFilterType = SymmetricEigenAnalysisImageFilter<InputImageType, EigenValueImageType>;

  /** Run-time type information ( and related methods ).   */
  itkOverrideGetNameOfClassMacro(SheetnessMeasureImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get macros for alpha */
  itkSetMacro(Alpha, double);
  itkGetConstMacro(Alpha, double);

  /** Set/Get macros for Beta. */
  itkSetMacro(Beta, double);
  itkGetConstMacro(Beta, double);

  /** Set/Get macros for Cfactor. */
  itkSetMacro(Cfactor, double);
  itkGetConstMacro(Cfactor, double);

  /** Set/Get DetectBrightSheets */
  itkBooleanMacro(DetectBrightSheets);
  itkSetMacro(DetectBrightSheets, bool);


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(DoubleConvertibleToOutputCheck, (Concept::Convertible<double, OutputPixelType>));
  /** End concept checking */
#endif
protected:
  SheetnessMeasureImageFilter(void);
  ~SheetnessMeasureImageFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Generate Data */
  void
  GenerateData(void) override;

private:
  SheetnessMeasureImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  typename EigenAnalysisFilterType::Pointer m_SymmetricEigenValueFilter;

  double m_Alpha;
  double m_Beta;
  double m_Cfactor;
  bool   m_DetectBrightSheets;

}; // End class SheetnessMeasureImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeSheetnessMeasureImageFilter.hxx"
#endif

#endif // End !defined( __itktubeSheetnessMeasureImageFilter_h )
