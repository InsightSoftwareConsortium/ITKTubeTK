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

#ifndef itktubeSingleValuedCostFunctionImageSource_h
#define itktubeSingleValuedCostFunctionImageSource_h

#include <itkImageSource.h>

namespace itk
{

namespace tube
{

/**
 * \class SingleValuedCostFunctionImageSource
 *
 * \brief Create an image by evaluating a single valued cost function on a
 * uniform grid evaluated over parameter ranges.
 *
 */
template <class TCostFunction, unsigned int VNumberOfParameters>
class SingleValuedCostFunctionImageSource
  : public ImageSource<Image<typename TCostFunction::MeasureType, VNumberOfParameters>>
{
public:
  using CostFunctionImageType = Image<typename TCostFunction::MeasureType, VNumberOfParameters>;

  using Self = SingleValuedCostFunctionImageSource;
  using Superclass = ImageSource<CostFunctionImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using CostFunctionType = TCostFunction;
  using MeasureType = typename CostFunctionType::MeasureType;
  using ParametersType = typename CostFunctionType::ParametersType;

  using OutputImageType = typename Superclass::OutputImageType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  /** Dimensionality of the output image. */
  itkStaticConstMacro(NumberOfParameters, unsigned int, VNumberOfParameters);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(SingleValuedCostFunctionImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set the single valued cost function from which an Image will be created by
   * by evaluating over the range of parameters given.  The cost function
   * already be initialized with all its inputs, etc. */
  itkSetObjectMacro(CostFunction, CostFunctionType);
  itkGetConstObjectMacro(CostFunction, CostFunctionType);

  /** Set the lower bound of parameter values to examine.  This determines the
   * Origin of the output Image. */
  virtual void
  SetParametersLowerBound(const ParametersType & lowerBound);
  ParametersType
  GetParametersLowerBound() const;

  /** Set the upper bound of parameter values to examine. */
  virtual void
  SetParametersUpperBound(const ParametersType & upperBound);
  ParametersType
  GetParametersUpperBound() const;

  /** Set the step in parameter values to examine. This parameters are evaluated
   * over the range given by the ParametersLowerBound, ParametersUpperBound at
   * equal intervals determined by the ParametersStep.  The ParametersStep
   * defines the Spacing in the output Image. */
  virtual void
  SetParametersStep(const ParametersType & step);
  ParametersType
  GetParametersStep() const;

protected:
  SingleValuedCostFunctionImageSource();
  virtual ~SingleValuedCostFunctionImageSource() {}

  virtual void
  GenerateOutputInformation() override;

  virtual void
  BeforeThreadedGenerateData() override;

  virtual void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

private:
  SingleValuedCostFunctionImageSource(const SingleValuedCostFunctionImageSource &); // purposely not implemented
  void
  operator=(const SingleValuedCostFunctionImageSource &); // purposely not implemented

  typename CostFunctionType::Pointer m_CostFunction;

  ParametersType m_ParametersLowerBound;
  ParametersType m_ParametersUpperBound;
  ParametersType m_ParametersStep;
}; // End class SingleValuedCostFunctionImageSource

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeSingleValuedCostFunctionImageSource.hxx"
#endif

#endif // !defined( __itktubeSingleValuedCostFunctionImageSource_h )
