/*=========================================================================

Library:   TubeTKLib

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

#ifndef __itktubeContrastCostFunction_hxx
#define __itktubeContrastCostFunction_hxx


namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template <class TPixel, unsigned int Dimension>
ContrastCostFunction<TPixel, Dimension>::ContrastCostFunction(void)
{
  m_InputMean = 0.0;
  m_MaskObjectValue = 0;
  m_MaskBackgroundValue = 0;
  m_CallsToGetValue = 0;
}

template <class TPixel, unsigned int Dimension>
void
ContrastCostFunction<TPixel, Dimension>::SetScales(ParametersType & scales)
{
  m_Scales = scales;
}

template <class TPixel, unsigned int Dimension>
unsigned int
ContrastCostFunction<TPixel, Dimension>::GetNumberOfParameters(void) const
{
  return 3;
}

template <class TPixel, unsigned int Dimension>
double
ContrastCostFunction<TPixel, Dimension>::GetValue(const ParametersType & params) const
{
  typename BlurFilterType::Pointer filterInputObj = BlurFilterType::New();
  filterInputObj->SetInput(m_InputImage);
  double sigmaObj = params[0];
  if (sigmaObj > 0.3 && sigmaObj < 100)
  {
    filterInputObj->SetSigma(sigmaObj);
  }
  else
  {
    return 100;
  }
  filterInputObj->Update();
  typename ImageType::Pointer imgObj = filterInputObj->GetOutput();

  typename BlurFilterType::Pointer filterInputBkg = BlurFilterType::New();
  filterInputBkg->SetInput(m_InputImage);
  double sigmaBkg = params[1];
  if (sigmaBkg > sigmaObj && sigmaBkg < 100)
  {
    filterInputBkg->SetSigma(sigmaBkg);
  }
  else
  {
    return 100;
  }
  filterInputBkg->Update();
  typename ImageType::Pointer imgBkg = filterInputBkg->GetOutput();

  double alpha = params[2];

  double meanObj = 0;
  double stdDevObj = 0;
  double countObj = 0;
  double meanBkg = 0;
  double stdDevBkg = 0;
  double countBkg = 0;

  double sumObj = 0;
  double sumsObj = 0;
  double sumBkg = 0;
  double sumsBkg = 0;

  typedef ImageRegionIterator<ImageType>      ImageIteratorType;
  typedef ImageRegionConstIterator<ImageType> ConstImageIteratorType;

  ConstImageIteratorType iterObj(imgObj, imgObj->GetLargestPossibleRegion());
  ConstImageIteratorType iterBkg(imgBkg, imgBkg->GetLargestPossibleRegion());
  ConstImageIteratorType iterMask(m_InputMask, m_InputMask->GetLargestPossibleRegion());
  ImageIteratorType      iterOut(m_OutputImage, m_OutputImage->GetLargestPossibleRegion());

  double meanRawBkg = 0;
  double countRawBkg = 0;
  while (!iterBkg.IsAtEnd())
  {
    meanRawBkg += iterBkg.Get();
    ++countRawBkg;
    ++iterBkg;
  }
  meanRawBkg /= countRawBkg;

  iterBkg.GoToBegin();
  while (!iterObj.IsAtEnd())
  {
    double tf = iterObj.Get() * (1 + alpha * (iterBkg.Get() - meanRawBkg));
    if (iterMask.Get() == m_MaskObjectValue)
    {
      sumObj += tf;
      sumsObj += tf * tf;
      ++countObj;
    }
    else if (iterMask.Get() == m_MaskBackgroundValue)
    {
      sumBkg += tf;
      sumsBkg += tf * tf;
      ++countBkg;
    }
    iterOut.Set(tf);
    ++iterObj;
    ++iterBkg;
    ++iterMask;
    ++iterOut;
  }

  if (countObj > 0)
  {
    meanObj = sumObj / countObj;
    stdDevObj = std::sqrt(sumsObj / countObj - meanObj * meanObj);
  }
  if (countBkg > 0)
  {
    meanBkg = sumBkg / countBkg;
    stdDevBkg = std::sqrt(sumsBkg / countBkg - stdDevBkg * stdDevBkg);
  }

  double dp = std::fabs(meanObj - meanBkg) / (stdDevObj * stdDevBkg);

  std::cout << ++m_CallsToGetValue << " : " << params[0] << ", " << params[1] << ", " << params[2] << ": " << meanObj
            << " ( " << stdDevObj << " ) " << meanBkg << " ( " << stdDevBkg << " ) "
            << " : " << dp << std::endl;

  return dp;
}

template <class TPixel, unsigned int Dimension>
void
ContrastCostFunction<TPixel, Dimension>::GetDerivative(const ParametersType & params, DerivativeType & deriv) const
{
  ParametersType tmpP = params;
  deriv = params;

  for (unsigned int i = 0; i < this->GetNumberOfParameters(); i++)
  {
    tmpP[i] = params[i] - 0.5 / m_Scales[i];
    double tf = this->GetValue(tmpP);
    tmpP[i] = params[i] + 0.5 / m_Scales[i];
    deriv[i] = this->GetValue(tmpP) - tf;
    tmpP[i] = params[i];
  }
}

template <class TPixel, unsigned int Dimension>
void
ContrastCostFunction<TPixel, Dimension>::Initialize(void)
{
  m_CallsToGetValue = 0;
}

} // End namespace tube

} // End namespace itk

#endif
