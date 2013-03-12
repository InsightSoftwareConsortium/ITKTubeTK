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
#ifndef __itkTubeNJetImageFunction_txx
#define __itkTubeNJetImageFunction_txx

#include "itkTubeNJetImageFunction.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vnl/algo/vnl_symmetric_eigensystem.h"


namespace itk
{

namespace tube
{

/**
 * Set the input Image
 */
template <class TInputImage>
NJetImageFunction<TInputImage>
::NJetImageFunction()
{
  m_InputImage = 0;
  m_InputImageMask = 0;

  m_UseInputImageMask = false;

  m_Extent = 3;

  m_ValidStats = false;
  m_InputImageMinX.Fill(0);
  m_InputImageMaxX.Fill(0);
  m_InputImageSize.Fill(0);
  m_InputImageSpacing.Fill(1);
  m_InputImageSpacingSquared.Fill(1);
  m_UseProjection = true;

  m_MostRecentIntensity = 0;
  m_MostRecentDerivative.Fill(0);
  m_MostRecentHessian.Fill(0);
  m_MostRecentRidgeness = 0;
  m_MostRecentRidgeRoundness = 0;
  m_MostRecentRidgeLevelness = 0;
  m_MostRecentRidgeCurvature = 0;
  m_MostRecentRidgeTangent.Fill(0);
}

/**
 * Set the input Image
 */
template <class TInputImage>
void
NJetImageFunction<TInputImage>
::SetInputImage( const InputImageType * ptr )
{
  m_InputImage = ptr;

  if(ptr != 0)
    {
    m_InputImageMinX = m_InputImage->GetLargestPossibleRegion().GetIndex();
    m_InputImageSize = m_InputImage->GetLargestPossibleRegion().GetSize();
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      m_InputImageMaxX[i] = m_InputImageMinX[i] + m_InputImageSize[i] - 1;
      }
    m_InputImageSpacing  = m_InputImage->GetSpacing();

    int i;
    for(i=0; i<ImageDimension; i++)
      {
      m_InputImageSpacingSquared[i] = m_InputImageSpacing[i]
                                      * m_InputImageSpacing[i];
      }
    }

  /* Values by default */
  this->SetExtent(3);

  m_UseInputImageMask = false;
  m_ValidStats = false;
}

/**
 * Set the input Image
 */
template <class TInputImage>
void
NJetImageFunction<TInputImage>
::SetInputImageMask( const InputImageType * ptr )
{
  m_InputImageMask = ptr;

  if(ptr != 0)
    {
    for(int i=0; i<ImageDimension; i++)
      {
      if(ptr->GetLargestPossibleRegion().GetSize()[i] != m_InputImageSize[i])
        {
        std::cout << "NJetImageFunction: ImageSize and MaskSize do not match!"
                  << std::endl
                  << "   Warning - InputImageMask not used!" << std::endl;
        return;
        }
      }

    m_UseInputImageMask = true;
    }
  else
    {
    m_UseInputImageMask = false;
    }

  m_ValidStats = false;
}

/**
 * Print
 */
template<class TInputImage>
void
NJetImageFunction<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "m_UseProjection = " << m_UseProjection << std::endl;
  os << indent << "m_UseInputImageMask = " << m_UseInputImageMask
     << std::endl;
  if( m_UseInputImageMask )
    {
    os << indent << "m_InputImageMask = " << m_InputImageMask << std::endl;
    }
  else
    {
    os << indent << "m_InputImageMask = NULL" << std::endl;
    }

  os << indent << "m_InputImageMinX = " << m_InputImageMinX
     << std::endl;
  os << indent << "m_InputImageMaxX = " << m_InputImageMaxX
     << std::endl;
  os << indent << "m_InputImageSize = " << m_InputImageSize
     << std::endl;
  os << indent << "m_InputImageSpacing = " << m_InputImageSpacing
     << std::endl;
  os << indent << "m_InputImageSpacingSquared = "
     << m_InputImageSpacingSquared << std::endl;
  if( m_InputImage.IsNotNull() )
    {
    os << indent << "m_InputImage = " << m_InputImage << std::endl;
    }
  else
    {
    os << indent << "m_InputImage = NULL" << std::endl;
    }
  os << indent << "m_Extent = " << m_Extent << std::endl;
  os << indent << "m_MostRecentIntensity = " << m_MostRecentIntensity
    << std::endl;
  os << indent << "m_MostRecentDerivative = " << m_MostRecentDerivative
    << std::endl;
  os << indent << "m_MostRecentHessian = " << m_MostRecentHessian
    << std::endl;
  os << indent << "m_MostRecentRidgeness = " << m_MostRecentRidgeness
    << std::endl;
}


/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
void
NJetImageFunction<TInputImage>
::ComputeStatistics(void)
{
  if(m_InputImage)
    {
    m_ValidStats = true;
    ImageRegionConstIterator<InputImageType> imageIt
    (m_InputImage, m_InputImage->GetLargestPossibleRegion());

    if(m_UseInputImageMask)
      {
      ImageRegionConstIterator<InputImageType>
      maskIt(m_InputImageMask,
                 m_InputImageMask->GetLargestPossibleRegion());

      imageIt.GoToBegin();
      maskIt.GoToBegin();
      while(!maskIt.IsAtEnd())
        {
        if(maskIt.Get() != 0)
          {
          m_StatsMin = imageIt.Get();
          m_StatsMax = imageIt.Get();
          break;
          }
        ++imageIt;
        ++maskIt;
        }
      while(!maskIt.IsAtEnd())
        {
        if(maskIt.Get() != 0)
          {
          if(imageIt.Get() < m_StatsMin)
            {
            m_StatsMin = imageIt.Get();
            }
          else if(imageIt.Get() > m_StatsMax)
            {
            m_StatsMax = imageIt.Get();
            }
          }
        ++imageIt;
        ++maskIt;
        }
      }
    else
      {
      imageIt.GoToBegin();
      m_StatsMin = imageIt.Get();
      m_StatsMax = imageIt.Get();
      ++imageIt;
      while(!imageIt.IsAtEnd())
        {
        if(imageIt.Get() < m_StatsMin)
          {
          m_StatsMin = imageIt.Get();
          }
        else if(imageIt.Get() > m_StatsMax)
          {
          m_StatsMax = imageIt.Get();
          }
        ++imageIt;
        }
      }
    }
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::GetMin(void) const
{
  return m_StatsMin;
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::GetMax(void) const
{
  return m_StatsMax;
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Evaluate(const PointType& point, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    //itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return EvaluateAtContinuousIndex(cIndex, scale);
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Evaluate(const PointType& point, const VectorType & v1, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    //itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return EvaluateAtContinuousIndex(cIndex, v1, scale);
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Evaluate(const PointType& point,
           const VectorType & v1, const VectorType & v2, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    //itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return EvaluateAtContinuousIndex(cIndex, v1, v2, scale);
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::EvaluateAtIndex(const IndexType& index, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return EvaluateAtContinuousIndex(cIndex, scale);
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::EvaluateAtIndex(const IndexType& index, const VectorType & v1, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return EvaluateAtContinuousIndex(cIndex, v1, scale);
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::EvaluateAtIndex(const IndexType& index,
                  const VectorType & v1, const VectorType & v2,
                  double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return EvaluateAtContinuousIndex(cIndex, v1, v2, scale);
}

/**
 * Evaluate the fonction at the specified point
 */
template <class TInputImage>
double
NJetImageFunction<TInputImage>
::EvaluateAtContinuousIndex(const ContinuousIndexType & cIndex,
                            double scale) const
{
  // EVALUATE
  double physGaussFactor = -0.5/(scale*scale);
  double physKernelRadiusSquared = scale*m_Extent * scale*m_Extent;

  double physDist;
  double pixelValue;
  double expValue;

  double v = 0;
  double vTotal = 0;

  Index<ImageDimension> xMin;
  Index<ImageDimension> xMax;
  Index<ImageDimension> xShift;

  for(unsigned int i=0; i<ImageDimension; i++)
    {
    xMin[i] = (int)vnl_math_floor(cIndex[i]
                          - (scale * m_Extent / m_InputImageSpacing[i]));
    if(xMin[i]<m_InputImageMinX[i])
      {
      xMin[i]=m_InputImageMinX[i];
      }
    xShift[i] = xMin[i];

    xMax[i] = (int)vnl_math_ceil(cIndex[i]
                         + (scale * m_Extent / m_InputImageSpacing[i]));
    if(xMax[i] > (int) m_InputImageMaxX[i])
      {
      xMax[i]= m_InputImageMaxX[i];
      }
    }

  bool done = false;
  while(!done)
    {
    if(!m_UseInputImageMask
       || (m_UseInputImageMask && m_InputImageMask->GetPixel(xShift)>0))
      {
      physDist = 0;
      for(unsigned int i=0; i< ImageDimension; i++)
        {
        physDist += (cIndex[i]-xShift[i]) * (cIndex[i]-xShift[i])
                                          * m_InputImageSpacingSquared[i];
        }

      if(physDist <= physKernelRadiusSquared)
        {
        pixelValue = m_InputImage->GetPixel( xShift );
        expValue = vcl_exp(physGaussFactor*physDist);

        vTotal += vnl_math_abs(expValue);
        v += pixelValue * expValue;
        }
      }

    unsigned int dimI = 0;
    xShift[dimI]++;
    while( !done && xShift[dimI]>xMax[dimI] )
      {
      xShift[dimI] = xMin[dimI];
      dimI++;
      if( dimI < ImageDimension )
        {
        xShift[dimI]++;
        }
      else
        {
        done = true;
        }
      }
    }

  if(vTotal == 0)
    {
    //itkWarningMacro(<< "wTotal = 0 : only zero-value pixels encountered");
    m_MostRecentIntensity = 0.0;
    }
  else
    {
    m_MostRecentIntensity = v/vTotal;
    }

  return m_MostRecentIntensity;
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::EvaluateAtContinuousIndex(const ContinuousIndexType & cIndex,
                            const VectorType & v1,
                            double scale) const
{
  // EVALUATE
  if(m_UseProjection)
    {
    return EvaluateAtContinuousIndex(cIndex, scale);
    }
  else
    {
    scale = scale / 2;

    double val0 = this->EvaluateAtContinuousIndex(cIndex, scale);

    double step = 2.0673*scale;
    ContinuousIndexType tempI;
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v1[i];
      }
    double val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v1[i];
      }
    double val2 = this->EvaluateAtContinuousIndex(tempI, scale);

    m_MostRecentIntensity =  (val0 + 0.5455*val1 + 0.5455*val2)
                             / (1+2*0.5455) / 2;

    return m_MostRecentIntensity;
    }
}


template <class TInputImage>
double
NJetImageFunction<TInputImage>
::EvaluateAtContinuousIndex(const ContinuousIndexType & cIndex,
                            const VectorType & v1, const VectorType & v2,
                            double scale) const
{
  // EVALUATE
  if(m_UseProjection)
    {
    return EvaluateAtContinuousIndex(cIndex, scale);
    }
  else
    {
    scale = scale / 2;

    double val0 = this->EvaluateAtContinuousIndex(cIndex, scale);

    double step = 2.0673*scale;
    ContinuousIndexType tempI;
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v1[i];
      }
    double val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v1[i];
      }
    double val2 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v2[i];
      }
    double val3 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v2[i];
      }
    double val4 = this->EvaluateAtContinuousIndex(tempI, scale);

    m_MostRecentIntensity = (2*val0 + 0.5455*(val1+val2+val3+val4))/(2+4*0.5455) / 2;

    return m_MostRecentIntensity;
    }
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Derivative(const PointType& point, double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return DerivativeAtContinuousIndex(cIndex, scale, d);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Derivative(const PointType& point, const VectorType & v1,
  double scale, typename NJetImageFunction<TInputImage>::VectorType & d
  ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    d.Fill( 0 );
    return 0.0;
    }

  return DerivativeAtContinuousIndex(cIndex, v1, scale, d);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Derivative(const PointType& point, const VectorType & v1,
  const VectorType & v2, double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    d.Fill( 0 );
    return 0.0;
    }

  return DerivativeAtContinuousIndex(cIndex, v1, v2, scale, d);
}


template <class TInputImage>
double
NJetImageFunction<TInputImage>
::DerivativeAtIndex(const IndexType& index, double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return DerivativeAtContinuousIndex(cIndex, scale, d);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::DerivativeAtIndex(const IndexType& index, const VectorType & v1,
  double scale, typename NJetImageFunction<TInputImage>::VectorType & d
  ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return DerivativeAtContinuousIndex(cIndex, v1, scale, d);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::DerivativeAtIndex(const IndexType& index,
  const VectorType & v1, const VectorType & v2, double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return DerivativeAtContinuousIndex(cIndex, v1, v2, scale, d);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::DerivativeAtContinuousIndex(const ContinuousIndexType & cIndex,
  double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  double val = 0;

  // VALUE AND DERIVATIVE
  double physGaussFactor = -0.5/(scale*scale);
  double physKernelRadiusSquared = scale*m_Extent * scale*m_Extent;

  double physDist = 0;
  double pixelValue = 0;
  double expValue = 0;
  double expValueD = 0;

  double v = 0;
  double vTotal = 0;

  itk::Vector<double, TInputImage::ImageDimension> dTotal;
  d.Fill(0);
  dTotal.Fill(0);

  Index<ImageDimension> xMin;
  Index<ImageDimension> xMax;
  Index<ImageDimension> xShift;
  int xRadius;

  for(unsigned int i=0; i<ImageDimension; i++)
    {
    xMin[i] = (int) vnl_math_floor(cIndex[i] - (scale * m_Extent
                                             / m_InputImageSpacing[i]));
    if(xMin[i]<m_InputImageMinX[i])
      {
      xMin[i]=m_InputImageMinX[i];
      }
    xShift[i] = xMin[i];
    xRadius = (int) vnl_math_floor(cIndex[i] - xMin[i]);

    xMax[i] = (int) vnl_math_ceil(cIndex[i] + xRadius);
    if(xMax[i] > (int) m_InputImageMaxX[i])
      {
      xMax[i]= m_InputImageMaxX[i];
      xRadius = (int) vnl_math_floor(xMax[i] - cIndex[i]);
      xMin[i] = (int) vnl_math_floor(cIndex[i] - xRadius);
      xShift[i] = xMin[i];
      }
    }

  bool done = false;
  while(!done)
    {
    if(!m_UseInputImageMask
       || (m_UseInputImageMask && m_InputImageMask->GetPixel(xShift)>0))
      {
      physDist = 0;
      for(unsigned int i=0; i< ImageDimension; i++)
        {
        physDist += (cIndex[i]-xShift[i]) * (cIndex[i]-xShift[i])
                                          * m_InputImageSpacingSquared[i];
        }

      if(physDist <= physKernelRadiusSquared)
        {
        pixelValue = m_InputImage->GetPixel( xShift );
        expValue = vcl_exp(physGaussFactor*physDist);

        v += pixelValue*expValue;
        vTotal += vnl_math_abs(expValue);

        for(unsigned int i=0; i< ImageDimension; i++)
          {
          expValueD = 2 * (cIndex[i]-xShift[i]) * m_InputImageSpacing[i]
                        * physGaussFactor
                        * expValue;
          d[i] += pixelValue * expValueD;
          dTotal[i] += vnl_math_abs(expValueD);
          }
        }
      }

    xShift[0]++;
    unsigned int i = 0;
    while( !done && xShift[i]>xMax[i] )
      {
      xShift[i] = xMin[i];
      i++;
      if( i < ImageDimension )
        {
        xShift[i]++;
        }
      else
        {
        done = true;
        }
      }
    }

  double dMag = 0;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    if(dTotal[i] != 0)
      {
      d[i] = d[i] / dTotal[i];
      }
    m_MostRecentDerivative[i] = d[i];
    dMag += d[i]*d[i];
    }
  if( dMag != 0 )
    {
    dMag = vcl_sqrt( dMag );
    }

  if(vTotal != 0)
    {
    val = v/vTotal;
    }
  m_MostRecentIntensity = val;

  return dMag;
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::DerivativeAtContinuousIndex(const ContinuousIndexType & cIndex,
  const VectorType & v1, double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  double val = 0;
  // VALUE AND DERIVATIVE
  if(m_UseProjection)
    {
    val = DerivativeAtContinuousIndex(cIndex, val, scale, d);
    double dp0 = 0;
    for(int i=0; i<ImageDimension; i++)
      {
      dp0 += v1[i] * d[i];
      }
    d.Fill(0);
    d[0] = dp0;
    }
  else
    {
    double val0;
    double val1;
    double val2;

    scale = scale / 2;

    val0 = this->EvaluateAtContinuousIndex(cIndex, scale);

    double step = 2.0673 * scale;
    ContinuousIndexType tempI;
    for(int i=0; i<ImageDimension; i++)
      {
      d[i] = 0;
      tempI[i] = cIndex[i] - step*v1[i];
      }
    val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v1[i];
      }
    val2 = this->EvaluateAtContinuousIndex(tempI, scale);

    d.Fill( 0 );
    d[0] = (val2-val1)/ 2 / 2;

    val = (val0 + 0.5455*val1 + 0.5455*val2)/(1+2*0.5455) / 2;
    }

  m_MostRecentIntensity = val;

  m_MostRecentDerivative.Fill(0);
  m_MostRecentDerivative[0] = d[0];

  return vnl_math_abs( d[0] );
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::DerivativeAtContinuousIndex(const ContinuousIndexType & cIndex,
  const VectorType & v1, const VectorType & v2, double scale,
  typename NJetImageFunction<TInputImage>::VectorType & d ) const
{
  double val = 0;
  // VALUE AND DERIVATIVE

  if(m_UseProjection)
    {
    val = DerivativeAtContinuousIndex(cIndex, scale, d);
    double dp0 = 0;
    double dp1 = 0;
    for(int i=0; i<ImageDimension; i++)
      {
      dp0 += v1[i] * d[i];
      dp1 += v2[i] * d[i];
      }
    d.Fill(0);
    d[0] = dp0;
    d[1] = dp1;
    }
  else
    {
    double val0;
    double val1;
    double val2;
    double val3;
    double val4;

    scale = scale / 2;

    val0 = this->EvaluateAtContinuousIndex(cIndex, scale);

    double step = 2.0673 * scale;
    ContinuousIndexType tempI;
    for(int i=0; i<ImageDimension; i++)
      {
      d[i] = 0;
      tempI[i] = cIndex[i] - step*v1[i];
      }
    val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v1[i];
      }
    val2 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v2[i];
      }
    val3 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v2[i];
      }
    val4 = this->EvaluateAtContinuousIndex(tempI, scale);

    d.Fill(0);
    d[0] = (val2-val1)/ 2 / 2;
    d[1] = (val4-val3)/ 2 / 2;

    val = (2*val0 + 0.5455*(val1+val2+val3+val4))/(2+0.5455*4) / 2;
    }

  m_MostRecentIntensity = val;

  m_MostRecentDerivative.Fill(0);
  m_MostRecentDerivative[0] = d[0];
  m_MostRecentDerivative[1] = d[1];

  double dMag = (d[0]*d[0] + d[1]*d[1]);
  if( dMag != 0 )
    {
    dMag = vcl_sqrt( dMag );
    }

  return dMag;
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Hessian(const PointType& point, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & h ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    h.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    h.Fill( 0 );
    return 0.0;
    }

  return HessianAtContinuousIndex(cIndex, scale, h);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Hessian(const PointType& point, const VectorType & v1, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & h ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    h.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    h.Fill( 0 );
    return 0.0;
    }

  return HessianAtContinuousIndex(cIndex, v1, scale, h);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Hessian(const PointType& point,
  const VectorType & v1, const VectorType & v2, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & h ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    h.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    h.Fill( 0 );
    return 0.0;
    }

  return HessianAtContinuousIndex(cIndex, v1, v2, scale, h);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::HessianAtIndex(const IndexType& index, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & h ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    h.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return HessianAtContinuousIndex(cIndex, scale, h);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::HessianAtIndex(const IndexType& index, const VectorType & v1,
  double scale, typename NJetImageFunction<TInputImage>::MatrixType & h
  ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    h.Fill( 0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return HessianAtContinuousIndex(cIndex, v1, scale, h);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::HessianAtIndex(const IndexType& index,
  const VectorType & v1, const VectorType & v2, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & h ) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    h.Fill( 0.0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return HessianAtContinuousIndex(cIndex, v1, v2, scale, h);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::HessianAtContinuousIndex(const ContinuousIndexType & cIndex,
  double scale, typename NJetImageFunction<TInputImage>::MatrixType & h
  ) const
{
  VectorType d;
  JetAtContinuousIndex(cIndex, d, h, scale);
  double mag = 0;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    mag += h[i][i] * h[i][i];
    }
  // Should return the determinant
  return vcl_sqrt( mag );
}


template <class TInputImage>
double
NJetImageFunction<TInputImage>
::HessianAtContinuousIndex(const ContinuousIndexType & cIndex,
  const VectorType & v1, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & m ) const
{
  // HESSIAN
  if(m_UseProjection)
    {
    HessianAtContinuousIndex(cIndex, scale, m);
    vnl_symmetric_eigensystem< double > eigSys(m.GetVnlMatrix());

    assert( eigSys.get_eigenvalue(0) <= eigSys.get_eigenvalue(1) );

    double dp = 0;
    for(int i=0; i<ImageDimension; i++)
      {
      dp += dot_product(v1.GetVnlVector(), eigSys.get_eigenvector(i))
              * eigSys.get_eigenvalue(i);
      }

    m.Fill( 0 );
    m[0][0] = dp;
    }
  else
    {
    scale = scale / 2;
    double step = 3.7 * scale;
    ContinuousIndexType tempI;
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v1[i];
      }
    double val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v1[i];
      }
    double val2 = this->EvaluateAtContinuousIndex(tempI, scale);
    double val3 = this->EvaluateAtContinuousIndex(cIndex, scale);

    m.Fill(0);
    m[0][0] = (0.5*val1+0.5*val2-val3)/2;
    }
  m_MostRecentHessian.Fill( 0 );
  m_MostRecentHessian[0][0] = m[0][0];

  return vnl_math_abs( m[0][0] );
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::HessianAtContinuousIndex(const ContinuousIndexType & cIndex,
  const VectorType & v1, const VectorType & v2, double scale,
  typename NJetImageFunction<TInputImage>::MatrixType & m ) const
{
  // HESSIAN
  if(m_UseProjection)
    {
    HessianAtContinuousIndex(cIndex, scale, m);
    vnl_symmetric_eigensystem< double > eigSys(m.GetVnlMatrix());

    assert( eigSys.get_eigenvalue(0) <= eigSys.get_eigenvalue(1) );

    double dp0 = 0;
    double dp1 = 0;
    for(int i=0; i<ImageDimension; i++)
      {
      dp0 += dot_product(v1.GetVnlVector(), eigSys.get_eigenvector(i))
              * eigSys.get_eigenvalue(i);
      dp1 += dot_product(v2.GetVnlVector(), eigSys.get_eigenvector(i))
              * eigSys.get_eigenvalue(i);
      }

    m[0][0] = dp0;
    m[1][1] = dp1;
    }
  else
    {
    scale = scale / 2;
    double step = 3.7 * scale;
    ContinuousIndexType tempI;
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v1[i];
      }
    double val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v1[i];
      }
    double val2 = this->EvaluateAtContinuousIndex(tempI, scale);
    double val3 = this->EvaluateAtContinuousIndex(cIndex, scale);

    m.Fill(0);
    m[0][0] = (0.5*val1+0.5*val2-val3)/2;

    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] - step*v2[i];
      }
    val1 = this->EvaluateAtContinuousIndex(tempI, scale);
    for(int i=0; i<ImageDimension; i++)
      {
      tempI[i] = cIndex[i] + step*v2[i];
      }
    val2 = this->EvaluateAtContinuousIndex(tempI, scale);

    m[1][1] = (0.5*val1+0.5*val2-val3)/2;
    }
  m_MostRecentHessian.Fill( 0 );
  m_MostRecentHessian[0][0] = m[0][0];
  m_MostRecentHessian[1][1] = m[1][1];
  return vcl_sqrt( m[0][0]*m[0][0] + m[1][1]*m[1][1] );
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Jet(const PointType& point, VectorType & d, MatrixType & h,
      double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0.0 );
    h.Fill( 0.0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    d.Fill( 0.0 );
    h.Fill( 0.0 );
    return 0.0;
    }

  return JetAtContinuousIndex(cIndex, d, h, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::JetAtIndex(const IndexType& index, VectorType & d, MatrixType & h,
             double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    d.Fill( 0.0 );
    h.Fill( 0.0 );
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return JetAtContinuousIndex(cIndex, d, h, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::JetAtContinuousIndex(const ContinuousIndexType & cIndex,
                                        VectorType & d,
                                        MatrixType & h,
                                        double scale) const
{
  // JET
  double physGaussFactor = -0.5/(scale*scale);
  double physKernelRadiusSquared = scale*m_Extent * scale*m_Extent;

  double physDist = 0;
  double pixelValue = 0;
  double expValue = 0;
  double expValueD = 0;

  double v = 0;
  double vTotal = 0;

  itk::Vector<double, TInputImage::ImageDimension> dTotal;
  d.Fill(0);
  dTotal.Fill(0);

  MatrixType hTotal;
  h.Fill(0);
  hTotal.Fill(0);

  Index<ImageDimension> xMin;
  Index<ImageDimension> xMax;
  Index<ImageDimension> xShift;
  int xRadius;

  for(unsigned int i=0; i<ImageDimension; i++)
    {
    xMin[i] = (int) vnl_math_floor(cIndex[i]
                          - (scale * m_Extent / m_InputImageSpacing[i]));
    if(xMin[i]<m_InputImageMinX[i])
      {
      xMin[i] = m_InputImageMinX[i];
      }
    xShift[i] = xMin[i];
    xRadius = (int) vnl_math_floor( cIndex[i] - xMin[i] );

    xMax[i] = (int) vnl_math_ceil( cIndex[i] + xRadius );
    if(xMax[i] > (int) m_InputImageMaxX[i])
      {
      xMax[i] = m_InputImageMaxX[i];
      xRadius = (int) vnl_math_floor(xMax[i] - cIndex[i]);
      xMin[i] = (int) vnl_math_floor(cIndex[i] - xRadius);
      xShift[i] = xMin[i];
      }
    }

  bool done = false;
  while(!done)
    {
    if(!m_UseInputImageMask
       || (m_UseInputImageMask && m_InputImageMask->GetPixel(xShift)>0))
      {
      physDist = 0;
      for(unsigned int i=0; i< ImageDimension; i++)
        {
        physDist += (cIndex[i]-xShift[i]) * (cIndex[i]-xShift[i])
                                          * m_InputImageSpacingSquared[i];
        }

      if(physDist <= physKernelRadiusSquared)
        {
        pixelValue = m_InputImage->GetPixel( xShift );
        expValue = vcl_exp(physGaussFactor*physDist);

        v += pixelValue*expValue;
        vTotal += vnl_math_abs(expValue);

        for(unsigned int i=0; i< ImageDimension; i++)
          {
          expValueD = 2 * (cIndex[i]-xShift[i]) * m_InputImageSpacing[i]
                        * physGaussFactor
                        * expValue;
          d[i] += pixelValue * expValueD;
          dTotal[i] += vnl_math_abs(expValueD);

          expValueD = (2 * (cIndex[i]-xShift[i]) * (xShift[i] - cIndex[i])
                         * m_InputImageSpacingSquared[i]
                         * physGaussFactor
                         + 1) * 2
                      * physGaussFactor
                      * expValue;
          h[i][i] += pixelValue * expValueD;
          hTotal[i][i] += vnl_math_abs(expValueD);
          for(unsigned int j=i+1; j< ImageDimension; j++)
            {
            expValueD = 4 * (cIndex[i]-xShift[i]) * m_InputImageSpacing[i]
                          * (cIndex[j]-xShift[j]) * m_InputImageSpacing[j]
                          * physGaussFactor * physGaussFactor
                          * expValue;
            h[i][j] += pixelValue * expValueD;
            hTotal[i][j] += vnl_math_abs(expValueD);
            }
          }
        }
      }

    xShift[0]++;
    unsigned int i = 0;
    while( !done && xShift[i]>xMax[i] )
      {
      xShift[i] = xMin[i];
      i++;
      if( i < ImageDimension )
        {
        xShift[i]++;
        }
      else
        {
        done = true;
        }
      }
    }

  if(vTotal != 0)
    {
    v = v / vTotal;
    }
  else
    {
    v = 0;
    }
  m_MostRecentIntensity = v;

  for(unsigned int i=0; i<ImageDimension; i++)
    {
    if(dTotal[i] != 0)
      {
      d[i] = d[i] / dTotal[i];
      }
    m_MostRecentDerivative[i] = d[i];

    if(hTotal[i][i] != 0)
      {
      h[i][i] = h[i][i] / hTotal[i][i];
      }
    else
      {
      h[i][i] = 0;
      }
    m_MostRecentHessian[i][i] = h[i][i];
    for(unsigned int j=i+1; j<ImageDimension; j++)
      {
      if(hTotal[i][j] != 0)
        {
        h[i][j] = h[i][j] / hTotal[i][j];
        h[j][i] = h[i][j];
        }
      else
        {
        h[i][j] = 0;
        h[j][i] = 0;
        }
      m_MostRecentHessian[i][j] = h[i][j];
      m_MostRecentHessian[j][i] = h[j][i];
      }
    }

  return v;
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Ridgeness(const PointType& point, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return RidgenessAtContinuousIndex(cIndex, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Ridgeness(const PointType& point, const VectorType & v1, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return RidgenessAtContinuousIndex(cIndex, v1, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::Ridgeness(const PointType& point,
            const VectorType & v1, const VectorType & v2, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  if(!m_InputImage->TransformPhysicalPointToContinuousIndex(point, cIndex))
    {
    itkWarningMacro(<< "Cannot convert point to continuous index");
    return 0.0;
    }

  return RidgenessAtContinuousIndex(cIndex, v1, v2, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::RidgenessAtIndex(const IndexType& index, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return RidgenessAtContinuousIndex(cIndex, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::RidgenessAtIndex(const IndexType& index,
                   const VectorType & v1, double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return RidgenessAtContinuousIndex(cIndex, v1, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::RidgenessAtIndex(const IndexType& index,
                   const VectorType & v1, const VectorType & v2,
                   double scale) const
{
  if( !m_InputImage )
    {
    itkWarningMacro(<< "Input image not set");
    return 0.0;
    }

  ContinuousIndexType cIndex;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    cIndex[i] = index[i];
    }

  return RidgenessAtContinuousIndex(cIndex, v1, v2, scale);
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::RidgenessAtContinuousIndex(const ContinuousIndexType & cIndex,
                             double scale) const
{
  // RIDGENESS
  VectorType d;
  MatrixType h;
  VectorType p;

  JetAtContinuousIndex(cIndex, d, h, scale);

  vnl_symmetric_eigensystem< double > eigSys(h.GetVnlMatrix());

  // Ensure ordering of eigenvalues; According to VNL documentation,
  // eigenvalues are in increasing order (with smallest eigenvalues
  // first)
  assert( eigSys.get_eigenvalue(0) <= eigSys.get_eigenvalue(1) );

  double dNorm = d.GetNorm();
  if( dNorm == 0 )
    {
    d.SetVnlVector( eigSys.get_eigenvector( ImageDimension-1 ) );
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      p[i] = 0;
      }
    }
  else
    {
    d.Normalize();
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      p[i] = 0;
      for( unsigned int j=0; j<ImageDimension; j++ )
        {
        p[i] += eigSys.get_eigenvector(i)[j] * d[j];
        }
      }
    }

  double sums = 0;
  double sumv = 0;
  int ridge = 1;
  for( unsigned int i=0; i<ImageDimension-1; i++ )
    {
    sums += p[i]*p[i];
    sumv += eigSys.get_eigenvalue(i)
            * eigSys.get_eigenvalue(i);
    if( eigSys.get_eigenvalue(i) >= 0 )
      {
      ridge = -1;
      }
    }
  sums /= (ImageDimension-1);

  double ridgeness = (1.0 - sums) * ridge;

  double roundness = 0;
  if( sumv != 0 )
    {
    if( ImageDimension > 2 )
      {
      if( eigSys.get_eigenvalue(0) != 0
          && eigSys.get_eigenvalue(ImageDimension-2) != 0 )
        {
        double r0 = ( eigSys.get_eigenvalue(0)
                      * eigSys.get_eigenvalue(0) )
                    / ( eigSys.get_eigenvalue(ImageDimension-2)
                        * eigSys.get_eigenvalue(ImageDimension-2) );
        double r1 = ( eigSys.get_eigenvalue(ImageDimension-2)
                        * eigSys.get_eigenvalue(ImageDimension-2) )
                    / ( eigSys.get_eigenvalue(0)
                        * eigSys.get_eigenvalue(0) );
        roundness = ( r0 < r1 ) ? r0 : r1;
        }
      }
    else
      {
      if( eigSys.get_eigenvalue(0) != 0
          && eigSys.get_eigenvalue(ImageDimension-1) != 0 )
        {
        double r0 = ( eigSys.get_eigenvalue(0)
                      * eigSys.get_eigenvalue(0) )
                    / ( eigSys.get_eigenvalue(ImageDimension-1)
                        * eigSys.get_eigenvalue(ImageDimension-1) );
        double r1 = ( eigSys.get_eigenvalue(ImageDimension-1)
                        * eigSys.get_eigenvalue(ImageDimension-1) )
                    / ( eigSys.get_eigenvalue(0)
                        * eigSys.get_eigenvalue(0) );
        roundness = ( r0 < r1 ) ? r0 : r1;
        roundness = ( 1 - roundness );
        }
      }
    }


  double levelness = 0;
  if( sumv != 0 )
    {
    levelness = sumv / ( sumv
                         + eigSys.get_eigenvalue(ImageDimension-1)
                          * eigSys.get_eigenvalue(ImageDimension-1) );
    }

  double curvature = vcl_sqrt( sumv / ( ImageDimension-1 ) );


  m_MostRecentRidgeness = ridgeness;
  m_MostRecentRidgeRoundness = roundness;
  m_MostRecentRidgeLevelness = levelness;
  m_MostRecentRidgeCurvature = curvature;
  m_MostRecentRidgeTangent.SetVnlVector( eigSys.get_eigenvector( ImageDimension-1 ) );

  return m_MostRecentRidgeness;
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::RidgenessAtContinuousIndex(const ContinuousIndexType & cIndex,
                             const VectorType & v1, double scale) const
{
  // RIDGENESS
  itk::Vector<double, TInputImage::ImageDimension> d;

  double val;
  MatrixType h;
  VectorType p;
  VectorType vv;

  val = JetAtContinuousIndex(cIndex, d, h, scale);

  vnl_symmetric_eigensystem< double > eigSys(h.GetVnlMatrix());

  assert( eigSys.get_eigenvalue(0) <= eigSys.get_eigenvalue(1) );

  if( d.GetNorm() != 0 )
    {
    d.Normalize();
    }
  else
    {
    d.SetVnlVector( eigSys.get_eigenvector(ImageDimension-1) );
    }

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    p[i] = 0;
    double dp = 0;
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      dp += eigSys.get_eigenvector(i)[j] * v1[j];
      }
    dp = vnl_math_abs( dp );
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      p[i] += dp * eigSys.get_eigenvector(i)[j] * d[j];
      }
    vv[i] = dp * eigSys.get_eigenvalue(i);
    }

  double sums = 0;
  double sumv = 0;
  int ridge = 1;
  for( unsigned int i=0; i<ImageDimension-1; i++ )
    {
    sums += p[i]*p[i];
    sumv += vv[i] * vv[i];
    if( vv[i] >= 0 )
      {
      ridge = -1;
      }
    }
  sums /= (ImageDimension-1);
  if( sumv != 0 )
    {
    sumv /= ( sumv + vv[ImageDimension-1] * vv[ImageDimension-1] );
    }

  val = (1.0 - sums) * sumv * ridge;

  m_MostRecentRidgeness = val;

  return val;
}

template <class TInputImage>
double
NJetImageFunction<TInputImage>
::RidgenessAtContinuousIndex(const ContinuousIndexType & cIndex,
                             const VectorType & v1, const VectorType & v2,
                             double scale) const
{
  // RIDGENESS
  itk::Vector<double, TInputImage::ImageDimension> d;

  double val;
  MatrixType h;
  VectorType p;
  VectorType vv;

  val = JetAtContinuousIndex(cIndex, d, h, scale);

  vnl_symmetric_eigensystem< double > eigSys(h.GetVnlMatrix());

  assert( eigSys.get_eigenvalue(0) <= eigSys.get_eigenvalue(1) );

  if( d.GetNorm() != 0 )
    {
    d.Normalize();
    }
  else
    {
    d.SetVnlVector( eigSys.get_eigenvector(ImageDimension-1) );
    }

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    p[i] = 0;
    double dp = 0;
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      dp += eigSys.get_eigenvector(i)[j] * v1[j];
      }
    dp = vnl_math_abs( dp );
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      p[i] += dp * eigSys.get_eigenvector(i)[j] * d[j];
      }
    vv[i] = dp * eigSys.get_eigenvalue(i);
    dp = 0;
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      dp += eigSys.get_eigenvector(i)[j] * v2[j];
      }
    dp = vnl_math_abs( dp );
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      p[i] += dp * eigSys.get_eigenvector(i)[j] * d[j];
      }
    vv[i] += dp * eigSys.get_eigenvalue(i);
    }

  double sums = 0;
  double sumv = 0;
  int ridge = 1;
  for( unsigned int i=0; i<ImageDimension-1; i++ )
    {
    sums += p[i]*p[i];
    sumv += vv[i] * vv[i];
    if( vv[i] >= 0 )
      {
      ridge = -1;
      }
    }
  sums /= (ImageDimension-1);
  if( sumv != 0 )
    {
    sumv /= (sumv + vv[ImageDimension-1] * vv[ImageDimension-1] );
    }

  val = (1.0 - sums) * sumv * ridge;

  m_MostRecentRidgeness = val;

  return val;
}

template <class TInputImage>
typename NJetImageFunction<TInputImage>::InputImagePointer
NJetImageFunction<TInputImage>
::ScaleSubsample(double factor)
{
  if(m_InputImage)
    {
    typedef typename InputImageType::PixelType PixelType;
    typename InputImageType::SizeType size;
    size = m_InputImage->GetLargestPossibleRegion().GetSize();
    typename InputImageType::SpacingType spacing;
    spacing = m_InputImage->GetSpacing();
    for(int i=0; i<ImageDimension; i++)
      {
      size[i] = (long unsigned int)(size[i] / factor);
      spacing[i] *= factor;
      }

    this->SetExtent(2);

    typename InputImageType::Pointer newImage = InputImageType::New();
    newImage->SetRegions( size );
    newImage->SetSpacing( spacing );
    newImage->SetOrigin( m_InputImage->GetOrigin() );
    newImage->Allocate();

    typedef ImageRegionIteratorWithIndex<InputImageType> IteratorType;
    IteratorType
        imageIt(newImage,
                newImage->GetLargestPossibleRegion());

    if(m_UseInputImageMask)
      {
      imageIt.GoToBegin();
      while(!imageIt.IsAtEnd())
        {
        typename InputImageType::PointType p;
        IndexType i;
        ContinuousIndexType ic;
        i = imageIt.GetIndex();
        newImage->TransformIndexToPhysicalPoint(i, p);
        m_InputImage->TransformPhysicalPointToIndex(p, i);
        m_InputImage->TransformPhysicalPointToContinuousIndex(p, ic);
        if(m_InputImageMask->GetPixel(i) != 0)
          {
          imageIt.Set( (PixelType)
                        this->EvaluateAtContinuousIndex(ic, spacing[0]/2) );
          }
        ++imageIt;
        }
      }
    else
      {
      imageIt.GoToBegin();
      while(!imageIt.IsAtEnd())
        {
        typename InputImageType::PointType p;
        IndexType i;
        ContinuousIndexType ic;
        i = imageIt.GetIndex();
        newImage->TransformIndexToPhysicalPoint(i, p);
        m_InputImage->TransformPhysicalPointToContinuousIndex(p, ic);
        imageIt.Set( (PixelType)
                      this->EvaluateAtContinuousIndex(ic, spacing[0]/2) );
        ++imageIt;
        }
      }

    return newImage;
    }

  return 0;
}

} // namespace tube

} // namespace itk

#endif
