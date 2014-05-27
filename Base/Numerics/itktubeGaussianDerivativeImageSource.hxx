/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itktubeGaussianDerivativeImageSource_hxx
#define __itktubeGaussianDerivativeImageSource_hxx

#include "itktubeGaussianDerivativeImageSource.h"
#include <itkGaussianSpatialFunction.h>
#include <itkImageRegionIterator.h>
#include <itkObjectFactory.h>
#include <itkProgressReporter.h>

namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template< typename TOutputImage >
GaussianDerivativeImageSource< TOutputImage >
::GaussianDerivativeImageSource()
{
  // Gaussian parameters, defined so that the gaussian
  // is centered in the default image
  m_Mean.Fill(32.0);
  m_Sigma.Fill(16.0);
  m_Scale = 255.0;

  m_Normalized = false;
  m_Order = 0;
  m_OrdersVector.Filled(0);
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Gaussian sigma: [";
  for ( unsigned int ii = 0; ii < NDimensions; ++ii )
    {
    os << m_Sigma[ii];
    if( ii != NDimensions - 1 )
      {
      os << ", ";
      }
    }
  os << "]" << std::endl;

  os << indent << "Gaussian mean: [";
  for ( unsigned int ii = 0; ii < NDimensions; ++ii )
    {
    os << m_Mean[ii];
    if( ii != NDimensions - 1 )
      {
      os << ", ";
      }
    }
  os << "]" << std::endl;

  os << indent << "Gaussian scale: " << m_Scale << std::endl;
  os << indent << "Normalized Gaussian?: " << m_Normalized << std::endl;
  os << indent << "Gaussian order: " << m_Order << std::endl;
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::SetParameters(const ParametersType & parameters)
{
  ArrayType sigma, mean;
  for ( unsigned int i = 0; i < ArrayType::Length; i++ )
    {
    sigma[i] = parameters[i];
    mean[i]  = parameters[i + ArrayType::Length];
    }
  this->SetSigma( sigma );
  this->SetMean( mean );

  double scale = parameters[2*ArrayType::Length];
  this->SetScale( scale );
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
typename GaussianDerivativeImageSource< TOutputImage >::ParametersType
GaussianDerivativeImageSource< TOutputImage >
::GetParameters() const
{
  ParametersType parameters( 2*ArrayType::Length + 1 );
  for ( unsigned int i = 0; i < ArrayType::Length; i++ )
    {
    parameters[i] = m_Sigma[i];
    parameters[i + ArrayType::Length] = m_Mean[i];
    }
  parameters[2*ArrayType::Length] = m_Scale;

  return parameters;
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
unsigned int
GaussianDerivativeImageSource< TOutputImage >
::GetNumberOfParameters() const
{
  return 2*ArrayType::Length + 1;
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::GenerateData()
{
  TOutputImage * outputPtr = this->GetOutput();
  std::cout <<"spacing-GenerateData"<<Superclass::GetSpacing();
  // allocate the output buffer
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();

  // Create an iterator that will walk the output region
  typedef itk::ImageRegionIterator< TOutputImage > OutputIterator;
  OutputIterator outIt = OutputIterator( outputPtr,
                                         outputPtr->GetRequestedRegion() );

  // The position at which the function is evaluated
  Point< double, TOutputImage::ImageDimension > evalPoint;

  ProgressReporter progress( this, 0,
                             outputPtr->GetRequestedRegion()
                             .GetNumberOfPixels() );
 //  Walk the output image, evaluating the spatial function at each pixel
  for (; !outIt.IsAtEnd(); ++outIt )
    {
    typename TOutputImage::IndexType index = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint(index, evalPoint);

    double prefixDenom = 1.0;

    if ( m_Normalized )
      {
      const double squareRootOfTwoPi = vcl_sqrt(2.0 * vnl_math::pi);

      for ( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
        {
        prefixDenom *= m_Sigma[i] * squareRootOfTwoPi;
        }
      }
    if( m_Order==1 || m_Order==2)
      {
      for ( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
        {
        prefixDenom *= vcl_pow(m_Sigma[i], 2*m_Order)/ (vcl_pow((- (evalPoint[i]
        - m_Mean[i])),m_Order) - (m_Order == 2 ? vcl_pow(m_Sigma[1], m_Order) : 0));
        }
      }

    for ( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
      {
      if(m_OrdersVector[i] != 0)
        {
        prefixDenom *= vcl_pow(m_Sigma[i], 2*m_OrdersVector[i])/ (vcl_pow((- (
        evalPoint[i]- m_Mean[i])),m_OrdersVector[i]) - (m_OrdersVector[i] == 2
        ? vcl_pow(m_Sigma[1], m_OrdersVector[i]) : 0));
        }
      }

    double suffixExp = 0;

    for ( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
      {
      suffixExp += ( evalPoint[i] - m_Mean[i] ) * ( evalPoint[i] - m_Mean[i] )
                   / ( 2 * m_Sigma[i] * m_Sigma[i] );
      }

    double value = m_Scale * ( 1 / prefixDenom ) * vcl_exp(-1 * suffixExp);

    // Set the pixel value to the function value
    outIt.Set( ( typename TOutputImage::PixelType )value );
    progress.CompletedPixel();
    }
}
} // End namespace tube

} // End namespace itk

#endif
