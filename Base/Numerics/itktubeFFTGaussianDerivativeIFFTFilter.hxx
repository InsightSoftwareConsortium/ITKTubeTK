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
#ifndef __itktubeFFTGaussianDerivativeIFFTFilter_hxx
#define __itktubeFFTGaussianDerivativeIFFTFilter_hxx

#include "itktubeFFTGaussianDerivativeIFFTFilter.h"
#include "itktubeGaussianDerivativeImageSource.h"
#include "itkParametricImageSource.h"
#include "itkFixedArray.h"
#include "itkSize.h"

namespace itk {

namespace tube {
//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::FFTGaussianDerivativeIFFTFilter()
{
  this->m_GaussianDerivative = GaussianDerivativeType::New();
  this->m_FFTFilter = FFTType::New();
  this->m_InverseFFTFilter = InverseFFTType::New();
  this->m_ImageData = this->GetInput();
  this->m_Orders = 0;
  this->m_Sigma = 0;
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::applyFFT()
{
  //Compute the FFT
  this->m_FFTFilter->SetInput(this->m_ImageData);
  this->m_FFTFilter->Update();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::applyIFFT()
{
  this->m_InverseFFTFilter->SetInput(this->m_GaussianDerivative->GetOutput());
  this->m_InverseFFTFilter->Update();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::applyGaussianDerivativeFilterIFFT()
{
  createGaussianDerivative();
  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
  fftShiftFilter->SetInput( this->m_GaussianDerivative->GetOutput() );

  MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
  multiplyFilter->SetInput1( this->m_FFTFilter->GetOutput() );
  multiplyFilter->SetInput2( fftShiftFilter->GetOutput() );

  this->m_InverseFFTFilter->SetInput(fftShiftFilter->GetOutput());
  this->m_InverseFFTFilter->Update();

}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::createGaussianDerivative()
{
  this->m_GaussianDerivative->SetNormalized( true );
  ComplexImageType::ConstPointer transformedInput
    = this->m_FFTFilter->GetOutput();
  const ComplexImageType::RegionType inputRegion(
    transformedInput->GetLargestPossibleRegion() );
  const ComplexImageType::SizeType inputSize
    = inputRegion.GetSize();
  const ComplexImageType::SpacingType inputSpacing =
    transformedInput->GetSpacing();
  const ComplexImageType::PointType inputOrigin =
    transformedInput->GetOrigin();
  const ComplexImageType::DirectionType inputDirection =
    transformedInput->GetDirection();

  this->m_GaussianDerivative->SetSize( inputSize );
  this->m_GaussianDerivative->SetSpacing( inputSpacing );
  this->m_GaussianDerivative->SetOrigin( inputOrigin );
  this->m_GaussianDerivative->SetDirection( inputDirection );
  GaussianDerivativeType::PointType mean;

  for( unsigned int ii = 0; ii < 3; ++ii )
    {
    const double halfLength = inputSize[ii]  / 2.0;
    mean[ii] = inputOrigin[ii] + halfLength;
    }
  mean = inputDirection * mean;
  this->m_GaussianDerivative->SetSigma( this->m_Sigma );
  this->m_GaussianDerivative->SetMean( mean );
  this->m_GaussianDerivative->SetOrdersVector( this->m_Orders);
  this->m_GaussianDerivative->Update();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::GenerateData()
{
  if(this->GetInput() != this->m_ImageData )
    {
    this->m_ImageData = this->GetInput();
    applyFFT();
    }
  applyGaussianDerivativeFilterIFFT();
  this->GraftOutput(this->m_InverseFFTFilter->GetOutput());
}

} // End namespace tube

} // End namespace itk

#endif
