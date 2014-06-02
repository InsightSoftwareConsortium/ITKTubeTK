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
#include "itktubePadImageFilter.h"
#include "itktubeRegionFromReferenceImageFilter.h"

#include "itkParametricImageSource.h"
#include "itkSize.h"

namespace itk {

namespace tube {
//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::FFTGaussianDerivativeIFFTFilter()
{
  m_FFTImage = NULL;
  m_IFFTImage = NULL;
  m_KernelImage = NULL;

  this->m_Orders.Fill(0);
  this->m_Sigmas.Fill(0);

  this->m_LastInputImage = NULL;
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::ApplyFFT()
{
  typedef PadImageFilter< TInputImage, TInputImage >   PadFilterType;
  typename PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput1( this->GetInput() );
  padFilter->SetInput2( this->GetInput() );
  padFilter->SetGreatestPrimeFactor( 5 );
  padFilter->SetPadMethod( PadFilterType::ZERO_FLUX_NEUMANN );
  padFilter->Update();

  typename FFTInputType::Pointer fftFilter = FFTInputType::New();
  fftFilter->SetInput( padFilter->GetOutput() );
  fftFilter->Update();

  m_FFTImage = fftFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::CreateGaussianDerivative()
{
  typename GaussianDerivativeImageSourceType::Pointer gaussSource =
    GaussianDerivativeImageSourceType::New();
  gaussSource->SetNormalized( true );

  const typename ComplexImageType::RegionType inputRegion(
    this->GetInput()->GetLargestPossibleRegion() );
  const typename ComplexImageType::SizeType inputSize
    = inputRegion.GetSize();

  const typename ComplexImageType::RegionType fftRegion(
    m_FFTImage->GetLargestPossibleRegion() );
  const typename ComplexImageType::SizeType fftSize
    = fftRegion.GetSize();
  const typename ComplexImageType::SpacingType fftSpacing =
    m_FFTImage->GetSpacing();
  const typename ComplexImageType::PointType fftOrigin =
    m_FFTImage->GetOrigin();
  const typename ComplexImageType::DirectionType fftDirection =
    m_FFTImage->GetDirection();

  gaussSource->SetIndex( fftRegion.GetIndex() );
  gaussSource->SetSize( fftSize );
  gaussSource->SetSpacing( fftSpacing );
  gaussSource->SetOrigin( fftOrigin );
  gaussSource->SetDirection( fftDirection );

  typename GaussianDerivativeImageSourceType::PointType mean;
  typename GaussianDerivativeImageSourceType::IndexType meanIndex;
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    const int halfLength = (inputSize[ii]  / 2.0);
    meanIndex[ii] = inputRegion.GetIndex()[ii] + halfLength;
    }
  this->GetInput()->TransformIndexToPhysicalPoint( meanIndex, mean );

  gaussSource->SetSigmas( m_Sigmas );
  gaussSource->SetMean( mean );
  gaussSource->SetOrders( this->m_Orders );

  gaussSource->Update();
  m_KernelImage = gaussSource->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::ApplyGaussianDerivativeIFFT()
{
  CreateGaussianDerivative();

  typename FFTShiftFilterType::Pointer fftShiftFilter =
    FFTShiftFilterType::New();
  fftShiftFilter->SetInput( m_KernelImage );
  fftShiftFilter->Update();

  typename FFTRealType::Pointer fftFilter = FFTRealType::New();
  fftFilter->SetInput( fftShiftFilter->GetOutput() );
  fftFilter->Update();
  typename ComplexImageType::Pointer fftKernel = fftFilter->GetOutput();

  typename MultiplyFilterType::Pointer multiplyFilter =
    MultiplyFilterType::New();
  multiplyFilter->SetInput1( m_FFTImage );
  multiplyFilter->SetInput2( fftKernel );
  multiplyFilter->Update();

  typename InverseFFTFilterType::Pointer 
    iFFTFilter = InverseFFTFilterType::New();
  iFFTFilter->SetInput( multiplyFilter->GetOutput() );
  iFFTFilter->Update();

  m_IFFTImage = iFFTFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::GenerateData()
{
  if( m_LastInputImage != this->GetInput() )
    {
    m_LastInputImage = this->GetInput();
    ApplyFFT();
    }

  ApplyGaussianDerivativeIFFT();

  typedef RegionFromReferenceImageFilter< TOutputImage > 
    RegionFromFilterType;
  typename RegionFromFilterType::Pointer regionFrom = 
    RegionFromFilterType::New();

  regionFrom->SetInput1( m_IFFTImage );
  regionFrom->SetInput2( this->GetInput() );
  regionFrom->Update();

  this->GraftOutput( regionFrom->GetOutput() );
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_FFTImage.IsNotNull() )
    {
    os << indent << "FFT Image           : " << m_FFTImage << std::endl;
    }
  else
    {
    os << indent << "FFT Image           : NULL" << std::endl;
    }
  if( m_IFFTImage.IsNotNull() )
    {
    os << indent << "Inverse FFT Image   : " << m_IFFTImage << std::endl;
    }
  else
    {
    os << indent << "Inverse FFT Image   : NULL" << std::endl;
    }
  if( m_KernelImage.IsNotNull() )
    {
    os << indent << "Kernel Image        : " << m_KernelImage << std::endl;
    }
  else
    {
    os << indent << "Kernel Image        : NULL" << std::endl;
    }
  os << indent << "Orders              : " << m_Orders << std::endl;
  os << indent << "Sigmas               : " << m_Sigmas << std::endl;
  os << indent << "Last Input Image    : " << m_LastInputImage << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif
