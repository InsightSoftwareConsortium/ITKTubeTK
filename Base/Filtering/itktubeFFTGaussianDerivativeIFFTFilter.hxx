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
  m_InputImageFFT = NULL;
  m_KernelImageFFT = NULL;
  m_ConvolvedImageFFT = NULL;
  m_ConvolvedImage = NULL;

  this->m_Orders.Fill(0);
  this->m_Sigmas.Fill(0);

  this->m_LastInputImage = NULL;
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::SetOrders( OrdersType & orders )
{
  m_Orders = orders;
  this->Modified();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::SetSigmas( SigmasType & sigmas )
{
  m_Sigmas = sigmas;
  this->Modified();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::ComputeInputImageFFT()
{
  std::cout << "FFTGaussianDerivativeIFFT: ComputeInputImageFFT" << std::endl;

  typedef PadImageFilter< InputImageType, RealImageType >
    PadFilterType;
  typename PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput( this->GetInput() );
  padFilter->SetGreatestPrimeFactor( 5 );
  padFilter->SetPadMethod( PadFilterType::ZERO_FLUX_NEUMANN );
  padFilter->Update();

  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  fftFilter->SetInput( padFilter->GetOutput() );
  fftFilter->Update();

  m_InputImageFFT = fftFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::ComputeKernelImageFFT()
{
  std::cout << "FFTGaussianDerivativeIFFT: ComputeKernelImageFFT" << std::endl;
  typename GaussianDerivativeImageSourceType::Pointer gaussSource =
    GaussianDerivativeImageSourceType::New();

  const typename ComplexImageType::RegionType inputRegion(
    this->GetInput()->GetLargestPossibleRegion() );
  const typename ComplexImageType::SizeType inputSize
    = inputRegion.GetSize();

  const typename ComplexImageType::RegionType fftRegion(
    m_InputImageFFT->GetLargestPossibleRegion() );
  const typename ComplexImageType::SizeType fftSize
    = fftRegion.GetSize();
  const typename ComplexImageType::SpacingType fftSpacing =
    m_InputImageFFT->GetSpacing();
  const typename ComplexImageType::PointType fftOrigin =
    m_InputImageFFT->GetOrigin();
  const typename ComplexImageType::DirectionType fftDirection =
    m_InputImageFFT->GetDirection();

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

  typename FFTShiftFilterType::Pointer fftShiftFilter =
    FFTShiftFilterType::New();
  fftShiftFilter->SetInput( gaussSource->GetOutput() );
  fftShiftFilter->Update();

  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  fftFilter->SetInput( fftShiftFilter->GetOutput() );
  fftFilter->Update();
  m_KernelImageFFT = fftFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::ComputeConvolvedImageFFT()
{
  std::cout << "FFTGaussianDerivativeIFFT: ComputeConvolvedImageFFT" << std::endl;
  typename MultiplyFilterType::Pointer multiplyFilter =
    MultiplyFilterType::New();
  multiplyFilter->SetInput1( m_InputImageFFT );
  multiplyFilter->SetInput2( m_KernelImageFFT );
  multiplyFilter->Update();

  m_ConvolvedImageFFT = multiplyFilter->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::ComputeConvolvedImage()
{
  std::cout << "FFTGaussianDerivativeIFFT: ComputeConvolvedImage" << std::endl;

  typename InverseFFTFilterType::Pointer 
    iFFTFilter = InverseFFTFilterType::New();
  iFFTFilter->SetInput( m_ConvolvedImageFFT );
  iFFTFilter->Update();

  typedef RegionFromReferenceImageFilter< RealImageType, TOutputImage >
    RegionFromFilterType;
  typename RegionFromFilterType::Pointer regionFrom = 
    RegionFromFilterType::New();

  regionFrom->SetInput1( iFFTFilter->GetOutput() );
  regionFrom->SetInput2( this->GetInput() );
  regionFrom->Update();

  m_ConvolvedImage = regionFrom->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::GenerateData()
{
  if( m_LastInputImage != this->GetInput() )
    {
    m_LastInputImage = this->GetInput();
    ComputeInputImageFFT();
    }

  ComputeKernelImageFFT();
  ComputeConvolvedImageFFT();
  ComputeConvolvedImage();

  this->SetNthOutput( 0, m_ConvolvedImage );
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::GenerateNJet( typename TOutputImage::Pointer & D,
  std::vector< typename TOutputImage::Pointer > & dX,
  std::vector< typename TOutputImage::Pointer > & dXX )
{
  if( m_LastInputImage != this->GetInput() )
    {
    m_LastInputImage = this->GetInput();
    ComputeInputImageFFT();
    }

  typedef RegionFromReferenceImageFilter< RealImageType, TOutputImage >
    RegionFromFilterType;

  if( dX.size() != ImageDimension )
    {
    dX.resize( ImageDimension );
    }
  std::vector< typename ComplexImageType::Pointer > 
    dXKernelImageFFT( ImageDimension );

  if( dXX.size() != ImageDimension * (ImageDimension-1) )
    {
    dXX.resize( ImageDimension * (ImageDimension-1) );
    }

  this->m_Orders.Fill( 0 );
  this->ComputeKernelImageFFT();
  this->ComputeConvolvedImageFFT();
  this->ComputeConvolvedImage();
  D = m_ConvolvedImage;

  for( unsigned int i = 0; i<ImageDimension; ++i )
    {
    this->m_Orders[i] = 1;
    this->ComputeKernelImageFFT();
    dXKernelImageFFT[i] = m_KernelImageFFT;
    this->ComputeConvolvedImageFFT();
    this->ComputeConvolvedImage();
    dX[i] = m_ConvolvedImage;
    this->m_Orders[i] = 0;
    }

  typename ComplexImageType::Pointer tmpInputImageFFT = m_InputImageFFT;
  typename ComplexImageType::Pointer tmpFirstConvolutionFFT;
  unsigned int count = 0;
  for( unsigned int i = 0; i<ImageDimension; ++i )
    {
    m_InputImageFFT = tmpInputImageFFT;
    m_KernelImageFFT = dXKernelImageFFT[i];
    this->ComputeConvolvedImageFFT();
    tmpFirstConvolutionFFT = m_ConvolvedImageFFT;
    
    for( unsigned int j = i; j<ImageDimension; ++j )
      {
      m_InputImageFFT = tmpFirstConvolutionFFT;
      m_KernelImageFFT = dXKernelImageFFT[j];
      this->ComputeConvolvedImageFFT();
      this->ComputeConvolvedImage();
      dXX[ count++ ] = m_ConvolvedImage;
      this->m_Orders[i] = 0;
      this->m_Orders[j] = 0;
      }
    }

  m_InputImageFFT = tmpInputImageFFT;

  this->SetNthOutput( 0, D );
}

template< typename TInputImage, typename TOutputImage >
void
FFTGaussianDerivativeIFFTFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_InputImageFFT.IsNotNull() )
    {
    os << indent << "FFT Image           : " << m_InputImageFFT << std::endl;
    }
  else
    {
    os << indent << "FFT Image           : NULL" << std::endl;
    }
  if( m_KernelImageFFT.IsNotNull() )
    {
    os << indent << "Kernel Image        : " << m_KernelImageFFT << std::endl;
    }
  else
    {
    os << indent << "Kernel Image        : NULL" << std::endl;
    }
  if( m_ConvolvedImageFFT.IsNotNull() )
    {
    os << indent << "Convolved Image FFT : " << m_ConvolvedImageFFT << std::endl;
    }
  else
    {
    os << indent << "Convolved Image FFT : NULL" << std::endl;
    }
  if( m_ConvolvedImage.IsNotNull() )
    {
    os << indent << "Convolved Image   : " << m_ConvolvedImage << std::endl;
    }
  else
    {
    os << indent << "Convolved Image   : NULL" << std::endl;
    }
  os << indent << "Orders              : " << m_Orders << std::endl;
  os << indent << "Sigmas               : " << m_Sigmas << std::endl;
  os << indent << "Last Input Image    : " << m_LastInputImage << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif
