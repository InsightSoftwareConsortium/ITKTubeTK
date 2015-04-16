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
#ifndef __itktubeGPUArrayFireGaussianDerivativeFilter_hxx
#define __itktubeGPUArrayFireGaussianDerivativeFilter_hxx

#include "itktubeGPUArrayFireGaussianDerivativeFilter.h"
#include "itktubeGaussianDerivativeImageSource.h"
#include "itktubePadImageFilter.h"
#include "itktubeRegionFromReferenceImageFilter.h"

#include "itkParametricImageSource.h"
#include "itkSize.h"

#include "itktubeArrayFireGlueUtilities.h"

namespace itk
{

namespace tube
{
//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::GPUArrayFireGaussianDerivativeFilter()
{
  m_ConvolvedImage = NULL;
  m_PaddedInputImage = NULL;
  this->m_LastInputImage = NULL;
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeInputImageFFT()
{
  // std::cout << "ArrayFireGPUGaussianDerivative: ComputeInputImageFFT" << std::endl;

  typedef PadImageFilter< InputImageType, RealImageType >
  PadFilterType;
  typename PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput ( this->GetInput() );
  padFilter->SetGreatestPrimeFactor ( 5 );
  padFilter->SetPadMethod ( PadFilterType::ZERO_FLUX_NEUMANN );
  padFilter->Update();

  m_PaddedInputImage = padFilter->GetOutput();

  af::array imPaddedInputAfArr;
  itk::convertITKImageToArrayFire< RealImageType > ( m_PaddedInputImage,
      imPaddedInputAfArr );

  switch ( imPaddedInputAfArr.numdims() )
    {
    case 1:

      m_InputImageFFTAfArr = af::fft ( imPaddedInputAfArr );
      break;

    case 2:

      m_InputImageFFTAfArr = af::fft2 ( imPaddedInputAfArr );
      break;

    case 3:

      m_InputImageFFTAfArr = af::fft3 ( imPaddedInputAfArr );
      break;

    default:

      itk::InvalidArgumentError e ( __FILE__, __LINE__ );
      e.SetDescription ( "Only Dimensions upto 3 is not supported" );
      e.SetLocation ( "convertITKImageToArrayFire" );
      throw e;
    }

}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeKernelImageFFT()
{
  // std::cout << "ArrayFireGPUGaussianDerivative: ComputeKernelImageFFT" << std::endl;

  typename GaussianDerivativeImageSourceType::Pointer gaussSource =
    GaussianDerivativeImageSourceType::New();

  const typename RealImageType::RegionType inputRegion (
    this->GetInput()->GetLargestPossibleRegion() );
  const typename RealImageType::SizeType inputSize
    = inputRegion.GetSize();

  const typename RealImageType::RegionType fftRegion (
    m_PaddedInputImage->GetLargestPossibleRegion() );
  const typename RealImageType::SizeType fftSize
    = fftRegion.GetSize();
  const typename RealImageType::SpacingType fftSpacing =
    m_PaddedInputImage->GetSpacing();
  const typename RealImageType::PointType fftOrigin =
    m_PaddedInputImage->GetOrigin();
  const typename RealImageType::DirectionType fftDirection =
    m_PaddedInputImage->GetDirection();

  gaussSource->SetIndex ( fftRegion.GetIndex() );
  gaussSource->SetSize ( fftSize );
  gaussSource->SetSpacing ( fftSpacing );
  gaussSource->SetOrigin ( fftOrigin );
  gaussSource->SetDirection ( fftDirection );

  typename GaussianDerivativeImageSourceType::PointType mean;
  typename GaussianDerivativeImageSourceType::IndexType meanIndex;
  for ( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    const int halfLength = ( inputSize[ii]  / 2.0 );
    meanIndex[ii] = inputRegion.GetIndex() [ii] + halfLength;
    }
  this->GetInput()->TransformIndexToPhysicalPoint ( meanIndex, mean );

  gaussSource->SetSigmas ( this->m_Sigmas );
  gaussSource->SetMean ( mean );
  gaussSource->SetOrders ( this->m_Orders );
  gaussSource->Update();

  typename FFTShiftFilterType::Pointer fftShiftFilter =
    FFTShiftFilterType::New();
  fftShiftFilter->SetInput ( gaussSource->GetOutput() );
  fftShiftFilter->Update();

  af::array imPaddedKernelAfArr;
  itk::convertITKImageToArrayFire< RealImageType > ( fftShiftFilter->GetOutput(),
      imPaddedKernelAfArr );

  switch ( imPaddedKernelAfArr.numdims() )
    {
    case 1:

      m_KernelImageFFTAfArr = af::fft ( imPaddedKernelAfArr );
      break;

    case 2:

      m_KernelImageFFTAfArr = af::fft2 ( imPaddedKernelAfArr );
      break;

    case 3:

      m_KernelImageFFTAfArr = af::fft3 ( imPaddedKernelAfArr );
      break;

    default:

      itk::InvalidArgumentError e ( __FILE__, __LINE__ );
      e.SetDescription ( "Only Dimensions upto 3 are supported" );
      e.SetLocation ( "GPUArrayFireGaussianDerivativeFilter: "\
                      "ComputeKernelImageFFT" );
      throw e;
    }


  /* The following code constructs the gaussian kernel directly in k-space
   * and computes the derivative (at requested order) using an analytical
   * formula. Compared to the code above which computes the gauss-derivative
   * kernel on the cpu, transfers to the gpu, and applies fft to it -- the
   * following code saves on transfer of kernel from cpu to gpu and computing
   * its fft. Though this wont give a speed up in order of magnitude this will
   * be faster depending on the cpu-gpu data transfer rate.
   *
   * As of now this code has issues with some corner cases either due to
   * normalization or in the calculation of frequency values of each element in
   * the output of fft. I've posted it on stackoverflow and Arrayfire and will
   * fix this once i hear from them.
   *
  // std::cout << "\nSigmas = " << m_Sigmas << ", Orders = " << m_Orders << std::endl;

  const typename RealImageType::SizeType fftSize
     = m_PaddedInputImage->GetLargestPossibleRegion().GetSize();

  // std::cout << "FFT Size = " << fftSize << std::endl;

  const typename RealImageType::RegionType inputRegion(
    this->GetInput()->GetLargestPossibleRegion() );

  const typename RealImageType::SizeType inputSize
    = this->GetInput()->GetLargestPossibleRegion().GetSize();

  typename RealImageType::PointType mean;
  typename RealImageType::IndexType meanIndex;
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    const int halfLength = (inputSize[ii]  / 2.0);
    meanIndex[ii] = inputRegion.GetIndex()[ii] + halfLength;
    }
  this->GetInput()->TransformIndexToPhysicalPoint( meanIndex, mean );

  // compute output frequency list of FFT for each dimension
  const int DIMENSION = TInputImage::ImageDimension;

  af::array kdim[ DIMENSION ];
  int N[ DIMENSION ];

  for(unsigned int i = 0; i < DIMENSION; i++)
    {
    N[i] = fftSize[i];
    if( N[i] % 2 == 0 )
      {
      kdim[i] = af::join(0, af::seq(0, (N[i]/2)-1), af::seq(-(N[i]/2), -1));
     }
    else // Bug in seq? with seq(-1, -1) throws out of memory. Report it.
      {
      kdim[i] = af::join(0, af::seq(0, (N[i]-1)/2), af::seq(-(N[i]-1)/2, -1));
      }
    // af_print( kdim[i] );
    kdim[i] *= (2.0 * af::Pi / N[i]);
    }

  // build meshgrid of frequency values for each dimension
  af::array kgrid[  DIMENSION ];
  switch( DIMENSION )
    {
    case 1:

      kgrid[0] = kdim[0];
      break;

    case 2:

      kdim[1] = af::reorder(kdim[1], 1, 0);

      kgrid[0] = af::tile(kdim[0], 1, N[1]);
      kgrid[1] = af::tile(kdim[1], N[0], 1);

      break;

    case 3:

      kdim[1] = af::reorder(kdim[1], 2, 0, 1);
      kdim[2] = af::reorder(kdim[2], 1, 2, 0);

      kgrid[0] = af::tile(kdim[0], 1, N[1], N[2]);
      kgrid[1] = af::tile(kdim[1], N[0], 1, N[2]);
      kgrid[2] = af::tile(kdim[2], N[0], N[1], 1);

      break;

    default:

      itk::InvalidArgumentError e( __FILE__, __LINE__ );
      e.SetDescription( "Only Dimensions upto 3 are supported" );
      e.SetLocation( "GPUArrayFireGaussianDerivativeFilter: "\
         "ComputeKernelImageFFT" );
      throw e;
    }

  // build FFT of gaussian
  af::array gaussFFT = af::complex( af::constant(1.0, kgrid[0].dims()) );

  for(unsigned i = 0; i < DIMENSION; i++)
    {
    // exp( -k^2 sigma^2 )
    gaussFFT *= af::exp( -( 0.5 * kgrid[i] * kgrid[i]
            * m_Sigmas[i] * m_Sigmas[i] ) );
    }

  // af_print( gaussFFT );

  // build FFT of derivative of gaussian
  af::cdouble cj(0.0, 1.0);
  m_KernelImageFFTAfArr = gaussFFT;

  for(unsigned i = 0; i < DIMENSION; i++)
    {
    if( m_Orders[i] > 0 )
      {
      // (1i * k)^order * FFT
      m_KernelImageFFTAfArr *= std::pow(cj, m_Orders[i]) *
               af::complex( af::pow(kgrid[i], m_Orders[i]) );
      }
    }
   // af_print( m_KernelImageFFTAfArr );
   */
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeConvolvedImageFFT()
{
  // std::cout << "ArrayFireGPUGaussianDerivative: ComputeConvolvedImageFFT" << std::endl;

  m_ConvolvedImageFFTAfArr = m_InputImageFFTAfArr * m_KernelImageFFTAfArr;
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeConvolvedImage()
{
  // std::cout << "ArrayFireGPUGaussianDerivative: ComputeConvolvedImage" << std::endl;

  af::array imConvolvedImageAfArr;

  switch ( m_ConvolvedImageFFTAfArr.numdims() )
    {
    case 1:

      imConvolvedImageAfArr = af::real ( af::ifft ( m_ConvolvedImageFFTAfArr ) );
      break;

    case 2:

      imConvolvedImageAfArr = af::real ( af::ifft2 ( m_ConvolvedImageFFTAfArr ) );
      break;

    case 3:

      imConvolvedImageAfArr = af::real ( af::ifft3 ( m_ConvolvedImageFFTAfArr ) );
      break;

    default:

      itk::InvalidArgumentError e ( __FILE__, __LINE__ );
      e.SetDescription ( "Only Dimensions upto 3 is not supported" );
      e.SetLocation ( "GPUArrayFireGaussianDerivativeFilter: "\
                      "ComputeConvolvedImage" );
      throw e;
    }

  RealImagePointerType pPaddedConvolvedImage;
  itk::convertArrayFireImageToITK< RealImageType > ( imConvolvedImageAfArr,
      pPaddedConvolvedImage,
      m_PaddedInputImage );

  typedef RegionFromReferenceImageFilter< RealImageType, TOutputImage >
  RegionFromFilterType;
  typename RegionFromFilterType::Pointer regionFrom =
    RegionFromFilterType::New();

  regionFrom->SetInput1 ( pPaddedConvolvedImage );
  regionFrom->SetInput2 ( this->GetInput() );
  regionFrom->Update();

  m_ConvolvedImage = regionFrom->GetOutput();
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::GenerateData()
{
  if ( m_LastInputImage != this->GetInput() )
    {
    m_LastInputImage = this->GetInput();
    ComputeInputImageFFT();
    }

  ComputeKernelImageFFT();
  ComputeConvolvedImageFFT();
  ComputeConvolvedImage();

  this->SetNthOutput ( 0, m_ConvolvedImage );
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::GenerateNJet ( typename TOutputImage::Pointer & D,
                 std::vector< typename TOutputImage::Pointer > & dX,
                 std::vector< typename TOutputImage::Pointer > & dXX )
{
  if ( m_LastInputImage != this->GetInput() )
    {
    m_LastInputImage = this->GetInput();
    ComputeInputImageFFT();
    }

  typedef RegionFromReferenceImageFilter< RealImageType, TOutputImage >
  RegionFromFilterType;

  if ( dX.size() != ImageDimension )
    {
    dX.resize ( ImageDimension );
    }

  unsigned int ddxSize = 0;
  for ( unsigned int i = 1; i<=ImageDimension; ++i )
    {
    ddxSize += i;
    }
  if ( dXX.size() != ddxSize )
    {
    dXX.resize ( ddxSize );
    }

  // Compute G_0
  this->m_Orders.Fill ( 0 );

  this->ComputeKernelImageFFT();
  this->ComputeConvolvedImageFFT();
  this->ComputeConvolvedImage();
  D = m_ConvolvedImage;

  // Compute G_1
  af::array dXKernelImageFFT[ ImageDimension ];
  for ( unsigned int i = 0; i<ImageDimension; ++i )
    {
    this->m_Orders[i] = 1;
    this->ComputeKernelImageFFT();
    dXKernelImageFFT[i] = m_KernelImageFFTAfArr;
    this->ComputeConvolvedImageFFT();
    this->ComputeConvolvedImage();
    dX[i] = m_ConvolvedImage;
    this->m_Orders[i] = 0;
    }

  // compute G_2
  unsigned int count = 0;

  /* This will be used to compute second-order derivatives instead of the code
   * below when computation of kernel directly in frequency domain is resolved.
   *
  this->m_Orders.Fill( 0.0 );
  for( unsigned int i = 0; i<ImageDimension; ++i )
    {
    this->m_Orders[i] += 1;
    for( unsigned int j = i; j<ImageDimension; ++j )
      {
      this->m_Orders[j] += 1;
      this->ComputeKernelImageFFT();
      this->ComputeConvolvedImageFFT();
      this->ComputeConvolvedImage();
      dXX[ count++ ] = m_ConvolvedImage;
      this->m_Orders[j] -= 1;
      }
    this->m_Orders[i] -= 1;
    }
  */

  // compute G2
  af::array tmpInputImageFFTAfArr = m_InputImageFFTAfArr;
  af::array tmpFirstConvolutionFFTAfArr;

  for ( unsigned int i = 0; i<ImageDimension; ++i )
    {
    m_InputImageFFTAfArr = tmpInputImageFFTAfArr;
    m_KernelImageFFTAfArr = dXKernelImageFFT[i];
    this->ComputeConvolvedImageFFT();
    tmpFirstConvolutionFFTAfArr = m_ConvolvedImageFFTAfArr;

    for ( unsigned int j = i; j<ImageDimension; ++j )
      {
      m_InputImageFFTAfArr = tmpFirstConvolutionFFTAfArr;
      m_KernelImageFFTAfArr = dXKernelImageFFT[j];
      this->ComputeConvolvedImageFFT();
      this->ComputeConvolvedImage();
      dXX[ count++ ] = m_ConvolvedImage;
      this->m_Orders[i] = 0;
      this->m_Orders[j] = 0;
      }
    }
  m_InputImageFFTAfArr = tmpInputImageFFTAfArr;

  this->SetNthOutput ( 0, D );
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::PrintSelf ( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf ( os, indent );

  if ( m_ConvolvedImage.IsNotNull() )
    {
    os << indent << "Convolved Image   : " << m_ConvolvedImage << std::endl;
    }
  else
    {
    os << indent << "Convolved Image   : NULL" << std::endl;
    }

  os << indent << "Last Input Image    : " << m_LastInputImage << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif
