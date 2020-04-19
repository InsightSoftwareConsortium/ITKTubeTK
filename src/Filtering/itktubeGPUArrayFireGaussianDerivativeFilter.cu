/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
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

namespace itk
{

namespace tube
{
//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::GPUArrayFireGaussianDerivativeFilter()
{
  m_PaddedInputImage = NULL;
  m_PaddedKernelImage = NULL;
  m_ConvolvedImage = NULL;
  m_PaddedTotalSize = 0;
  //m_Plan;

  m_GaussSource = GaussianDerivativeImageSourceType::New();

  m_ComplexPaddedInputImage = NULL;
  m_OnGPUPaddedInputImage = NULL;
  m_ComplexPaddedKernelImage = NULL;
  m_OnGPUPaddedKernelImage = NULL;

  m_LastInputImage = NULL;
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeInputImageFFT()
{
  typedef PadImageFilter< InputImageType, RealImageType >
  PadFilterType;
  typename PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput( this->GetInput() );
  padFilter->SetGreatestPrimeFactor( 5 );
  padFilter->SetPadMethod( PadFilterType::ZERO_FLUX_NEUMANN );
  padFilter->Update();

  m_PaddedInputImage = padFilter->GetOutput();

  typename InputImageType::SizeType paddedSize = 
    m_PaddedInputImage.GetSize();
  m_PaddedTotalSize = 1;
  for( unsigned int d=0; d<ImageDimension; ++d )
    {
    m_PaddedTotalSize *= paddedSize[d];
    }

  //
  // Allocate complex image
  //
  if( m_ComplexPaddedInputImage != NULL )
    {
    delete m_ComplexPaddedInputImage;
    }
  m_ComplexPaddedInputImage = reinterpret_cast<Complex *>(
    malloc(sizeof(Complex) * m_PaddedTotalSize) );
  itk::ImageRegionIterator<InputImageType> iter( m_PaddedInputImage,
    m_PaddedInputImage->GetLargestPossibleRegion() );
  unsigned int count = 0;
  while( !iter.IsAtEnd() )
    {
    m_ComplexPaddedInputImage[count].x = iter.Get();
    m_ComplexPaddedInputImage[count].y = 0;
    ++iter;
    ++count;
    }

  //
  // Allocate complex kernel
  //
  if( m_ComplexPaddedKernelImage != NULL )
    {
    delete m_ComplexPaddedKernelImage;
    }
  m_ComplexPaddedKernelImage = reinterpret_cast<Complex *>(
    malloc(sizeof(Complex) * m_PaddedTotalSize) );

  // 
  // Identify NVidia device
  //
  int deviceCount = 0;
  checkCudaErrors( cudaGetDeviceCount(&deviceCount) );

  int targetDevice = 0;
  targetDevice = findCudaDevice();
  cudaCheckErrors( cudaSetDevice(targetDevice) );

  //
  // Allocate image memory on GPU
  //
  int mem_size = sizeof(Complex) * m_PaddedTotalSize;

  if( m_OnGPUPaddedInputImage != NULL )
    {
    delete m_OnGPUPaddedInputImage;
    }
  checkCudaErrors( cudaMalloc(
    reinterpret_cast<void **>(&m_OnGPUPaddedInputImage), mem_size));

  //
  // Copy image to GPU
  //
  checkCudaErrors( cudaMemcpy(
    m_OnGPUPaddedInputImage, m_ComplexPaddedInputImage, mem_size,
    cudaMemcpyHostToDevice));

  //
  // Describe GPU workspace for image
  //
  switch( ImageDimension )
    {
    case 1:
      checkCudaErrors(cufftPlan1d(&m_Plan, inputSize[0],
          CUFFT_C2C, 1));
      break;
    case 2:
      checkCudaErrors(cufftPlan2d(&m_Plan, inputSize[1],
          inputSize[0], CUFFT_C2C));
      break;
    case 3:
      checkCudaErrors(cufftPlan3d(&m_Plan, inputSize[2],
          inputSize[1], inputSize[0], CUFFT_C2C));
      break;
    default:
      itk::InvalidArgumentError e( __FILE__, __LINE__ );
      e.SetDescription( "Only Dimensions up to 3 are supported" );
      e.SetLocation( "CudaFFT" );
      throw e;
    }
  
  //
  // Perform FFT of image
  //
  checkCudaErrors(cufftExecC2C(m_Plan,
    reinterpret_cast<cufftComplex *>(m_OnGPUPaddedInputImage),
    reinterpret_cast<cufftComplex *>(m_OnGPUPaddedInputImage),
    CUFFT_FORWARD));

  //
  // Setup Gaussian kernel
  //
  const typename RealImageType::RegionType inputRegion(
    this->GetInput()->GetLargestPossibleRegion() );
  const typename RealImageType::SizeType inputSize
    = inputRegion.GetSize();

  const typename RealImageType::RegionType fftRegion(
    m_PaddedInputImage->GetLargestPossibleRegion() );
  const typename RealImageType::SizeType fftSize
    = fftRegion.GetSize();
  const typename RealImageType::SpacingType fftSpacing =
    m_PaddedInputImage->GetSpacing();
  const typename RealImageType::PointType fftOrigin =
    m_PaddedInputImage->GetOrigin();
  const typename RealImageType::DirectionType fftDirection =
    m_PaddedInputImage->GetDirection();

  m_GaussSource->SetIndex( fftRegion.GetIndex() );
  m_GaussSource->SetSize( fftSize );
  m_GaussSource->SetSpacing( fftSpacing );
  m_GaussSource->SetOrigin( fftOrigin );
  m_GaussSource->SetDirection( fftDirection );

  typename GaussianDerivativeImageSourceType::PointType mean;
  typename GaussianDerivativeImageSourceType::IndexType meanIndex;
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    const int halfLength = ( inputSize[ii]  / 2.0 );
    meanIndex[ii] = inputRegion.GetIndex() [ii] + halfLength;
    }
  this->GetInput()->TransformIndexToPhysicalPoint( meanIndex,
    mean );
  m_GaussSource->SetMean( mean );
}


template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeKernelImageFFT()
{
  // 
  // Define kernel
  //
  m_GaussSource->SetSigmas( this->m_Sigmas );
  m_GaussSource->SetOrders( this->m_Orders );
  m_GaussSource->Update();

  typename FFTShiftFilterType::Pointer fftShiftFilter =
    FFTShiftFilterType::New();
  fftShiftFilter->SetInput( m_GaussSource->GetOutput() );
  fftShiftFilter->Update();

  m_PaddedKernelImage = fftShiftFilter->GetOutput();

  itk::ImageRegionIterator<InputImageType> iter(
    m_PaddedKernelImage,
    m_PaddedKernelImage->GetLargestPossibleRegion() );
  unsigned int count = 0;
  while( !iter.IsAtEnd() )
    {
    m_ComplexPaddedKernelImage[count].x = iter.Get();
    m_ComplexPaddedKernelImage[count].y = 0;
    ++iter;
    ++count;
    }

  //
  // Identify NVidia device
  //
  int deviceCount = 0;
  checkCudaErrors( cudaGetDeviceCount(&deviceCount) );

  int targetDevice = 0;
  targetDevice = findCudaDevice();
  cudaCheckErrors( cudaSetDevice(targetDevice) );

  //
  // Allocate GPU memory for kernel and copy kernel to GPU
  //
  int mem_size = sizeof(Complex) * m_PaddedTotalSize;
  if( m_OnGPUPaddedKernelImage == NULL )
    {
    checkCudaErrors( cudaMalloc(
      reinterpret_cast<void **>(&m_OnGPUPaddedKernelImage),
      mem_size));
    }
  checkCudaErrors( cudaMemcpy(
    m_OnGPUPaddedKernelImage, m_ComplexPaddedKernelImage, mem_size,
    cudaMemcpyHostToDevice));

  //
  // Describe GPU workspace for kernel
  //
  cufftHandle plan_adv;
  size_t workSize;
  long long int new_size_long[3];
  for( unsigned int d=0; d<ImageDimension; ++d)
    {
    new_size_long[d] = paddedSize[d];
    }
  checkCudaErrors(cufftCreate(&plan_adv));
  switch( ImageDimension )
    {
    case 1:
      checkCudaErrors(cufftXtMakePlanMany(plan_adv, 1,
          new_size_long, NULL, 1, 1, CUDA_C_32F,
          NULL, 1, 1, CUDA_C_32F, 1, &workSize, CUDA_C_32F));
      break;
    case 2:
      checkCudaErrors(cufftXtMakePlanMany(plan_adv, 2,
          &new_size_long, NULL, 1, 1, CUDA_C_32F,
          NULL, 1, 1, CUDA_C_32F, 1, &workSize, CUDA_C_32F));
      break;
    case 3:
      checkCudaErrors(cufftXtMakePlanMany(plan_adv, 3,
          &new_size_long, NULL, 1, 1, CUDA_C_32F,
          NULL, 1, 1, CUDA_C_32F, 1, &workSize, CUDA_C_32F));
      break;
    default:
      itk::InvalidArgumentError e( __FILE__, __LINE__ );
      e.SetDescription( "Only Dimensions up to 3 are supported" );
      e.SetLocation( "CudaFFT" );
      throw e;
    }

  //
  // Perform FFT of kernel
  //
  checkCudaErrors(cufftExecC2C(plan_adv,
    reinterpret_cast<cufftComplex *>(m_OnGPUPaddedKernelImage),
    reinterpret_cast<cufftComplex *>(m_OnGPUPaddedKernelImage),
    CUFFT_FORWARD));
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeConvolvedImageFFT()
{
  //
  // Allocate GPU memory for convolution results and init memory
  //
  if( m_OnGPUPaddedConvolutionImage == NULL )
    {
    checkCudaErrors( cudaMalloc(
      reinterpret_cast<void **>(&m_OnGPUPaddedConvolutionImage),
      mem_size));
    }

  //
  // Compute convolution
  //
  GPUFFTComplexPointwiseMulAndScale<<<32, 256>>>(
    m_OnGPUPaddedConvolutionImage, m_OnGPUPaddedInputImage,
    m_OnGPUPaddedKernelImage, m_PaddedTotalSize,
    1.0f / m_PaddedTotalSize);
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::ComputeConvolvedImage()
{
  checkCudaErrors(cufftExecC2C( m_Plan,
    reinterpret_cast<cufftComplex *>(m_OnGPUPaddedConvolutionImage),
    reinterpret_cast<cufftComplex *>(m_OnGPUPaddedConvolutionImage),
    CUFFT_INVERSE));

  Complex * tmpComplexConvolutionImage =
    reinterpret_cast<Complex *>(
      malloc(sizeof(Complex) * m_PaddedTotalSize) );
  checkCudaErrors( cudaMemcpy(
    tmpComplexConvolutionImage, m_OnGPUPaddedConvolutionImage,
    mem_size, cudaMemcpyDeviceToHost));

  m_ConvolvedImage = TOutputImage::New();
  m_ConvolvedImage->SetRegions(
    m_InputImage->GetLargestPossibleRegion() );
  m_ConvolvedImage->CopyInformation( m_InputImage );
  m_ConvolvedImage->Allocate();
  itk::ImageRegionIterator<OutputImageType> iterPad(
    m_PaddedInputImage,
    m_PaddedInputImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<OutputImageType> iter(
    m_ConvolvedImage,
    m_ConvolvedImage->GetLargestPossibleRegion() );
  unsigned int count = 0;
  while( !iter.IsAtEnd() )
    {
    if( iterPad.GetIndex() == iter.GetIndex() )
      {
      iter.Set( tmpComplexConvolutionImage[count].x );
      ++iter;
      }
    ++iterPad;
    ++count;
    }

  delete tmpComplexConvolutionImage;
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
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
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::GenerateNJet( typename TOutputImage::Pointer & D,
  std::vector< typename TOutputImage::Pointer > & dX,
  std::vector< typename TOutputImage::Pointer > & dXX )
{
  if( m_LastInputImage != this->GetInput() )
    {
    m_LastInputImage = this->GetInput();
    ComputeInputImageFFT();
    }

  if( dX.size() != ImageDimension )
    {
    dX.resize( ImageDimension );
    }

  unsigned int ddxSize = 0;
  for( unsigned int i = 1; i<=ImageDimension; ++i )
    {
    ddxSize += i;
    }
  if( dXX.size() != ddxSize )
    {
    dXX.resize( ddxSize );
    }

  // Compute G_0
  this->m_Orders.Fill( 0 );

  this->ComputeKernelImageFFT();
  this->ComputeConvolvedImageFFT();
  this->ComputeConvolvedImage();
  D = m_ConvolvedImage;

  // Compute G_1
  std::vector<Complex *> dXKernelImageFFT( ImageDimension );
  for( unsigned int i = 0; i<ImageDimension; ++i )
    {
    this->m_Orders[i] = 1;
    this->ComputeKernelImageFFT();
    this->ComputeConvolvedImageFFT();
    this->ComputeConvolvedImage();
    dXKernelImageFFT[i] = m_OnGPUPaddedKernelImage;
    m_OnGPUPaddedKernelImage = NULL;
    dX[i] = m_ConvolvedImage;
    this->m_Orders[i] = 0;
    }

  // compute G_2
  unsigned int count = 0;
  for( unsigned int i = 0; i<ImageDimension; ++i )
    {
    m_OnGPUPaddedKernelImage = dXKernelImageFFT[i];
    this->ComputeConvolvedImageFFT();
    Complex * tmpConvolvedFFT = m_OnGPUPaddedConvolvedImage;
    Complex * tmpInputFFT = m_OnGPUPaddedInputImage;
    m_OnGPUPaddedConvolvedImage = NULL;
    this->m_Orders[i] = 1;
    for( unsigned int j = i; j<ImageDimension; ++j )
      {
      this->m_Orders[j] = 1;
      m_OnGPUPaddedInputImage = tmpConvolvedFFT;
      m_OnGPUPaddedKernelImage = dXKernelImageFFT[j];
      this->ComputeConvolvedImageFFT();
      this->ComputeConvolvedImage();
      dXX[ count++ ] = m_ConvolvedImage;
      this->m_Orders[j] = 0;
      }
    cudaFree( tmpConvolvedFFT );
    m_OnGPUPaddedInputImage = tmpInputFFT;
    this->m_Orders[i] = 0;
    }

  delete * dXKernelImgeFFT;

  this->SetNthOutput( 0, D );
}

template< typename TInputImage, typename TOutputImage >
void
GPUArrayFireGaussianDerivativeFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_ConvolvedImage.IsNotNull() )
    {
    os << indent << "Convolved Image   : " << m_ConvolvedImage << std::endl;
    }
  else
    {
    os << indent << "Convolved Image   : NULL" << std::endl;
    }

  if( m_LastInputImage != NULL )
    {
    os << indent << "Last Input Image    : Initialized" << std::endl;
    }
  else
    {
    os << indent << "Last Input Image    : NULL" << std::endl;
    }
}

} // End namespace tube

} // End namespace itk

#endif
