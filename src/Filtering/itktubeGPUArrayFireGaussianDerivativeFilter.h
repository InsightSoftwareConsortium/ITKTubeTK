/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0( the "License" );
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
#ifndef __itktubeGPUArrayFireGaussianDerivativeFilter_h
#define __itktubeGPUArrayFireGaussianDerivativeFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkParametricImageSource.h"
#include "itkSize.h"

#include "itktubeGaussianDerivativeFilter.h"
#include "itkFFTShiftImageFilter.h"

// includes, project
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <helper_cuda.h>
#include <helper_functions.h>

namespace itk
{
namespace tube
{

template< typename TInputImage, typename TOutputImage =
Image< float, TInputImage::ImageDimension >  >
class GPUArrayFireGaussianDerivativeFilter :
  public GaussianDerivativeFilter< TInputImage, TOutputImage >
{
public:

  typedef GPUArrayFireGaussianDerivativeFilter                  Self;
  typedef GaussianDerivativeFilter< TInputImage,
          TOutputImage >                                        Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( GPUArrayFireGaussianDerivativeFilter,
                 GaussianDerivativeFilter );

  itkStaticConstMacro( ImageDimension, unsigned int,
                        TInputImage::ImageDimension );

  typedef TInputImage                      InputImageType;
  typedef TOutputImage                     OutputImageType;
  typedef float                            RealImagePixelType;
  typedef Image< float, ImageDimension >   RealImageType;
  typedef typename RealImageType::Pointer  RealImagePointerType;
  typedef GaussianDerivativeImageSource< RealImageType >
                                           GaussianDerivativeImageSourceType;
  typedef typename Superclass::OrdersType  OrdersType;
  typedef typename Superclass::SigmasType  SigmasType;

  void GenerateNJet( typename OutputImageType::Pointer & D,
                      std::vector< typename TOutputImage::Pointer > & Dx,
                      std::vector< typename TOutputImage::Pointer > & Dxx );

protected:

  typedef FFTShiftImageFilter< RealImageType, RealImageType >
  FFTShiftFilterType;

  GPUArrayFireGaussianDerivativeFilter( void );
  virtual ~GPUArrayFireGaussianDerivativeFilter( void ) {}

  void ComputeInputImageFFT();
  void RestoreInputImageFFT();
  void ComputeKernelImageFFT();
  void ComputeConvolvedImageFFT();
  void ComputeConvolvedImage();


  void GenerateData() override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  // Purposely not implemented
  GPUArrayFireGaussianDerivativeFilter( const Self & );
  void operator =( const Self & );

  typename RealImageType::Pointer     m_PaddedInputImage;
  Complex                           * m_ComplexPaddedInputImage;
  Complex                           * m_OnGPUPaddedInputImage;

  typename RealImageType::Pointer     m_PaddedKernelImage;
  Complex                           * m_ComplexPaddedKernelImage;
  Complex                           * m_OnGPUPaddedKernelImage;

  typename TOutputImage::Pointer      m_ConvolvedImage;

  unsigned int                        m_PaddedTotalSize;
  cufftHandle                         m_Plan;

  typename GaussianDerivativeImageSourceType::Pointer
                                      m_GaussSource;


  const InputImageType              * m_LastInputImage;
};

// Complex addition
static __device__ __host__ inline Complex
GPUFFTComplexAdd( Complex a, Complex b )
{
  Complex c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  return c;
}

// Complex scale
static __device__ __host__ inline Complex
GPUFFTComplexScale( Complex a, float s )
{
  Complex c;
  c.x = s * a.x;
  c.y = s * a.y;
  return c;
}

// Complex multiplication
static __device__ __host__ inline Complex
GPUFFTComplexMul( Complex a, Complex b )
{
  Complex c;
  c.x = a.x * b.x - a.y * b.y;
  c.y = a.x * b.y + a.y * b.x;
  return c;
}

// Complex pointwise multiplication
static __global__ void GPUFFTComplexPointwiseMulAndScale(
  Complex *a, const Complex *b, const Complex *c,
  int size, float scale)
{
  const int numThreads = blockDim.x * gridDim.x;
  const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

  for (int i = threadID; i < size; i += numThreads)
  {
    a[i] = GPUFFTComplexScale(GPUFFTComplexMul(b[i], c[i]), scale);
  }
}
// End class GPUArrayFireGaussianDerivativeFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeGPUArrayFireGaussianDerivativeFilter.cu"
#endif

#endif
