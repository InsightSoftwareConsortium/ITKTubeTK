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

#include <arrayfire.h>

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
  typedef Image< double, ImageDimension >  RealImageType;
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
  void ComputeKernelImageFFT();
  void ComputeConvolvedImageFFT();
  void ComputeConvolvedImage();

  void GenerateData() override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  // Purposely not implemented
  GPUArrayFireGaussianDerivativeFilter( const Self & );
  void operator =( const Self & );

  typename RealImageType::Pointer                     m_PaddedInputImage;

  af::array                                           m_InputImageFFTAfArr;

  af::array                                           m_GaussianFFTAfArr;

  af::array                                           m_KernelImageFFTAfArr;

  af::array                                           m_ConvolvedImageFFTAfArr;

  typename TOutputImage::Pointer                      m_ConvolvedImage;

  const InputImageType *                              m_LastInputImage;
};


// End class GPUArrayFireGaussianDerivativeFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeGPUArrayFireGaussianDerivativeFilter.hxx"
#endif

#endif
