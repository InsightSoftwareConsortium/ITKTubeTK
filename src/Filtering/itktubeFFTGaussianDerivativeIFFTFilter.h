/*=========================================================================

Library:   TubeTK

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

#ifndef itktubeFFTGaussianDerivativeIFFTFilter_h
#define itktubeFFTGaussianDerivativeIFFTFilter_h

#include "itktubeGaussianDerivativeFilter.h"

#include "itkFFTShiftImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkInverseFFTImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkParametricImageSource.h"
#include "itkSize.h"

namespace itk
{
namespace tube
{

template< typename TInputImage, typename TOutputImage =
  Image< float, TInputImage::ImageDimension >  >
class FFTGaussianDerivativeIFFTFilter :
    public GaussianDerivativeFilter< TInputImage, TOutputImage >
{
public:

  using Self = FFTGaussianDerivativeIFFTFilter;
  using Superclass = GaussianDerivativeFilter< TInputImage,
    TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  itkNewMacro( Self );

  itkTypeMacro( FFTGaussianDerivativeIFFTFilter, GaussianDerivativeFilter );

  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using RealImageType = Image< double, ImageDimension >;
  using RealImagePointerType = typename RealImageType::Pointer;
  typedef typename Superclass::GaussianDerivativeImageSourceType
                                              GaussianDerivativeImageSourceType;
  using OrdersType = typename Superclass::OrdersType;
  using SigmasType = typename Superclass::SigmasType;

  void GenerateNJet( typename OutputImageType::Pointer & D,
    std::vector< typename TOutputImage::Pointer > & Dx,
    std::vector< typename TOutputImage::Pointer > & Dxx ) override;

  void GenerateData() override;

protected:
  using FFTFilterType = ForwardFFTImageFilter< RealImageType >;

  using ComplexImageType = typename FFTFilterType::OutputImageType;

  using FFTShiftFilterType = FFTShiftImageFilter< RealImageType, RealImageType >;

  using InverseFFTFilterType = InverseFFTImageFilter< ComplexImageType, RealImageType >;

  using MultiplyFilterType = MultiplyImageFilter< ComplexImageType, ComplexImageType,
    ComplexImageType >;

  FFTGaussianDerivativeIFFTFilter( void );
  virtual ~FFTGaussianDerivativeIFFTFilter( void ) {}

  void ComputeInputImageFFT();
  void ComputeKernelImageFFT();
  void ComputeConvolvedImageFFT();
  void ComputeConvolvedImage();


  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  // Purposely not implemented
  FFTGaussianDerivativeIFFTFilter( const Self & );
  void operator = ( const Self & );

  typename ComplexImageType::Pointer                  m_InputImageFFT;

  typename ComplexImageType::Pointer                  m_KernelImageFFT;

  typename ComplexImageType::Pointer                  m_ConvolvedImageFFT;

  typename TOutputImage::Pointer                      m_ConvolvedImage;

  const InputImageType *                              m_LastInputImage;
};


// End class FFTGaussianDerivativeIFFTFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeFFTGaussianDerivativeIFFTFilter.hxx"
#endif

#endif
