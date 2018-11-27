/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeFFTGaussianDerivativeIFFTFilter_h
#define __itktubeFFTGaussianDerivativeIFFTFilter_h

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

  typedef FFTGaussianDerivativeIFFTFilter               Self;
  typedef GaussianDerivativeFilter< TInputImage,
    TOutputImage >                                      Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( FFTGaussianDerivativeIFFTFilter, GaussianDerivativeFilter );

  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  typedef TInputImage                         InputImageType;
  typedef TOutputImage                        OutputImageType;
  typedef Image< double, ImageDimension >     RealImageType;
  typedef typename RealImageType::Pointer     RealImagePointerType;
  typedef typename Superclass::GaussianDerivativeImageSourceType
                                              GaussianDerivativeImageSourceType;
  typedef typename Superclass::OrdersType     OrdersType;
  typedef typename Superclass::SigmasType     SigmasType;

  void GenerateNJet( typename OutputImageType::Pointer & D,
    std::vector< typename TOutputImage::Pointer > & Dx,
    std::vector< typename TOutputImage::Pointer > & Dxx );

protected:
  typedef ForwardFFTImageFilter< RealImageType >          FFTFilterType;

  typedef typename FFTFilterType::OutputImageType         ComplexImageType;

  typedef FFTShiftImageFilter< RealImageType, RealImageType >
                                                          FFTShiftFilterType;

  typedef InverseFFTImageFilter< ComplexImageType, RealImageType >
                                                          InverseFFTFilterType;

  typedef MultiplyImageFilter< ComplexImageType, ComplexImageType,
    ComplexImageType >                                    MultiplyFilterType;

  FFTGaussianDerivativeIFFTFilter( void );
  virtual ~FFTGaussianDerivativeIFFTFilter( void ) {}

  void ComputeInputImageFFT();
  void ComputeKernelImageFFT();
  void ComputeConvolvedImageFFT();
  void ComputeConvolvedImage();

  void GenerateData();

  void PrintSelf( std::ostream & os, Indent indent ) const;

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
