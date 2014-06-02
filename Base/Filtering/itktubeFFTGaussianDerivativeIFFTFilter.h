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
#ifndef __itktubeFFTGaussianDerivativeIFFTFilter_h
#define __itktubeFFTGaussianDerivativeIFFTFilter_h

#include "itktubeGaussianDerivativeImageSource.h"

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

template< typename TInputImage, typename TOutputImage = TInputImage >
class FFTGaussianDerivativeIFFTFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  typedef FFTGaussianDerivativeIFFTFilter                 Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( FFTGaussianDerivativeIFFTFilter, ImageToImageFilter );

  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  typedef TInputImage                         InputImageType;

  typedef TOutputImage                        OutputImageType;

  typedef Image< float, ImageDimension >      RealImageType;

  typedef GaussianDerivativeImageSource< RealImageType >
                                         GaussianDerivativeImageSourceType;

  typedef typename GaussianDerivativeImageSourceType::OrdersType
                                         OrdersType;

  typedef typename GaussianDerivativeImageSourceType::SigmasType
                                         SigmasType;

  itkSetMacro( Orders, OrdersType );
  itkGetConstReferenceMacro( Orders, OrdersType );

  itkSetMacro( Sigmas, SigmasType );
  itkGetConstReferenceMacro( Sigmas, SigmasType );

protected:
  typedef itk::ForwardFFTImageFilter< InputImageType > FFTInputType;
  typedef typename FFTInputType::OutputImageType       ComplexImageType;

  typedef itk::ForwardFFTImageFilter< RealImageType, ComplexImageType >
                                                       FFTRealType;

  typedef itk::FFTShiftImageFilter< RealImageType, RealImageType >
                                                       FFTShiftFilterType;

  typedef itk::InverseFFTImageFilter< ComplexImageType,
    OutputImageType >                                  InverseFFTFilterType;

  typedef itk::MultiplyImageFilter< ComplexImageType, ComplexImageType,
    ComplexImageType >                                 MultiplyFilterType;

  FFTGaussianDerivativeIFFTFilter( void );
  virtual ~FFTGaussianDerivativeIFFTFilter( void ) {}

  void ApplyFFT();

  void CreateGaussianDerivative();

  void ApplyGaussianDerivativeIFFT();

  void GenerateData();

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:
  // Purposely not implemented
  FFTGaussianDerivativeIFFTFilter( const Self & );
  void operator = ( const Self & );

  typename ComplexImageType::Pointer                  m_FFTImage;

  typename OutputImageType::Pointer                   m_IFFTImage;

  typename RealImageType::Pointer                     m_KernelImage;

  OrdersType                                          m_Orders;

  SigmasType                                          m_Sigmas;

  const InputImageType *                              m_LastInputImage;
};


// End class FFTGaussianDerivativeIFFTFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeFFTGaussianDerivativeIFFTFilter.hxx"
#endif

#endif
