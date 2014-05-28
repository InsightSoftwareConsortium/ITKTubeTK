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

#include "itkFFTShiftImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkInverseFFTImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkParametricImageSource.h"
#include "itkFixedArray.h"
#include "itkSize.h"

namespace itk
{

namespace tube
{

template< typename TInputImage, typename TOutputImage >
class FFTGaussianDerivativeIFFTFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef FFTGaussianDerivativeIFFTFilter               Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(FFTGaussianDerivativeIFFTFilter, ImageToImageFilter);
  typedef typename TInputImage::PixelType                InputPixelType;
  typedef typename TOutputImage::PixelType               OutputPixelType;
  typedef itk::Image< InputPixelType, 3 >                InputImageType;
  typedef itk::Image< OutputPixelType, 3 >               OutputImageTypeFilter;
  typedef itk::Image< unsigned char, 3 >                 UnsignedCharImageType;
  typedef GaussianDerivativeImageSource<InputImageType>  GaussianDerivativeType;
  typedef itk::ForwardFFTImageFilter<InputImageType>     FFTType;
  typedef typename FFTType::OutputImageType              ComplexImageType;
  typedef itk::FFTShiftImageFilter< InputImageType,
  InputImageType > FFTShiftFilterType;
  typedef itk::InverseFFTImageFilter<ComplexImageType,
  OutputImageTypeFilter> InverseFFTType;
  typedef itk::MultiplyImageFilter< ComplexImageType,
  InputImageType, ComplexImageType > MultiplyFilterType;


  itkSetMacro(Orders, typename GaussianDerivativeType::VectorType);
  itkGetConstReferenceMacro(Orders, typename GaussianDerivativeType::VectorType);
  itkSetMacro(Sigma, typename GaussianDerivativeType::ArrayType);
  itkGetConstReferenceMacro(Sigma, typename GaussianDerivativeType::ArrayType);

protected:
  FFTGaussianDerivativeIFFTFilter();
  void applyFFT();
  void applyIFFT(ComplexImageType* image);
  void applyGaussianDerivativeFilterIFFT();
  void createGaussianDerivative();
  void GenerateData();

private:
  typename GaussianDerivativeType::Pointer      m_GaussianDerivative;
  typename FFTType::Pointer                     m_FFTFilter;
  typename InverseFFTType::Pointer              m_InverseFFTFilter;
  typename GaussianDerivativeType::VectorType   m_Orders;
  typename GaussianDerivativeType::ArrayType    m_Sigma;
};


// End class FFTGaussianDerivativeIFFTFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeFFTGaussianDerivativeIFFTFilter.hxx"
#endif

#endif
