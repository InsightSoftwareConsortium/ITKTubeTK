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
#ifndef __itktubeRidgeFFTFilter_h
#define __itktubeRidgeFFTFilter_h

#include "itktubeFFTGaussianDerivativeIFFTFilter.h"

#include "itkImageToImageFilter.h"

namespace itk
{

namespace tube
{

template< typename TInputImage >
class RidgeFFTFilter :
  public ImageToImageFilter< TInputImage, Image< float,
    TInputImage::ImageDimension > >
{
public:

  typedef TInputImage                                       InputImageType;

  typedef Image< float, TInputImage::ImageDimension >       OutputImageType;

  typedef RidgeFFTFilter                                     Self;
  typedef ImageToImageFilter< TInputImage, OutputImageType > Superclass;
  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( RidgeFFTFilter, ImageToImageFilter );

  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  itkSetMacro( Scale, double );
  itkGetMacro( Scale, double );

  itkSetMacro( CurvatureExpectedMax, double );
  itkGetMacro( CurvatureExpectedMax, double );

  itkGetConstReferenceMacro( Intensity, typename OutputImageType::Pointer );
  itkGetConstReferenceMacro( Ridgeness, typename OutputImageType::Pointer );
  itkGetConstReferenceMacro( Curvature, typename OutputImageType::Pointer );
  itkGetConstReferenceMacro( Levelness, typename OutputImageType::Pointer );
  itkGetConstReferenceMacro( Roundness, typename OutputImageType::Pointer );

protected:
  RidgeFFTFilter( void );
  virtual ~RidgeFFTFilter( void ) {}

  void GenerateData();

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:
  // Purposely not implemented
  RidgeFFTFilter( const Self & );
  void operator = ( const Self & );

  typename OutputImageType::Pointer                     m_Intensity;
  typename OutputImageType::Pointer                     m_Ridgeness;
  typename OutputImageType::Pointer                     m_Curvature;
  typename OutputImageType::Pointer                     m_Levelness;
  typename OutputImageType::Pointer                     m_Roundness;

  double                                                m_Scale;

  double                                                m_CurvatureExpectedMax;

  const InputImageType *                                m_LastInputImage;
};


// End class RidgeFFTFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeFFTFilter.hxx"
#endif

#endif
