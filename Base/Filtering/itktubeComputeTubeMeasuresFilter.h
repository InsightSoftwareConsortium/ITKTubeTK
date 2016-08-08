/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeComputeTubeMeasuresFilter_h
#define __itktubeComputeTubeMeasuresFilter_h

#include <itkImageToImageFilter.h>
#include <itktubeRidgeFFTFilter.h>
#include <itkRescaleIntensityImageFilter.h>

namespace itk
{

namespace tube
{

/** \class ComputeTubeMeasuresFilter
 */

template< class TPixel, unsigned int Dimension >
class ComputeTubeMeasuresFilter
  : public ImageToImageFilter< Image< TPixel, Dimension >,
  Image< float, Dimension > >
{
public:

  /** Tube class typedef */
  typedef Image< TPixel, Dimension >                     InputImageType;
  typedef Image< float, Dimension >                      OutputImageType;

  /** Standard class typedefs. */
  typedef ComputeTubeMeasuresFilter                      Self;
  typedef ImageToImageFilter
    < InputImageType, OutputImageType >                  SuperClass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  typedef itk::RescaleIntensityImageFilter
    < InputImageType, OutputImageType >                  RescaleFilterType;
  typedef itk::tube::RidgeFFTFilter< OutputImageType >   RidgeFilterType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeTubeMeasuresFilter,
                ImageToImageFilter );

  /** Set/Get scale */
  itkSetMacro( Scale, int );
  itkGetMacro( Scale, int );

  /** Set/Get input Image */
  itkSetConstObjectMacro( InputImage, InputImageType );
  itkGetConstObjectMacro( InputImage, InputImageType );

  /** Get output Ridge Image */
  itkGetObjectMacro( Ridgeness, OutputImageType );

  /** Get output Round Image */
  itkGetObjectMacro( Roundness, OutputImageType );

  /** Get output Curve Image */
  itkGetObjectMacro( Curvature, OutputImageType );

  /** Get output Level Image */
  itkGetObjectMacro( Levelness, OutputImageType );

protected:

  ComputeTubeMeasuresFilter( void );
  ~ComputeTubeMeasuresFilter( void ) {};

  void GenerateData( void );
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  int                                         m_Scale;
  typename InputImageType::ConstPointer       m_InputImage;
  typename OutputImageType::Pointer           m_Ridgeness;
  typename OutputImageType::Pointer           m_Roundness;
  typename OutputImageType::Pointer           m_Curvature;
  typename OutputImageType::Pointer           m_Levelness;
}; // End class ComputeTubeMeasuresFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeTubeMeasuresFilter.hxx"
#endif

#endif // End !defined(__itktubeComputeTubeMeasuresFilter_h)
