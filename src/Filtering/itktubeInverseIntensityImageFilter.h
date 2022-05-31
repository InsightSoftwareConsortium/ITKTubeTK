/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeInverseIntensityImageFilter_h
#define __itktubeInverseIntensityImageFilter_h

#include <itkImageToImageFilter.h>

#ifndef Tdimension
#define Tdimension 3
#endif

namespace itk
{

namespace tube
{

template< class TInputImage >
class InverseIntensityImageFilter
  : public ImageToImageFilter< TInputImage, TInputImage >
{
public:

  /** Standard class typedefs. */
  typedef InverseIntensityImageFilter                     Self;
  typedef ImageToImageFilter< TInputImage, TInputImage>   SuperClass;

  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::PixelType      InputPixelType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;

  typedef TInputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType     OutputPixelType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;

  typedef typename InputImageType::RegionType     RegionType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  itkSetMacro( InverseMaximumIntensity, InputPixelType );
  itkGetMacro( InverseMaximumIntensity, InputPixelType );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubeNetworkSpatialObjectToImageFilter, ImageToImageFilter );

protected:

  InverseIntensityImageFilter( void );
  ~InverseIntensityImageFilter( void ) {}

  /** GenerateData produce the main work */
  void GenerateData( void ) override;

private:

  void PrintSelf( std::ostream& os, Indent indent ) const override
    { SuperClass::PrintSelf( os, indent );   }

  InputPixelType                 m_InverseMaximumIntensity;

}; // End class InverseIntensityImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeInverseIntensityImageFilter.hxx"
#endif

#endif // End !defined( __itktubeInverseIntensityImageFilter_h )
