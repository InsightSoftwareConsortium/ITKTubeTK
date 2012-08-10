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
#ifndef __itkLabelMapToAcousticImpedanceImageFilter_h
#define __itkLabelMapToAcousticImpedanceImageFilter_h

#include <itkUnaryFunctorImageFilter.h>
#include <itkLabelMapToAcousticImpedanceFunctor.h>

namespace itk
{

/** \class LabelMapToAcousticImpedanceImageFilter
 *
 * \brief Creates an image of approximate acoustic impedance from a label map of
 * classified tissues.
 *
 */
template< class TInputImage, class TOutputImage, class TLookupTable >
class LabelMapToAcousticImpedanceImageFilter:
  public UnaryFunctorImageFilter< TInputImage, TOutputImage,
    Functor::LabelMapToAcousticImpedanceFunctor< typename TInputImage::PixelType,
      typename TOutputImage::PixelType, TLookupTable > >
{
public:
  /** Standard class typedefs. */
  typedef LabelMapToAcousticImpedanceImageFilter                            Self;
  typedef UnaryFunctorImageFilter< TInputImage,
            TOutputImage,
            Functor::LabelMapToAcousticImpedanceFunctor<
              typename TInputImage::PixelType,
              typename TOutputImage::PixelType,
              TLookupTable > >                                              Superclass;
  typedef SmartPointer< Self >                                              Pointer;
  typedef SmartPointer< const Self >                                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime type information. */
  itkTypeMacro( LabelMapToAcousticImpedanceImageFilter,
                UnaryFunctorImageFilter );

  typedef typename Superclass::FunctorType FunctorType;

protected:
  LabelMapToAcousticImpedanceImageFilter() {}
  virtual ~LabelMapToAcousticImpedanceImageFilter() {}

  virtual void BeforeThreadedGenerateData();

private:
  LabelMapToAcousticImpedanceImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & ); // purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelMapToAcousticImpedanceImageFilter.txx"
#endif

#endif
