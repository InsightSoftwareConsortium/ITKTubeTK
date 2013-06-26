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

#ifndef __itkAcousticImpulseResponseImageFilterSerializer_h
#define __itkAcousticImpulseResponseImageFilterSerializer_h

#include <itkImageToImageFilterSerializer.h>

namespace itk
{

namespace tube
{

/** \class AcousticImpulseResponseImageFilterSerializer
 * \brief Parameter serializer for AcousticImpulseResponseImageFilter.
 * \sa ParameterSerializer
 */
template< class TAcousticImpulseResponseImageFilter >
class AcousticImpulseResponseImageFilterSerializer :
  public ImageToImageFilterSerializer< typename TAcousticImpulseResponseImageFilter::Superclass >
{
public:
  /** Standard class typedefs. */
  typedef AcousticImpulseResponseImageFilterSerializer  Self;
  typedef ImageToImageFilterSerializer
    < typename TAcousticImpulseResponseImageFilter::Superclass >
                                                        Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AcousticImpulseResponseImageFilterSerializer,
    ImageToImageFilterSerializer );

  typedef TAcousticImpulseResponseImageFilter AcousticImpulseResponseImageFilterType;

  virtual void Serialize( void );
  virtual void DeSerialize( void );

protected:
  AcousticImpulseResponseImageFilterSerializer( void );
  virtual ~AcousticImpulseResponseImageFilterSerializer( void );

  typedef ImageToImageFilterSerializer
    < typename AcousticImpulseResponseImageFilterType::GradientMagnitudeFilterType >
      GradientMagnitudeFilterSerializerType;

  DoubleValue *                 m_AngleDependence;

  ParameterSerializerValue *    m_GradientMagnitudeFilter;
  typename GradientMagnitudeFilterSerializerType::Pointer
    m_GradientMagnitudeFilterSerializer;

  /** Set the TargetObject on m_GradientMagnitudeFilter after
   * making sure m_GradientMagnitudeFilterSerializer is the correct
   * type. */
  virtual void AssignGradientMagnitudeFilter( AcousticImpulseResponseImageFilterType * filter );

private:
  AcousticImpulseResponseImageFilterSerializer( const Self & );
  void operator=( const Self & ); // purposely not implemented

}; // End class AcousticImpulseResponseImageFilterSerializer

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAcousticImpulseResponseImageFilterSerializer.hxx"
#endif

#endif // End !defined(__itkAcousticImpulseResponseImageFilterSerializer_h)
