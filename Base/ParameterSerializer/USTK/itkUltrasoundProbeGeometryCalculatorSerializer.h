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

#ifndef __itkUltrasoundProbeGeometryCalculatorSerializer_h
#define __itkUltrasoundProbeGeometryCalculatorSerializer_h

#include <itkDoubleValue.h>
#include <itkObjectSerializer.h>
#include <itkUnsignedIntegerValue.h>

namespace itk
{

namespace tube
{

/** \class UltrasoundProbeGeometryCalculatorSerializer
 * \brief Parameter serializer for UltrasoundProbeGeometryCalculator.
 * \sa ParameterSerializer
 */
template< class TUltrasoundProbeGeometryCalculator >
class UltrasoundProbeGeometryCalculatorSerializer
  : public ObjectSerializer
{
public:
  /** Standard class typedefs. */
  typedef UltrasoundProbeGeometryCalculatorSerializer<
    TUltrasoundProbeGeometryCalculator > Self;
  typedef ObjectSerializer               Superclass;
  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( UltrasoundProbeGeometryCalculatorSerializer,
    ObjectSerializer );

  typedef TUltrasoundProbeGeometryCalculator UltrasoundProbeGeometryCalculatorType;

  virtual void Serialize( void );
  virtual void DeSerialize( void );

protected:
  UltrasoundProbeGeometryCalculatorSerializer( void );
  virtual ~UltrasoundProbeGeometryCalculatorSerializer( void );

  UnsignedIntegerValue * m_GeneralBeamDirection;
  DoubleValue *          m_BackgroundValue;

private:
  UltrasoundProbeGeometryCalculatorSerializer( const Self & );
  void operator=( const Self & ); // purposely not implemented

}; // End class UltrasoundProbeGeometryCalculatorSerializer

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkUltrasoundProbeGeometryCalculatorSerializer.hxx"
#endif

#endif // End !defined(__itkUltrasoundProbeGeometryCalculatorSerializer_h)
