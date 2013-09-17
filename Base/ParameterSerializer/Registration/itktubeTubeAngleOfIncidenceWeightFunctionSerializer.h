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

#ifndef __itktubeTubeAngleOfIncidenceWeightFunctionSerializer_h
#define __itktubeTubeAngleOfIncidenceWeightFunctionSerializer_h

#include <itkObjectSerializer.h>

#include <itkDoubleValue.h>
#include <itkDoubleArrayValue.h>

namespace itk
{

namespace tube
{

/** \class TubeAngleOfIncidenceWeightFunctionSerializer
 * \brief Parameter serializer for TubeAngleOfIncidenceWeightFunction.
 * \sa ParameterSerializer
 */
template< class TTubeAngleOfIncidenceWeightFunction >
class TubeAngleOfIncidenceWeightFunctionSerializer
  : public ObjectSerializer
{
public:
  /** Standard class typedefs. */
  typedef TubeAngleOfIncidenceWeightFunctionSerializer Self;
  typedef ObjectSerializer                             Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TubeAngleOfIncidenceWeightFunctionSerializer,
    ObjectSerializer );

  typedef TTubeAngleOfIncidenceWeightFunction TubeAngleOfIncidenceWeightFunctionType;

  virtual void Serialize( void );
  virtual void DeSerialize( void );

protected:
  TubeAngleOfIncidenceWeightFunctionSerializer( void );
  virtual ~TubeAngleOfIncidenceWeightFunctionSerializer( void );

  DoubleValue *      m_FractionalImportance;
  DoubleValue *      m_AngleDependence;
  DoubleArrayValue * m_UltrasoundProbeOrigin;

private:
  TubeAngleOfIncidenceWeightFunctionSerializer( const Self & );
  void operator=( const Self & ); // purposely not implemented

}; // End class TubeAngleOfIncidenceWeightFunctionSerializer

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeAngleOfIncidenceWeightFunctionSerializer.hxx"
#endif

#endif // End !defined(__itktubeTubeAngleOfIncidenceWeightFunctionSerializer_h)
