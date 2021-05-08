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
#ifndef __itkSEMModuleParameterSerializer_h
#define __itkSEMModuleParameterSerializer_h

#include "ModuleParameter.h"

#include "itkParameterSerializer.h"

#include "itkObjectFactory.h"

namespace itk
{

/** \class ModuleParameter
 *
 * \brief An itk::LightObject proxy for the SEM ModuleParameter class.
 *
 */
class ModuleDescriptionParser_EXPORT SEMModuleParameter: public LightObject
{
public:
  typedef SEMModuleParameter         Self;
  typedef LightObject                Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );

  void SetModuleParameter( const ModuleParameter & parameter )
    {
    this->m_ModuleParameter = parameter;
    }

  const ModuleParameter & GetModuleParameter() const
    {
    return this->m_ModuleParameter;
    }

  ModuleParameter & GetModuleParameter()
    {
    return const_cast< ModuleParameter & >(
      static_cast< const SEMModuleParameter & >(*this)
        .GetModuleParameter() );
    }

protected:
  SEMModuleParameter() {}
  ~SEMModuleParameter() {}

  ModuleParameter m_ModuleParameter;

private:
  SEMModuleParameter( const SEMModuleParameter & ); // purposely not implemented
  void operator=( const SEMModuleParameter & ); // purposely not implemented
};

/** \class SEMModuleParameterSerializer
 *
 * \brief Parameter serializer for ITKModuleParameter.
 *
 * \sa ParameterSerializer
 *
 */
class ModuleDescriptionParser_EXPORT SEMModuleParameterSerializer:
  public ParameterSerializer
{
public:
  /** Standard class typedefs. */
  typedef SEMModuleParameterSerializer  Self;
  typedef ParameterSerializer           Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SEMModuleParameterSerializer,
    ParameterSerializer );

  virtual void Serialize();
  virtual void DeSerialize();

protected:
  SEMModuleParameterSerializer();
  virtual ~SEMModuleParameterSerializer();

  StringValue *      m_Channel;
  StringValue *      m_CoordinateSystem;
  StringValue *      m_CPPType;
  StringValue *      m_Value;
  StringValue *      m_Description;
  StringValue *      m_Flag;
  StringValue *      m_Index;
  StringValue *      m_Label;
  StringValue *      m_LongFlag;
  StringValue *      m_Maximum;
  StringValue *      m_Minimum;
  StringValue *      m_Multiple;
  StringValue *      m_Name;
  StringValue *      m_Step;
  StringValue *      m_Tag;
  StringArrayValue * m_Elements;

private:
  SEMModuleParameterSerializer( const Self & );
  void operator=( const Self & ); // purposely not implemented
};

} // end namespace itk

#endif
