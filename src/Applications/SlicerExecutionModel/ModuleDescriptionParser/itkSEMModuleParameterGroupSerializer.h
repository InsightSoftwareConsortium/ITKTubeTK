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
#ifndef __itkSEMModuleParameterGroupSerializer_h
#define __itkSEMModuleParameterGroupSerializer_h

#include "ModuleParameterGroup.h"

#include "itkParameterSerializer.h"

#include "itkObjectFactory.h"

namespace itk
{

/** \class ModuleParameterGroup
 *
 * \brief An itk::LightObject proxy for the SEM ModuleParameterGroup class.
 *
 */
class ModuleDescriptionParser_EXPORT SEMModuleParameterGroup: public LightObject
{
public:
  typedef SEMModuleParameterGroup    Self;
  typedef LightObject                Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );

  void SetModuleParameterGroup( const ModuleParameterGroup & group )
    {
    this->m_ModuleParameterGroup = group;
    }

  const ModuleParameterGroup & GetModuleParameterGroup() const
    {
    return this->m_ModuleParameterGroup;
    }

  ModuleParameterGroup & GetModuleParameterGroup()
    {
    return const_cast< ModuleParameterGroup & >(
      static_cast< const SEMModuleParameterGroup & >(*this)
        .GetModuleParameterGroup() );
    }

protected:
  SEMModuleParameterGroup() {}
  ~SEMModuleParameterGroup() {}

  ModuleParameterGroup m_ModuleParameterGroup;

private:
  SEMModuleParameterGroup( const SEMModuleParameterGroup & ); // purposely not implemented
  void operator=( const SEMModuleParameterGroup & ); // purposely not implemented
};

/** \class SEMModuleParameterGroupSerializer
 *
 * \brief Parameter serializer for ITKModuleParameterGroup.
 *
 * \sa ParameterSerializer
 *
 */
class ModuleDescriptionParser_EXPORT SEMModuleParameterGroupSerializer:
  public ParameterSerializer
{
public:
  /** Standard class typedefs. */
  typedef SEMModuleParameterGroupSerializer  Self;
  typedef ParameterSerializer                Superclass;
  typedef SmartPointer< Self >               Pointer;
  typedef SmartPointer< const Self >         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SEMModuleParameterGroupSerializer,
    ParameterSerializer );

  virtual void Serialize();
  virtual void DeSerialize();

protected:
  SEMModuleParameterGroupSerializer();
  virtual ~SEMModuleParameterGroupSerializer();

  StringValue * m_Description;
  StringValue * m_Label;

  ParameterSerializerArrayValue * m_ModuleParameters;

private:
  SEMModuleParameterGroupSerializer( const Self & );
  void operator=( const Self & ); // purposely not implemented
};

} // end namespace itk

#endif
