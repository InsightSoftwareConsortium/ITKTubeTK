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
#ifndef __itkSEMModuleDescriptionSerializer_h
#define __itkSEMModuleDescriptionSerializer_h

#include "ModuleDescription.h"

#include "itkParameterSerializer.h"

#include "itkObjectFactory.h"

namespace itk
{

/** \class ModuleDescription
 *
 * \brief An itk::LightObject proxy for the SEM ModuleDescription class.
 *
 */
class ModuleDescriptionParser_EXPORT SEMModuleDescription: public LightObject
{
public:
  typedef SEMModuleDescription       Self;
  typedef LightObject                Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );

  void SetModuleDescription( const ModuleDescription & module )
    {
    this->m_ModuleDescription = module;
    }

  const ModuleDescription & GetModuleDescription() const
    {
    return this->m_ModuleDescription;
    }

  ModuleDescription & GetModuleDescription()
    {
    return const_cast< ModuleDescription & >(
      static_cast< const SEMModuleDescription & >(*this)
        .GetModuleDescription() );
    }

protected:
  SEMModuleDescription() {}
  ~SEMModuleDescription() {}

  ModuleDescription m_ModuleDescription;

private:
  SEMModuleDescription( const SEMModuleDescription & ); // purposely not implemented
  void operator=( const SEMModuleDescription & ); // purposely not implemented
};

/** \class SEMModuleDescriptionSerializer
 *
 * \brief Parameter serializer for ITKModuleDescription.
 *
 * \sa ParameterSerializer
 *
 */
class ModuleDescriptionParser_EXPORT SEMModuleDescriptionSerializer:
  public ParameterSerializer
{
public:
  /** Standard class typedefs. */
  typedef SEMModuleDescriptionSerializer  Self;
  typedef ParameterSerializer             Superclass;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SEMModuleDescriptionSerializer,
    ParameterSerializer );

  virtual void Serialize();
  virtual void DeSerialize();

protected:
  SEMModuleDescriptionSerializer();
  virtual ~SEMModuleDescriptionSerializer();

  StringValue * m_Category;
  StringValue * m_Contributor;
  StringValue * m_Description;
  StringValue * m_DocumentationURL;
  StringValue * m_License;
  StringValue * m_Title;
  StringValue * m_Version;

  ParameterSerializerArrayValue * m_ParameterGroups;

private:
  SEMModuleDescriptionSerializer( const Self & );
  void operator=( const Self & ); // purposely not implemented
};

} // end namespace itk

#endif
