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
#include "itkSEMModuleParameterGroupSerializer.h"

#include "itkSEMModuleParameterSerializer.h"
#include "itkParameterSerializerValue.h"
#include "itkParameterSerializerArrayValue.h"

namespace itk
{

SEMModuleParameterGroupSerializer
::SEMModuleParameterGroupSerializer()
{
  m_Description = new StringValue;
  this->m_Parameters["Description"] = m_Description;
  m_Label = new StringValue;
  this->m_Parameters["Label"] = m_Label;

  m_ModuleParameters = new ParameterSerializerArrayValue;
  this->m_Parameters["Parameters"] = m_ModuleParameters;

  this->SetTargetObject( SEMModuleParameterGroup::New() );
}


SEMModuleParameterGroupSerializer
::~SEMModuleParameterGroupSerializer()
{
  delete m_Description;
  delete m_Label;
  delete m_ModuleParameters;
}


void
SEMModuleParameterGroupSerializer
::Serialize()
{
  SEMModuleParameterGroup * parameterGroup =
    dynamic_cast< SEMModuleParameterGroup * >
      ( this->GetTargetObject() );
  if( parameterGroup == NULL )
    {
    itkWarningMacro("SEMModuleParameterGroupSerializer target object not set");
    }
  else
    {
    const ModuleParameterGroup & groups =
      parameterGroup->GetModuleParameterGroup();
    m_Description->SetValue( groups.GetDescription() );
    this->m_Parameters["Description"] = m_Description;
    m_Label->SetValue( groups.GetLabel() );
    this->m_Parameters["Label"] = m_Label;

    const std::vector< ModuleParameter > & parameters =
      groups.GetParameters();
    ParameterSerializerArrayValue::ParameterSerializerSmartArrayType
      parameterSerializers( parameters.size() );
    for( size_t ii = 0; ii < parameters.size(); ++ii )
      {
      SEMModuleParameter::Pointer parameter =
        SEMModuleParameter::New();
      parameter->SetModuleParameter( parameters[ii] );

      SEMModuleParameterSerializer::Pointer parameterSerializer =
        SEMModuleParameterSerializer::New();
      parameterSerializer->SetTargetObject( parameter.GetPointer() );

      parameterSerializers[ii] = parameterSerializer;
      }
    m_ModuleParameters->SetValue( parameterSerializers );
    this->m_Parameters["Parameters"] = m_ModuleParameters;
    }

  Superclass::Serialize();
}


void
SEMModuleParameterGroupSerializer
::DeSerialize()
{
  SEMModuleParameterGroup * parameterGroup =
    dynamic_cast< SEMModuleParameterGroup * >
      ( this->GetTargetObject() );

  if( parameterGroup != NULL )
    {
    // First, we must instantiate the ParameterSerializerArray
    // with one element of the correct ParameterSerializer type
    // (ParameterSerializer), so it can be replicated as needed
    ParameterSerializerArrayValue::ParameterSerializerSmartArrayType
      serializerArray( 1 );
    serializerArray[0] = SEMModuleParameterSerializer::New();
    this->m_ModuleParameters->SetValue( serializerArray );
    }

  Superclass::DeSerialize();

  if( parameterGroup != NULL )
    {
    ModuleParameterGroup & group = parameterGroup->GetModuleParameterGroup();
    group.SetDescription( this->m_Description->GetValue() );
    group.SetLabel( this->m_Label->GetValue() );

    const ParameterSerializerArrayValue::ParameterSerializerSmartArrayType &
      parameterSerializerArray = this->m_ModuleParameters->GetValue();
    std::vector< ModuleParameter > parameters( parameterSerializerArray.size() );
    for( size_t ii = 0; ii < parameters.size(); ++ii )
      {
      SEMModuleParameterSerializer * serializer =
        static_cast< SEMModuleParameterSerializer * >( parameterSerializerArray[ii].GetPointer() );
      const SEMModuleParameter * semParameter =
        dynamic_cast< const SEMModuleParameter * >( serializer->GetTargetObject() );
      parameters[ii] = semParameter->GetModuleParameter();
      }
    group.SetParameters( parameters );
    }
}

} // end namespace itk
