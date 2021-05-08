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
#include "itkSEMModuleParameterSerializer.h"

#include "itkParameterSerializerValue.h"

namespace itk
{

SEMModuleParameterSerializer
::SEMModuleParameterSerializer()
{
  m_Channel = new StringValue;
  this->m_Parameters["Channel"] = m_Channel;
  m_CoordinateSystem = new StringValue;
  this->m_Parameters["CoordinateSystem"] = m_CoordinateSystem;
  m_CPPType = new StringValue;
  this->m_Parameters["CPPType"] = m_CPPType;
  m_Value = new StringValue;
  this->m_Parameters["Value"] = m_Value;
  m_Description = new StringValue;
  this->m_Parameters["Description"] = m_Description;
  m_Flag = new StringValue;
  this->m_Parameters["Flag"] = m_Flag;
  m_Index = new StringValue;
  this->m_Parameters["Index"] = m_Index;
  m_Label = new StringValue;
  this->m_Parameters["Label"] = m_Label;
  m_LongFlag = new StringValue;
  this->m_Parameters["LongFlag"] = m_LongFlag;
  m_Maximum = new StringValue;
  this->m_Parameters["Maximum"] = m_Maximum;
  m_Minimum = new StringValue;
  this->m_Parameters["Minimum"] = m_Minimum;
  m_Multiple = new StringValue;
  this->m_Parameters["Multiple"] = m_Multiple;
  m_Name = new StringValue;
  this->m_Parameters["Name"] = m_Name;
  m_Step = new StringValue;
  this->m_Parameters["Step"] = m_Step;
  m_Tag = new StringValue;
  this->m_Parameters["Tag"] = m_Tag;
  m_Elements = new StringArrayValue;
  this->m_Parameters["Elements"] = m_Elements;

  this->SetTargetObject( SEMModuleParameter::New() );
}


SEMModuleParameterSerializer
::~SEMModuleParameterSerializer()
{
  delete m_Channel;
  delete m_CPPType;
  delete m_CoordinateSystem;
  delete m_Value;
  delete m_Description;
  delete m_Flag;
  delete m_Index;
  delete m_Label;
  delete m_LongFlag;
  delete m_Maximum;
  delete m_Minimum;
  delete m_Multiple;
  delete m_Name;
  delete m_Step;
  delete m_Tag;
  delete m_Elements;
}


void
SEMModuleParameterSerializer
::Serialize()
{
  SEMModuleParameter * moduleParameter =
    dynamic_cast< SEMModuleParameter * >
      ( this->GetTargetObject() );
  if( moduleParameter == NULL )
    {
    itkWarningMacro("SEMModuleParameterSerializer target object not set");
    }
  else
    {
    const ModuleParameter & parameter =
      moduleParameter->GetModuleParameter();
    m_Channel->SetValue( parameter.GetChannel() );
    this->m_Parameters["Channel"] = m_Channel;
    m_CoordinateSystem->SetValue( parameter.GetCoordinateSystem() );
    this->m_Parameters["CoordinateSystem"] = m_CoordinateSystem;
    m_CPPType->SetValue( parameter.GetCPPType() );
    this->m_Parameters["CPPType"] = m_CPPType;
    m_Value->SetValue( parameter.GetValue() );
    this->m_Parameters["Value"] = m_Value;
    m_Description->SetValue( parameter.GetDescription() );
    this->m_Parameters["Description"] = m_Description;
    m_Flag->SetValue( parameter.GetFlag() );
    this->m_Parameters["Flag"] = m_Flag;
    m_Index->SetValue( parameter.GetIndex() );
    this->m_Parameters["Index"] = m_Index;
    m_Label->SetValue( parameter.GetLabel() );
    this->m_Parameters["Label"] = m_Label;
    m_LongFlag->SetValue( parameter.GetLongFlag() );
    this->m_Parameters["LongFlag"] = m_LongFlag;
    m_Maximum->SetValue( parameter.GetMaximum() );
    this->m_Parameters["Maximum"] = m_Maximum;
    m_Minimum->SetValue( parameter.GetMinimum() );
    this->m_Parameters["Minimum"] = m_Minimum;
    m_Multiple->SetValue( parameter.GetMultiple() );
    this->m_Parameters["Multiple"] = m_Multiple;
    m_Name->SetValue( parameter.GetName() );
    this->m_Parameters["Name"] = m_Name;
    m_Step->SetValue( parameter.GetStep() );
    this->m_Parameters["Step"] = m_Step;
    m_Tag->SetValue( parameter.GetTag() );
    this->m_Parameters["Tag"] = m_Tag;
    m_Elements->SetValue( parameter.GetElements() );
    this->m_Parameters["Elements"] = m_Elements;
    }

  Superclass::Serialize();
}


void
SEMModuleParameterSerializer
::DeSerialize()
{
  SEMModuleParameter * moduleParameter =
    dynamic_cast< SEMModuleParameter * >
      ( this->GetTargetObject() );

  Superclass::DeSerialize();

  if( moduleParameter != NULL )
    {
    ModuleParameter & parameter =
      moduleParameter->GetModuleParameter();
    parameter.SetChannel( this->m_Channel->GetValue() );
    parameter.SetCPPType( this->m_CPPType->GetValue() );
    parameter.SetCoordinateSystem( this->m_CoordinateSystem->GetValue() );
    parameter.SetValue( this->m_Value->GetValue() );
    parameter.SetDescription( this->m_Description->GetValue() );
    parameter.SetFlag( this->m_Flag->GetValue() );
    parameter.SetIndex( this->m_Index->GetValue() );
    parameter.SetLabel( this->m_Label->GetValue() );
    parameter.SetLongFlag( this->m_LongFlag->GetValue() );
    parameter.SetMaximum( this->m_Maximum->GetValue() );
    parameter.SetMinimum( this->m_Minimum->GetValue() );
    parameter.SetMultiple( this->m_Multiple->GetValue() );
    parameter.SetName( this->m_Name->GetValue() );
    parameter.SetStep( this->m_Step->GetValue() );
    parameter.SetTag( this->m_Tag->GetValue() );
    parameter.SetElements( this->m_Elements->GetValue() );
    }
}

} // end namespace itk
