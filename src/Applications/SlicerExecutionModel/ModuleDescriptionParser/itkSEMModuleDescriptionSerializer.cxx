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
#include "itkSEMModuleDescriptionSerializer.h"

#include "itkSEMModuleParameterGroupSerializer.h"
#include "itkParameterSerializerValue.h"
#include "itkParameterSerializerArrayValue.h"

namespace itk
{

SEMModuleDescriptionSerializer
::SEMModuleDescriptionSerializer()
{
  m_Category = new StringValue;
  this->m_Parameters["Category"] = m_Category;
  m_Contributor = new StringValue;
  this->m_Parameters["Contributor"] = m_Contributor;
  m_Description = new StringValue;
  this->m_Parameters["Description"] = m_Description;
  m_DocumentationURL = new StringValue;
  this->m_Parameters["DocumentationURL"] = m_DocumentationURL;
  m_License = new StringValue;
  this->m_Parameters["License"] = m_License;
  m_Title = new StringValue;
  this->m_Parameters["Title"] = m_Title;
  m_Version = new StringValue;
  this->m_Parameters["Version"] = m_Version;

  m_ParameterGroups = new ParameterSerializerArrayValue;
  this->m_Parameters["ParameterGroups"] = m_ParameterGroups;

  this->SetTargetObject( SEMModuleDescription::New() );
}


SEMModuleDescriptionSerializer
::~SEMModuleDescriptionSerializer()
{
  delete m_Category;
  delete m_Contributor;
  delete m_Description;
  delete m_DocumentationURL;
  delete m_License;
  delete m_Title;
  delete m_Version;
  delete m_ParameterGroups;
}


void
SEMModuleDescriptionSerializer
::Serialize()
{
  SEMModuleDescription * description =
    dynamic_cast< SEMModuleDescription * >
    ( this->GetTargetObject() );
  if( description == NULL )
    {
    itkWarningMacro("SEMModuleDescriptionSerializer target object not set");
    }
  else
    {
    const ModuleDescription & moduleDescription =
      description->GetModuleDescription();
    m_Category->SetValue( moduleDescription.GetCategory() );
    this->m_Parameters["Category"] = m_Category;
    m_Contributor->SetValue( moduleDescription.GetContributor() );
    this->m_Parameters["Contributor"] = m_Contributor;
    m_Description->SetValue( moduleDescription.GetDescription() );
    this->m_Parameters["Description"] = m_Description;
    m_DocumentationURL->SetValue( moduleDescription.GetDocumentationURL() );
    this->m_Parameters["DocumentationURL"] = m_DocumentationURL;
    m_License->SetValue( moduleDescription.GetLicense() );
    this->m_Parameters["License"] = m_License;
    m_Title->SetValue( moduleDescription.GetTitle() );
    this->m_Parameters["Title"] = m_Title;
    m_Version->SetValue( moduleDescription.GetVersion() );
    this->m_Parameters["Version"] = m_Version;

    const std::vector< ModuleParameterGroup > & parameterGroups =
      moduleDescription.GetParameterGroups();
    ParameterSerializerArrayValue::ParameterSerializerSmartArrayType
      parameterGroupSerializers( parameterGroups.size() );
    for( size_t ii = 0; ii < parameterGroups.size(); ++ii )
      {
      SEMModuleParameterGroup::Pointer semGroup =
        SEMModuleParameterGroup::New();
      semGroup->SetModuleParameterGroup( parameterGroups[ii] );

      SEMModuleParameterGroupSerializer::Pointer parameterGroupSerializer =
        SEMModuleParameterGroupSerializer::New();
      parameterGroupSerializer->SetTargetObject( semGroup.GetPointer() );

      parameterGroupSerializers[ii] = parameterGroupSerializer;
      }
    m_ParameterGroups->SetValue( parameterGroupSerializers );
    this->m_Parameters["ParameterGroups"] = m_ParameterGroups;
    }

  Superclass::Serialize();
}


void
SEMModuleDescriptionSerializer
::DeSerialize()
{
  SEMModuleDescription * description =
    dynamic_cast< SEMModuleDescription * >
    ( this->GetTargetObject() );

  if( description != NULL )
    {
    // First, we must instantiate the ParameterSerializerArray
    // with one element of the correct ParameterSerializer type
    // (ParameterGroupSerializer), so it can be replicated as needed
    ParameterSerializerArrayValue::ParameterSerializerSmartArrayType
      serializerArray( 1 );
    serializerArray[0] = SEMModuleParameterGroupSerializer::New();
    this->m_ParameterGroups->SetValue( serializerArray );
    }

  Superclass::DeSerialize();

  if( description != NULL )
    {
    ModuleDescription & moduleDescription =
      description->GetModuleDescription();
    moduleDescription.SetCategory( this->m_Category->GetValue() );
    moduleDescription.SetContributor( this->m_Contributor->GetValue() );
    moduleDescription.SetDescription( this->m_Description->GetValue() );
    moduleDescription.SetDocumentationURL( this->m_DocumentationURL->GetValue() );
    moduleDescription.SetLicense( this->m_License->GetValue() );
    moduleDescription.SetTitle( this->m_Title->GetValue() );
    moduleDescription.SetVersion( this->m_Version->GetValue() );

    const ParameterSerializerArrayValue::ParameterSerializerSmartArrayType &
      groupsSerializerArray = this->m_ParameterGroups->GetValue();
    std::vector< ModuleParameterGroup > parameterGroups( groupsSerializerArray.size() );
    for( size_t ii = 0; ii < parameterGroups.size(); ++ii )
      {
      SEMModuleParameterGroupSerializer * groupSerializer =
        static_cast< SEMModuleParameterGroupSerializer * >( groupsSerializerArray[ii].GetPointer() );
      const SEMModuleParameterGroup * semGroup =
        dynamic_cast< const SEMModuleParameterGroup * >( groupSerializer->GetTargetObject() );
      parameterGroups[ii] = semGroup->GetModuleParameterGroup();
      }
    moduleDescription.SetParameterGroups( parameterGroups );
    }
}

} // end namespace itk
