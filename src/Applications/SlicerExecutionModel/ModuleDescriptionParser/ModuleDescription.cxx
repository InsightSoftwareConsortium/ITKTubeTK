/*=========================================================================

  Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   Module Description Parser
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#include "ModuleDescription.h"
#include "ModuleDescriptionUtilities.h"

#include <sstream>
#include <fstream>
#include <string>

//----------------------------------------------------------------------------
ModuleDescription::ModuleDescription()
{
  this->Title = "Unknown";
  this->Type = "Unknown";
  this->Description = "No description provided";
  this->Category = "Unspecified";
  this->Version = "Unspecified";
  this->DocumentationURL = "";
  this->License = "";
  this->Acknowledgements = "Thank you everyone.";
  this->Contributor = "Anonymous";
  this->Target = "";
  this->Location = "";

  std::stringstream ss;
  ss << (unsigned short) -1;
  ss >> this->Index;

  this->LibraryLoader = 0;
  this->TargetCallback = 0;
}

//----------------------------------------------------------------------------
ModuleDescription::ModuleDescription(const ModuleDescription &md)
{
  this->Title = md.Title;
  this->Category = md.Category;
  this->Index = md.Index;
  this->Description = md.Description;
  this->Version = md.Version;
  this->DocumentationURL = md.DocumentationURL;
  this->License = md.License;
  this->Acknowledgements = md.Acknowledgements;
  this->Contributor = md.Contributor;
  this->Type = md.Type;
  this->Target = md.Target;
  this->Location = md.Location;
  this->ParameterGroups = md.ParameterGroups;
  this->Logo = md.Logo;
  this->LibraryLoader = md.LibraryLoader;
  this->TargetCallback = md.TargetCallback;

  this->ProcessInformation.Initialize();
}

//----------------------------------------------------------------------------
void ModuleDescription::operator=(const ModuleDescription &md)
{
  this->Title = md.Title;
  this->Category = md.Category;
  this->Index = md.Index;
  this->Description = md.Description;
  this->Version = md.Version;
  this->DocumentationURL = md.DocumentationURL;
  this->License = md.License;
  this->Acknowledgements = md.Acknowledgements;
  this->Contributor = md.Contributor;
  this->Type= md.Type;
  this->Target = md.Target;
  this->Location = md.Location;
  this->ParameterGroups = md.ParameterGroups;
  this->ProcessInformation = md.ProcessInformation;
  this->Logo = md.Logo;
  this->LibraryLoader = md.LibraryLoader;
  this->TargetCallback = md.TargetCallback;
}

//----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream &os, const ModuleDescription &module)
{
  os << "Title: " << module.GetTitle() << std::endl;
  os << "Category: " << module.GetCategory() << std::endl;
  os << "Index: " << module.GetIndex() << std::endl;
  os << "Description: " << module.GetDescription() << std::endl;
  os << "Version: " << module.GetVersion() << std::endl;
  os << "DocumentationURL: " << module.GetDocumentationURL() << std::endl;
  os << "License: " << module.GetLicense() << std::endl;
  os << "Contributor: " << module.GetContributor() << std::endl;
  os << "Acknowledgements: " << module.GetAcknowledgements() << std::endl;
  os << "Type: " << module.GetType() << std::endl;
  os << "Target: " << module.GetTarget() << std::endl;
  os << "Location: " << module.GetLocation() << std::endl;
  //os << "Logo: " << module.GetLogo() << std::endl;

  os << "ProcessInformation: " << std::endl
     << *(module.GetProcessInformation());

  os << "ParameterGroups: " << std::endl;
  std::vector<ModuleParameterGroup>::const_iterator it = module.GetParameterGroups().begin();
  while (it != module.GetParameterGroups().end())
    {
    os << *it;
    ++it;
    }
  return os;
}


//----------------------------------------------------------------------------
bool ModuleDescription::HasParameter(const std::string& name) const
{
  // iterate over each parameter group
  std::vector<ModuleParameterGroup>::const_iterator pgbeginit
    = this->ParameterGroups.begin();
  std::vector<ModuleParameterGroup>::const_iterator pgendit
    = this->ParameterGroups.end();
  std::vector<ModuleParameterGroup>::const_iterator pgit;

  for (pgit = pgbeginit; pgit != pgendit; ++pgit)
    {
    // iterate over each parameter in this group
    std::vector<ModuleParameter>::const_iterator pbeginit
      = (*pgit).GetParameters().begin();
    std::vector<ModuleParameter>::const_iterator pendit
      = (*pgit).GetParameters().end();
    std::vector<ModuleParameter>::const_iterator pit;

    for (pit = pbeginit; pit != pendit; ++pit)
      {
      if ((*pit).GetName() == name)
        {
        return true;
        }
      }
    }

  return false;
}

//----------------------------------------------------------------------------
bool ModuleDescription::HasReturnParameters() const
{
  // iterate over each parameter group
  std::vector<ModuleParameterGroup>::const_iterator pgbeginit
    = this->ParameterGroups.begin();
  std::vector<ModuleParameterGroup>::const_iterator pgendit
    = this->ParameterGroups.end();
  std::vector<ModuleParameterGroup>::const_iterator pgit;

  for (pgit = pgbeginit; pgit != pgendit; ++pgit)
    {
    // iterate over each parameter in this group
    std::vector<ModuleParameter>::const_iterator pbeginit
      = (*pgit).GetParameters().begin();
    std::vector<ModuleParameter>::const_iterator pendit
      = (*pgit).GetParameters().end();
    std::vector<ModuleParameter>::const_iterator pit;

    for (pit = pbeginit; pit != pendit; ++pit)
      {
      if ((*pit).IsReturnParameter())
        {
        return true;
        }
      }
    }

  return false;
}


//----------------------------------------------------------------------------
std::vector<ModuleParameter> ModuleDescription
::FindParametersWithDefaultValue(const std::string& value) const
{
  return this->FindParametersWithValue(value);
}


//----------------------------------------------------------------------------
std::vector<ModuleParameter> ModuleDescription
::FindParametersWithValue(const std::string& value) const
{
  std::vector<ModuleParameter> parameters;
  // iterate over each parameter group
  std::vector<ModuleParameterGroup>::const_iterator pgbeginit
    = this->ParameterGroups.begin();
  std::vector<ModuleParameterGroup>::const_iterator pgendit
    = this->ParameterGroups.end();
  std::vector<ModuleParameterGroup>::const_iterator pgit;

  for (pgit = pgbeginit; pgit != pgendit; ++pgit)
    {
    // iterate over each parameter in this group
    std::vector<ModuleParameter>::const_iterator pbeginit
      = (*pgit).GetParameters().begin();
    std::vector<ModuleParameter>::const_iterator pendit
      = (*pgit).GetParameters().end();
    std::vector<ModuleParameter>::const_iterator pit;

    for (pit = pbeginit; pit != pendit; ++pit)
      {
      if ((*pit).GetValue() == value)
        {
        parameters.push_back(*pit);
        }
      }
    }

  return parameters;
}

//----------------------------------------------------------------------------
bool ModuleDescription::SetParameterDefaultValue(const std::string& name, const std::string& value)
{
  return this->SetParameterValue( name, value);
}

//----------------------------------------------------------------------------
bool ModuleDescription::SetParameterValue(const std::string& name, const std::string& value)
{
  // iterate over each parameter group
  std::vector<ModuleParameterGroup>::iterator pgbeginit
    = this->ParameterGroups.begin();
  std::vector<ModuleParameterGroup>::iterator pgendit
    = this->ParameterGroups.end();
  std::vector<ModuleParameterGroup>::iterator pgit;

  for (pgit = pgbeginit; pgit != pgendit; ++pgit)
    {
    // iterate over each parameter in this group
    std::vector<ModuleParameter>::iterator pbeginit
      = (*pgit).GetParameters().begin();
    std::vector<ModuleParameter>::iterator pendit
      = (*pgit).GetParameters().end();
    std::vector<ModuleParameter>::iterator pit;

    for (pit = pbeginit; pit != pendit; ++pit)
      {
      if ((*pit).GetName() == name)
        {
        (*pit).SetValue(value);
        return true;
        }
      }
    }

  return false;
}

//----------------------------------------------------------------------------
std::string ModuleDescription::GetParameterDefaultValue(const std::string& name) const
{
  return this->GetParameterValue(name);
}

//----------------------------------------------------------------------------
std::string ModuleDescription::GetParameterValue(const std::string& name) const
{
  // iterate over each parameter group
  std::vector<ModuleParameterGroup>::const_iterator pgbeginit
    = this->ParameterGroups.begin();
  std::vector<ModuleParameterGroup>::const_iterator pgendit
    = this->ParameterGroups.end();
  std::vector<ModuleParameterGroup>::const_iterator pgit;

  for (pgit = pgbeginit; pgit != pgendit; ++pgit)
    {
    // iterate over each parameter in this group
    std::vector<ModuleParameter>::const_iterator pbeginit
      = (*pgit).GetParameters().begin();
    std::vector<ModuleParameter>::const_iterator pendit
      = (*pgit).GetParameters().end();
    std::vector<ModuleParameter>::const_iterator pit;

    for (pit = pbeginit; pit != pendit; ++pit)
      {
      if ((*pit).GetName() == name)
        {
        return (*pit).GetValue();
        }
      }
    }

  return "";
}

//----------------------------------------------------------------------------
void ModuleDescription ::SetLogo(const ModuleLogo& logo)
{
  this->Logo = logo;
}

//----------------------------------------------------------------------------
const ModuleLogo& ModuleDescription::GetLogo() const
{
  return this->Logo;
}

//----------------------------------------------------------------------------
bool ModuleDescription ::ReadParameterFile(const std::string& filename)
{
  std::ifstream rtp;
  bool modified = false;

  rtp.open(filename.c_str());
  if (rtp.fail())
    {
    std::cout << "Parameter file " << filename << " could not be opened." << std::endl;
    return false;
    }

  std::string line;
  while (std::getline(rtp, line))
    {
    // split the line into key: value
    std::string key, value;

    std::string::size_type start = line.find_first_not_of(" \t");
    std::string::size_type stop = line.find_first_of("=", start);

    key = line.substr(start, stop-start);
    start = line.find_first_not_of(" \t", stop+1);
    if (start != std::string::npos)
      {
      value = line.substr(start, line.length() - start + 1);
      }

    trimLeadingAndTrailing(key);
    trimLeadingAndTrailing(value);

    // std::cout << "key=" << key << ", value=" << value << "!" << std::endl;

    if (this->HasParameter(key))
      {
      if (value != this->GetParameterValue(key))
        {
        this->SetParameterValue(key, value);
        modified = true;

        // multiple="true" may have to be handled differently
        }
      }
    }

  rtp.close();
  return modified;
}

//----------------------------------------------------------------------------
bool ModuleDescription::
WriteParameterFile(const std::string& filename, bool withHandlesToBulkParameters)
{
  std::ofstream rtp;

  rtp.open(filename.c_str());
  if (rtp.fail())
    {
    std::cout << "Parameter file " << filename << " could not be opened for writing." << std::endl;
    return false;
    }

  // iterate over each parameter group
  std::vector<ModuleParameterGroup>::const_iterator pgbeginit
    = this->ParameterGroups.begin();
  std::vector<ModuleParameterGroup>::const_iterator pgendit
    = this->ParameterGroups.end();
  std::vector<ModuleParameterGroup>::const_iterator pgit;

  for (pgit = pgbeginit; pgit != pgendit; ++pgit)
    {
    // iterate over each parameter in this group
    std::vector<ModuleParameter>::const_iterator pbeginit
      = (*pgit).GetParameters().begin();
    std::vector<ModuleParameter>::const_iterator pendit
      = (*pgit).GetParameters().end();
    std::vector<ModuleParameter>::const_iterator pit;

    for (pit = pbeginit; pit != pendit; ++pit)
      {
      // write out all parameters or just the ones that are not bulk parameters
      if (withHandlesToBulkParameters
          || (!withHandlesToBulkParameters
              && ((*pit).GetTag() != "image"
                  && (*pit).GetTag() != "geometry"
                  && (*pit).GetTag() != "transform"
                  && (*pit).GetTag() != "table"
                  && (*pit).GetTag() != "measurement"
                  && (*pit).GetTag() != "point"  // point and point file and region are special
                  && (*pit).GetTag() != "pointfile"
                  && (*pit).GetTag() != "region")))
        {
        rtp << (*pit).GetName() << " = "
            << (*pit).GetValue() << std::endl;

        // multiple="true" may have to be handled differently
        }
      }
    }

  rtp.close();
  return true;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetCategory(const std::string &cat)
{
  this->Category = cat;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetCategory() const
{
  return this->Category;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetIndex(const std::string &ind)
{
  this->Index = ind;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetIndex() const
{
  return this->Index;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetTitle(const std::string &title)
{
  this->Title = title;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetTitle() const
{
  return this->Title;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetDescription(const std::string &description)
{
  this->Description = description;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetDescription() const
{
  return this->Description;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetVersion(const std::string &version)
{
  this->Version = version;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetVersion() const
{
  return this->Version;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetDocumentationURL(const std::string &documentationURL)
{
  this->DocumentationURL = documentationURL;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetDocumentationURL() const
{
  return this->DocumentationURL;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetLicense(const std::string &license)
{
  this->License = license;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetLicense() const
{
  return this->License;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetAcknowledgements(const std::string &acknowledgements)
{
  this->Acknowledgements = acknowledgements;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetAcknowledgements() const
{
  return this->Acknowledgements;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetContributor(const std::string &contributor)
{
  this->Contributor = contributor;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetContributor() const
{
  return this->Contributor;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetType(const std::string &type)
{
  if (type == "SharedObjectModule"
      || type == "CommandLineModule")
    {
    this->Type = type;
    }
  else
    {
    this->Type = "Unknown";
    }
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetType() const
{
  return this->Type;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetTarget(const std::string &target)
{
  this->Target = target;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetTarget() const
{
  // Lazyily set the target
  if (this->Target.empty() && this->LibraryLoader && this->TargetCallback)
    {
    (*TargetCallback)(this->LibraryLoader, *const_cast<ModuleDescription*>(this));
    }
  return this->Target;
}

//----------------------------------------------------------------------------
void ModuleDescription::SetLocation(const std::string &target)
{
  this->Location = target;
}

//----------------------------------------------------------------------------
const std::string& ModuleDescription::GetLocation() const
{
  return this->Location;
}

//----------------------------------------------------------------------------
void ModuleDescription::AddParameterGroup(const ModuleParameterGroup &group)
{
  this->ParameterGroups.push_back(group);
}

//----------------------------------------------------------------------------
const std::vector<ModuleParameterGroup>&
ModuleDescription::GetParameterGroups() const
{
  return this->ParameterGroups;
}

//----------------------------------------------------------------------------
std::vector<ModuleParameterGroup>& ModuleDescription::GetParameterGroups()
{
  return this->ParameterGroups;
}

//----------------------------------------------------------------------------
void ModuleDescription::
SetParameterGroups(const std::vector<ModuleParameterGroup>& groups)
{
    this->ParameterGroups = groups;
}

//----------------------------------------------------------------------------
const ModuleProcessInformation* ModuleDescription::GetProcessInformation() const
{
  return &ProcessInformation;
}

//----------------------------------------------------------------------------
ModuleProcessInformation* ModuleDescription::GetProcessInformation()
{
  return &ProcessInformation;
}

//----------------------------------------------------------------------------
void ModuleDescription::
SetTargetCallback(void* libraryLoader, TargetCallbackType targetCallback)
{
  this->LibraryLoader = libraryLoader;
  this->TargetCallback = targetCallback;
}

//----------------------------------------------------------------------------
void* ModuleDescription::GetLibraryLoader()const
{
  return this->LibraryLoader;
}

//----------------------------------------------------------------------------
ModuleDescription::TargetCallbackType ModuleDescription::GetTargetCallback()
{
  return this->TargetCallback;
}
