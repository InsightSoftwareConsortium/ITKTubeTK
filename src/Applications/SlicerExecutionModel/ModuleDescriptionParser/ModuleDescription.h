/*=========================================================================

  Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   Module Description Parser
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __ModuleDescription_h
#define __ModuleDescription_h

#include "ModuleDescriptionParserExport.h"

#include "ModuleParameterGroup.h"

#include "ModuleProcessInformation.h"
#include "ModuleLogo.h"

#include <string>
#include <vector>

class ModuleDescriptionParser_EXPORT ModuleDescription
{
public:
  ModuleDescription();
  ModuleDescription(const ModuleDescription &md);

  void operator=(const ModuleDescription &md);

  void SetCategory(const std::string &cat);
  const std::string&  GetCategory() const;

  void SetIndex(const std::string &ind);
  const std::string& GetIndex() const;

  void SetTitle(const std::string &title);
  const std::string& GetTitle() const;

  void SetDescription(const std::string &description);
  const std::string& GetDescription() const;

  void SetVersion(const std::string &version);
  const std::string& GetVersion() const;

  void SetDocumentationURL(const std::string &documentationURL);
  const std::string& GetDocumentationURL() const;

  void SetLicense(const std::string &license);
  const std::string& GetLicense() const;

  void SetAcknowledgements(const std::string &acknowledgements);
  const std::string& GetAcknowledgements() const;

  void SetContributor(const std::string &contributor);
  const std::string& GetContributor() const;

  /// Set the type of module: Unknown, SharedObjectModule, CommandLineModule
  void SetType(const std::string &type);

  /// Get the type of the module: Unknown, SharedObjectModule, CommandLineModule
  const std::string& GetType() const;

  /// Set the target for the module.  This is the entry point for a
  /// shared object module and the full command (with path) for an executable.
  void SetTarget(const std::string &target);

  /// Get the target for the module.  This is the entry point for a
  /// shared object module and the full command (with path) for an executable.
  const std::string& GetTarget() const;

  /// Set the location for the module.  This is path to the file (shared
  /// object or executable) for the module.
  void SetLocation(const std::string &target);

  /// Get the location for the module.  This is path to the file (shared
  /// object or executable) for the module.
  const std::string& GetLocation() const;

  void SetLogo(const ModuleLogo& logo);
  const ModuleLogo& GetLogo() const;

  void AddParameterGroup(const ModuleParameterGroup &group);

  const std::vector<ModuleParameterGroup>& GetParameterGroups() const;
  std::vector<ModuleParameterGroup>& GetParameterGroups();

  void SetParameterGroups(const std::vector<ModuleParameterGroup>& groups);

  /// Return true if the module has a parameter matching the \a name.
  /// \sa HasReturnParameters()
  bool HasParameter(const std::string& name) const;

  /// Does the module have any simple (primitive) return types?
  /// \sa HasParameter()
  bool HasReturnParameters() const;

  /// THIS FUNCTION SHOULD NOT BE USED
  /// SEE FindParametersWithValue INSTEAD
  std::vector<ModuleParameter> FindParametersWithDefaultValue(
    const std::string& value)const;

  /// Search the list of parameters and return a copy of the parameters
  /// that have the same \a Value.
  /// \sa HasParameter(), HasReturnParameters(), GetParameterValue()
  std::vector<ModuleParameter> FindParametersWithValue(
    const std::string& value)const;

  /// THIS FUNCTION SHOULD NOT BE USED
  /// SEE SetParameterValue INSTEAD
  bool SetParameterDefaultValue(const std::string& name,
                         const std::string& value);

  /// Set the value of the parameter \a name.
  /// Return true if the parameter is found and different than \a value,
  /// false otherwise.
  /// \sa GetParameterValue(),
  /// \sa FindParametersWithValue(), HasParameter()
  bool SetParameterValue(const std::string& name,
                         const std::string& value);

  /// THIS FUNCTION SHOULD NOT BE USED
  /// SEE GetParameterValue INSTEAD
  std::string GetParameterDefaultValue(const std::string& name) const;

  /// Return the parameter value and an empty string if the parameter
  /// can not be found.
  /// \sa SetParameterValue()
  std::string GetParameterValue(const std::string& name) const;

  const ModuleProcessInformation* GetProcessInformation() const;
  ModuleProcessInformation* GetProcessInformation();

  /// Read a parameter file. Syntax of file is "name: value" for each
  /// parameter. Returns a bool indicating whether any parameter value
  /// was modified.
  bool ReadParameterFile(const std::string& filename);

  /// Write a parameter file. By default, the method writes out all
  /// the parameters.  The "withHandlesToBulkParameters" parameter
  /// controls whether the handles to the bulk parameters (image,
  /// geometry, etc.) are writte to the file.
  bool WriteParameterFile(const std::string& filename, bool withHandlesToBulkParameters = true);

  /// Define a function used to lazily load the CLI shared library
  /// and resolve symbols.
  typedef void (*TargetCallbackType)(void* libraryLoader, ModuleDescription& md);

  /// Set callback function allowing \a libraryLoader to load the CLI shared
  /// library, resolve symbols and set ModuleDescription properties.
  ///
  /// The module descrition properties expected to be set are the \a Target
  /// and the \a Logo.
  ///
  /// The \a Target is the entrypoint, it is usually defined as a string
  /// of the form "slicer:%p" where "%p" is the address of the
  /// \a ModuleEntryPoint symbol.
  ///
  /// \sa TargetCallbackType, GetLibraryLoader(), GetTargetCallback()
  /// \sa SetTarget(const std::string &target), SetLogo(const ModuleLogo& logo)
  void SetTargetCallback(void* libraryLoader, TargetCallbackType targetCallback);

  void* GetLibraryLoader()const;
  TargetCallbackType GetTargetCallback();

private:
  std::string Title;
  std::string Category;
  std::string Index;
  std::string Description;
  std::string Version;
  std::string DocumentationURL;
  std::string License;
  std::string Acknowledgements;
  std::string Contributor;
  std::string Type;
  std::string Target;
  std::string Location;
  std::vector<ModuleParameterGroup> ParameterGroups;

  ModuleProcessInformation ProcessInformation;
  ModuleLogo Logo;

  void* LibraryLoader;
  TargetCallbackType TargetCallback;
};

ModuleDescriptionParser_EXPORT std::ostream & operator<<(std::ostream &os, const ModuleDescription &module);

#endif
