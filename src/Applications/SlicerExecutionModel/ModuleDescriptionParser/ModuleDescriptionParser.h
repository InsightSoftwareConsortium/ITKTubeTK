/*=========================================================================

  Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   Module Description Parser
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/
#ifndef __ModuleDescriptionParser_h
#define __ModuleDescriptionParser_h

#include <string>

#include "ModuleDescriptionParserExport.h"

class ModuleDescription;
class ModuleParameter;
class ModuleParameterGroup;
class ParserState;

class ModuleDescriptionParser_EXPORT ModuleDescriptionParser
{
public:
  ModuleDescriptionParser() {}
  ~ModuleDescriptionParser() {}

  int Parse( const std::string& xml, ModuleDescription& description);

  static bool processHiddenAttribute(const char* value, ModuleParameter* parameter, ParserState* ps);
};

#endif
