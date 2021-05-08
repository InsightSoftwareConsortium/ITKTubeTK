/*=========================================================================

  Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   Module Description Parser
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/
#include "ModuleDescriptionParser.h"

#include "ModuleDescription.h"
#include "ModuleDescriptionUtilities.h"

#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <stack>
#include "expat.h"

/*********************
 * Utility procedures for strings
 *********************/

//----------------------------------------------------------------------------
static bool
startsWithValidVariableChar(std::string &s)
{
  return (s.find_first_of("_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") == 0);
}

//----------------------------------------------------------------------------
static bool
validVariable(std::string &s)
{
  return (startsWithValidVariableChar(s) && s.find_first_not_of("_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") == std::string::npos);
}

//----------------------------------------------------------------------------
/* ParserState: A class to keep state information for the parser. This
 * is passed to the expat code as user data.
 */
class ParserState
{
public:
  XML_Parser Parser;                           /* The XML parser */
  std::vector<std::string> LastData;           /* The last tag
                                                * processed by expat
                                                * at each depth */
  ModuleDescription CurrentDescription;
  ModuleParameterGroup *CurrentGroup;          /* The parameter group */
  ModuleParameter *CurrentParameter;           /* The current parameter */
  std::stack<std::string> OpenTags;            /* The current tag */
  bool Debug;                                  /* Debug flag */
  bool Error;                                  /* Error detected */
  std::string ErrorDescription;                /* Error description */
  int ErrorLine;                               /* Error line number */
  int Depth;                                   /* The depth of the tag */

  ParserState():
    LastData(10),
    CurrentGroup(0),
    CurrentParameter(0),
    Debug(false),
    Error(false),
    ErrorLine(-1),
    Depth(-1)
  {}
};

//----------------------------------------------------------------------------
/***************************
 * expat callbacks to process the XML
 ***************************/
void
startElement(void *userData, const char *element, const char **attrs)
{
  ParserState *ps = reinterpret_cast<ParserState *>(userData);
  ModuleParameter *parameter = ps->CurrentParameter;
  ModuleParameterGroup *group = ps->CurrentGroup;
  std::string name(element);

  ps->Depth++;
  // Clear the last tag for this depth
  ps->LastData[ps->Depth].clear();

  // Check for a valid module description file
  //
  if (ps->Depth == 0 && (name != "executable") )
    {
    std::string error("ModuleDescriptionParser Error: <executable> must be the outer most tag. <" + name + std::string("> was found instead."));
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
    ps->CurrentParameter = 0;
    ps->CurrentGroup = 0;
    ps->OpenTags.push(name);
    return;
    }
  else if (ps->Depth != 0 && (name == "executable"))
    {
    std::string error("ModuleDescriptionParser Error: <executable> was found inside another tag <" + ps->OpenTags.top() + ">.");
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
    ps->CurrentParameter = 0;
    ps->CurrentGroup = 0;
    ps->OpenTags.push(name);
    return;
    }

  if (name == "parameters")
    {
    if (ps->OpenTags.top() != "executable")
      {
      std::string error("ModuleDescriptionParser Error: <parameters> can only be inside <executable> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->CurrentParameter = 0;
      ps->CurrentGroup = 0;
      ps->OpenTags.push(name);
      return;
      }
    else
      {
      group = new ModuleParameterGroup;
      int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
      if (attrCount == 2 && (strcmp(attrs[0], "advanced") == 0))
        {
        group->SetAdvanced(attrs[1]);
        }
      }
    }
  else if (name == "integer")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + std::string(name) + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
      ps->OpenTags.push(name);
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("int");
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      // check hidden
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      // check multiple
      if ((strcmp(attrs[2*attr], "multiple") == 0) &&
          (strcmp(attrs[2*attr+1], "true") == 0))
        {
        parameter->SetMultiple(attrs[2*attr+1]);
        parameter->SetCPPType("std::vector<int>");
        parameter->SetArgType("int");
        parameter->SetStringToType("atoi");
        }
      }
    }
  else if (name == "float")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + std::string(name) + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
      ps->OpenTags.push(name);
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("float");
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      // check hidden
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      // check multiple
      if ((strcmp(attrs[2*attr], "multiple") == 0) &&
          (strcmp(attrs[2*attr+1], "true") == 0))
        {
        parameter->SetMultiple(attrs[2*attr+1]);
        parameter->SetCPPType("std::vector<float>");
        parameter->SetArgType("float");
        parameter->SetStringToType("atof");
        }
      }
    }
  else if (name == "double")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + std::string(name) + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("double");
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      // check hidden
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      // check multiple
      if ((strcmp(attrs[2*attr], "multiple") == 0) &&
          (strcmp(attrs[2*attr+1], "true") == 0))
        {
        parameter->SetMultiple(attrs[2*attr+1]);
        parameter->SetCPPType("std::vector<double>");
        parameter->SetArgType("double");
        parameter->SetStringToType("atof");
        }
      }
    }
  else if (name == "string")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
      ps->OpenTags.push(name);
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("std::string");
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      // check hidden
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      // check multiple
      if ((strcmp(attrs[2*attr], "multiple") == 0) &&
          (strcmp(attrs[2*attr+1], "true") == 0))
        {
        parameter->SetMultiple(attrs[2*attr+1]);
        parameter->SetCPPType("std::vector<std::string>");
        parameter->SetArgType("std::string");
        parameter->SetStringToType("");
        }
      }
    }
  else if (name == "boolean")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
      ps->OpenTags.push(name);
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("bool");

    // Parse attribute pairs
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      // check hidden
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    }
  else if (name == "integer-vector")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("std::vector<int>");
    parameter->SetArgType("int");
    parameter->SetStringToType("atoi");
    }
  else if (name == "float-vector")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("std::vector<float>");
    parameter->SetArgType("float");
    parameter->SetStringToType("atof");
    }
  else if (name == "string-vector")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("std::vector<std::string>");
    parameter->SetArgType("std::string");
    parameter->SetStringToType("");
    }
  else if (name == "double-vector")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("std::vector<double>");
    parameter->SetArgType("double");
    parameter->SetStringToType("atof");
    }
  else if (name == "point")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    // Parse attribute pairs
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          parameter->SetCPPType("std::vector<std::vector<float> >");
          parameter->SetArgType("float");
          parameter->SetStringToType("atof");
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "coordinateSystem") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "ijk") == 0) ||
            (strcmp(attrs[2*attr+1], "lps") == 0) ||
            (strcmp(attrs[2*attr+1], "ras") == 0))
          {
          parameter->SetCoordinateSystem(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid coordinate system. Only \"ijk\", \"lps\" and \"ras\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else
        {
          std::string error("ModuleDescriptionParser Error: " + std::string(attrs[2*attr]) + " is not a valid attribute for the tag" + name);
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
        }
      }
    if (parameter->GetMultiple() != "true")
      {
      parameter->SetCPPType("std::vector<float>");
      parameter->SetArgType("float");
      parameter->SetStringToType("atof");
      }
    parameter->SetTag(name);
    }
  else if (name == "pointfile")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("point");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          if (strcmp(attrs[2*attr+1], "true") == 0)
            {
            parameter->SetCPPType("std::vector<std::string>");
            parameter->SetArgType("std::string");
            }
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "type") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "any") == 0) ||
            (strcmp(attrs[2*attr+1], "point") == 0) ||
            (strcmp(attrs[2*attr+1], "line") == 0) ||
            (strcmp(attrs[2*attr+1], "angle") == 0) ||
            (strcmp(attrs[2*attr+1], "curve") == 0) ||
            (strcmp(attrs[2*attr+1], "closedcurve")) )
          {
          parameter->SetType(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid value for the attribute \"" + "type" + "\". Only \"point\", \"line\" , \"angle\", \"curve\", \"closedcurve\", and \"any\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else if ((strcmp(attrs[2*attr], "coordinateSystem") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "ijk") == 0) ||
            (strcmp(attrs[2*attr+1], "lps") == 0) ||
            (strcmp(attrs[2*attr+1], "ras") == 0))
          {
          parameter->SetCoordinateSystem(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid coordinate system. Only \"ijk\", \"lps\" and \"ras\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\" and \"fileExtensions\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    }
  else if (name == "region")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    // Parse attribute pairs
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          parameter->SetCPPType("std::vector<std::vector<float> >");
          parameter->SetArgType("float");
          parameter->SetStringToType("atof");
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "coordinateSystem") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "ijk") == 0) ||
            (strcmp(attrs[2*attr+1], "lps") == 0) ||
            (strcmp(attrs[2*attr+1], "ras") == 0))
          {
          parameter->SetCoordinateSystem(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid coordinate system. Only \"ijk\", \"lps\" and \"ras\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else
        {
          std::string error("ModuleDescriptionParser Error: " + std::string(attrs[2*attr]) + " is not a valid attribute for the tag" + name);
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
        }
      }
    if (parameter->GetMultiple() != "true")
      {
      parameter->SetCPPType("std::vector<float>");
      parameter->SetArgType("float");
      parameter->SetStringToType("atof");
      }
    parameter->SetTag(name);
    }
  else if (name == "string-enumeration")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("std::string");
    // Parse attribute pairs
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    }
  else if (name == "integer-enumeration")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("int");
    // Parse attribute pairs
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    }
  else if (name == "float-enumeration")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("float");
    // Parse attribute pairs
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    }
  else if (name == "double-enumeration")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    parameter->SetTag(name);
    parameter->SetCPPType("double");
    // Parse attribute pairs
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      bool res = true;
      if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        res = ModuleDescriptionParser::processHiddenAttribute(attrs[2*attr+1], parameter, ps);
        }
      if (!res)
        {
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    }
  else if (name == "file")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("scalar");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          if (strcmp(attrs[2*attr+1], "true") == 0)
            {
            parameter->SetCPPType("std::vector<std::string>");
            parameter->SetArgType("std::string");
            }
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\" and \"fileExtensions\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    }
  else if (name == "directory")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    if (attrCount == 2 &&
        (strcmp(attrs[0], "multiple") == 0) &&
        (strcmp(attrs[1], "true") == 0))
      {
      parameter->SetMultiple(attrs[1]);
      parameter->SetCPPType("std::vector<std::string>");
      parameter->SetArgType("std::string");
      }
    else
      {
      parameter->SetCPPType("std::string");
      }
    parameter->SetTag(name);
    }
  else if (name == "transform")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("unknown");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          if (strcmp(attrs[2*attr+1], "true") == 0)
            {
            parameter->SetCPPType("std::vector<std::string>");
            parameter->SetArgType("std::string");
            }
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "type") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "linear") == 0) ||
            (strcmp(attrs[2*attr+1], "nonlinear") == 0) ||
            (strcmp(attrs[2*attr+1], "bspline") == 0))
          {
          parameter->SetType(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid value for the attribute \"" + "type" + "\". Only \"linear\", \"bspline\" and \"nonlinear\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else if ((strcmp(attrs[2*attr], "reference") == 0))
        {
        parameter->SetReference(attrs[2*attr+1]);
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\", \"fileExtensions\" and \"type\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    }
  else if (name == "image")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("scalar");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          if (strcmp(attrs[2*attr+1], "true") == 0)
            {
            parameter->SetCPPType("std::vector<std::string>");
            parameter->SetArgType("std::string");
            }
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "type") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "any") == 0) ||
            (strcmp(attrs[2*attr+1], "scalar") == 0) ||
            (strcmp(attrs[2*attr+1], "label") == 0) ||
            (strcmp(attrs[2*attr+1], "tensor") == 0) ||
            (strcmp(attrs[2*attr+1], "diffusion-weighted") == 0) ||
            (strcmp(attrs[2*attr+1], "vector") == 0) ||
            (strcmp(attrs[2*attr+1], "signal") == 0) ||
            (strcmp(attrs[2*attr+1], "multichannel") == 0) ||
            (strcmp(attrs[2*attr+1], "dynamic-contrast-enhanced") == 0))
          {
          parameter->SetType(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid value for the attribute \"" + "type" + "\". Only \"scalar\", \"label\" , \"tensor\", \"diffusion-weighted\", \"vector\", \"signal\", \"multichannel\", \"dynamic-contrast-enhanced\" and \"any\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetHidden(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"hidden\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "reference") == 0))
        {
        parameter->SetReference(attrs[2*attr+1]);
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\", \"fileExtensions\", \"type\", \"hidden\" and \"reference\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    }
  else if (name == "geometry")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("scalar");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "aggregate") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetAggregate(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"aggregate\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }      else if ((strcmp(attrs[2*attr], "type") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "fiberbundle") == 0) ||
            (strcmp(attrs[2*attr+1], "model") == 0))
          {
          parameter->SetType(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid value for the attribute \"" + "type" + "\". Only \"fiberbundle\", \"model\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else if ((strcmp(attrs[2*attr], "reference") == 0))
        {
        parameter->SetReference(attrs[2*attr+1]);
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\", \"type\", and \"fileExtensions\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    if (parameter->GetMultiple() == "true"
        && parameter->GetAggregate() != "true")
      {
      parameter->SetCPPType("std::vector<std::string>");
      parameter->SetArgType("std::string");
      }
    }
  else if (name == "table")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("scalar");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          if (strcmp(attrs[2*attr+1], "true") == 0)
            {
            parameter->SetCPPType("std::vector<std::string>");
            parameter->SetArgType("std::string");
            }
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "reference") == 0))
        {
        parameter->SetReference(attrs[2*attr+1]);
        }
      else if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetHidden(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"hidden\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "type") == 0))
        {
        // Need to add other tables with context (similar to color
        // tables) as well as add general CSV files
        if ((strcmp(attrs[2*attr+1], "color") == 0))
          {
          parameter->SetType(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid value for the attribute \"" + "type" + "\". Only \"color\" is currently accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\", \"hidden\", \"reference\", \"type\", and \"fileExtensions\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    }
  else if (name == "measurement")
    {
    if (!group || (ps->OpenTags.top() != "parameters"))
      {
      std::string error("ModuleDescriptionParser Error: <" + name + "> can only be used inside <parameters> but was found inside <" + ps->OpenTags.top() + ">");
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      ps->OpenTags.push(name);
      delete parameter;
      return;
      }
    parameter = new ModuleParameter;
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);

    // Parse attribute pairs
    parameter->SetCPPType("std::string");
    parameter->SetType("scalar");
    for (int attr=0; attr < (attrCount / 2); attr++)
      {
      if ((strcmp(attrs[2*attr], "multiple") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetMultiple(attrs[2*attr+1]);
          if (strcmp(attrs[2*attr+1], "true") == 0)
            {
            parameter->SetCPPType("std::vector<std::string>");
            parameter->SetArgType("std::string");
            }
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"multiple\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "reference") == 0))
        {
        parameter->SetReference(attrs[2*attr+1]);
        }
      else if ((strcmp(attrs[2*attr], "hidden") == 0))
        {
        if ((strcmp(attrs[2*attr+1], "true") == 0) ||
            (strcmp(attrs[2*attr+1], "false") == 0))
          {
          parameter->SetHidden(attrs[2*attr+1]);
          }
        else
          {
          std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid argument for the attribute \"hidden\". Only \"true\" and \"false\" are accepted.");
          if (ps->ErrorDescription.size() == 0)
            {
            ps->ErrorDescription = error;
            ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
            ps->Error = true;
            }
          ps->OpenTags.push(name);
          delete parameter;
          return;
          }
        }
      else if ((strcmp(attrs[2*attr], "type") == 0))
        {
        // "type" will refer to the error models (regression,
        // orthogonal regression, ...). We'll define these later. When
        // we do, put the error handling back in.
//         if ((strcmp(attrs[2*attr+1], "some_error_type") == 0))
//           {
          parameter->SetType(attrs[2*attr+1]);
//           }
//         else
//           {
//           std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr+1]) + "\" is not a valid value for the attribute \"" + "type" + "\". Only \"color\" is currently accepted.");
//           if (ps->ErrorDescription.size() == 0)
//             {
//             ps->ErrorDescription = error;
//             ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
//             ps->Error = true;
//             }
//           ps->OpenTags.push(name);
//           delete parameter;
//           return;
//           }
        }
      else if ((strcmp(attrs[2*attr], "fileExtensions") == 0))
        {
        parameter->SetFileExtensionsAsString(attrs[2*attr+1]);
        }
      else
        {
        std::string error("ModuleDescriptionParser Error: \"" + std::string(attrs[2*attr]) + "\" is not a valid attribute for \"" + name + "\". Only \"multiple\", \"hidden\", \"reference\", \"type\", and \"fileExtensions\" are accepted.");
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        ps->OpenTags.push(name);
        delete parameter;
        return;
        }
      }
    parameter->SetTag(name);
    }
  else if (parameter && (name == "flag"))
    {
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    if (attrCount > 0)
      {
      for (int i=0; i < attrCount/2; ++i)
        {
        if (strcmp(attrs[2*i], "alias") == 0)
          {
          parameter->SetFlagAliasesAsString(attrs[2*i+1]);
          }
        else if (strcmp(attrs[2*i], "deprecatedalias") == 0)
          {
          parameter->SetDeprecatedFlagAliasesAsString(attrs[2*i+1]);
          }
        }
      }
    }
  else if (parameter && (name == "longflag"))
    {
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    if (attrCount > 0)
      {
      for (int i=0; i < attrCount/2; ++i)
        {
        if (strcmp(attrs[2*i], "alias") == 0)
          {
          parameter->SetLongFlagAliasesAsString(attrs[2*i+1]);
          }
        else if (strcmp(attrs[2*i], "deprecatedalias") == 0)
          {
          parameter->SetDeprecatedLongFlagAliasesAsString(attrs[2*i+1]);
          }
        }
      }
    }
  else if (parameter && (name == "reference"))
    {
    int attrCount = XML_GetSpecifiedAttributeCount(ps->Parser);
    if (attrCount > 0)
      {
      std::string role("");
      std::string parameterName("");
      for (int i=0; i < attrCount/2; ++i)
        {
        if (strcmp(attrs[2*i], "role") == 0)
          {
          role = std::string(attrs[2*i+1]);
          }
        else if (strcmp(attrs[2*i], "parameter") == 0)
          {
          parameterName = std::string(attrs[2*i+1]);
          }
        }
      parameter->AddForwardReference(role, parameterName);
      }
    }
  ps->CurrentParameter = parameter;
  ps->CurrentGroup = group;
  ps->OpenTags.push(name);

}

//----------------------------------------------------------------------------
void
endElement(void *userData, const char *element)
{
  ParserState *ps = reinterpret_cast<ParserState *>(userData);
  ModuleParameter *parameter = ps->CurrentParameter;
  ModuleParameterGroup *group = ps->CurrentGroup;
  std::string name(element);

  if (name == "parameters" && ps->Depth == 1)
    {
    ps->CurrentDescription.AddParameterGroup(*ps->CurrentGroup);
    if (ps->CurrentGroup)
      {
      delete ps->CurrentGroup;
      }
    ps->CurrentGroup = 0;
    ps->CurrentParameter = 0;
    }
  else if (name == "parameters" && ps->Depth != 1)
    {
    ps->CurrentGroup = 0;
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "integer"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "float"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "double"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "string"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name =="boolean"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "file"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "directory"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "transform"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "image"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "geometry"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "table"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "measurement"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "integer-vector"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "float-vector"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "string-vector"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "double-vector"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "point"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "pointfile"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "region"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "string-enumeration"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "integer-enumeration"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "float-enumeration"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (group && parameter && (name == "double-enumeration"))
    {
    ps->CurrentGroup->AddParameter(*parameter);
    if (ps->CurrentParameter)
      {
      delete ps->CurrentParameter;
      }
    ps->CurrentParameter = 0;
    }
  else if (parameter && (name == "flag"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    trimLeading(temp, "-");
    if (temp.size() > 1)
      {
      std::string error("ModuleDescriptionParser Error: <"
                        + name
                        + std::string("> can only contain one character. \"")
                        + temp
                        + std::string("\" has more than one character."));
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      if (!ps->OpenTags.empty())
        {
        ps->OpenTags.pop();
        }
      ps->Depth--;
      return;
      }
    else
      {
      if (!parameter->GetIndex().empty())
        {
        std::string error("ModuleDescriptionParser Error: <"
                          + std::string(name)
                          + "> cannot be specified because an index has been specified for this parameter."
                          + std::string("\""));
        if (ps->ErrorDescription.size() == 0)
          {
          ps->ErrorDescription = error;
          ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
          ps->Error = true;
          }
        if (!ps->OpenTags.empty())
          {
          ps->OpenTags.pop();
          ps->Depth--;
          }
        return;
        }
      parameter->SetFlag(temp);
      }
    }
  else if (parameter && (name == "longflag"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    trimLeading(temp, "-");
    if (!validVariable(temp))
      {
      std::string error("ModuleDescriptionParser Error: <"
                        + std::string(name)
                        + "> can only contain letters, numbers and underscores and must start with a _ or letter. The offending name is \""
                        + temp
                        + std::string("\""));
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      if (!ps->OpenTags.empty())
        {
        ps->OpenTags.pop();
        ps->Depth--;
        }
      return;
      }
    if (!parameter->GetIndex().empty())
      {
      std::string error("ModuleDescriptionParser Error: <"
                        + std::string(name)
                        + "> cannot be specified because an index has been specified for this parameter."
                        + std::string("\""));
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
      if (!ps->OpenTags.empty())
        {
        ps->OpenTags.pop();
        ps->Depth--;
        }
      return;
      }
    parameter->SetLongFlag(temp);
    if (parameter->GetName().empty())
      {
      parameter->SetName(temp);
      }
    }
  else if (parameter && (name == "name"))
    {
    std::string temp = std::string(ps->LastData[ps->Depth]);
    trimLeadingAndTrailing(temp);
    if (!validVariable(temp))
      {
      std::string error("ModuleDescriptionParser Error: <"
                        + std::string(name)
                        + "> can only contain letters, numbers and underscores and must start with an _ or letter. The offending name is \""
                        + temp
                        + std::string("\""));
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
        if (!ps->OpenTags.empty())
          {
          ps->OpenTags.pop();
          ps->Depth--;
          }
      return;
      }
    parameter->SetName(temp);
    }
  else if ((group || parameter) && (name == "label"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    if (group && !parameter)
      {
      group->SetLabel(temp);
      }
    else
      {
      parameter->SetLabel(temp);
      }
    }
  else if (name == "category")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetCategory(temp);
    }
  else if (!parameter && (name == "index"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetIndex(temp);
    }
  else if (name == "title")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetTitle(temp);
    }
  else if (name == "version")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetVersion(temp);
    }
  else if (name == "documentation-url")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetDocumentationURL(temp);
    }
  else if (name == "license")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetLicense(temp);
    }
  else if (name == "acknowledgements")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetAcknowledgements(temp);
    }
  else if (name == "contributor")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetContributor(temp);
    }
  else if (name == "location")
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    ps->CurrentDescription.SetLocation(temp);
    }
  else if (name ==  "description")
    {
    std::string temp = ps->LastData[ps->Depth];
    replaceSubWithSub(temp, "\"", "'");
    replaceSubWithSub(temp, "\n", "\\n");
    trimLeadingAndTrailing(temp);
    if (!group && !parameter)
      {
      ps->CurrentDescription.SetDescription(temp);
      }
    else if (group && !parameter)
      {
      group->SetDescription(temp);
      }
    else
      {
      parameter->SetDescription(temp);
      }
    }
  else if (parameter && (name == "element"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->GetElements().push_back(temp);
    }
  else if (parameter && (name == "default"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->SetValue(temp);
    }
  else if (parameter && (name == "channel"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->SetChannel(temp);
    }
  else if (parameter && (name == "index"))
    {
    if (!parameter->GetLongFlag().empty() || !parameter->GetFlag().empty() )
      {
      std::string error("ModuleDescriptionParser Error: <"
                        + std::string(name)
                        + "> cannot be specified because a <longflag> and/or <flag> has been specified for this parameter."
                        + std::string("\""));
      if (ps->ErrorDescription.size() == 0)
        {
        ps->ErrorDescription = error;
        ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
        ps->Error = true;
        }
        if (!ps->OpenTags.empty())
          {
          ps->OpenTags.pop();
          ps->Depth--;
          }
      return;
      }
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->SetIndex(temp);
    }
  else if (parameter && (name == "constraints"))
    {
    parameter->SetConstraints(name);
    }
  else if (parameter && (name == "minimum"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->SetMinimum(temp);
    }
  else if (parameter && (name == "maximum"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->SetMaximum(temp);
    }
  else if (parameter && (name == "step"))
    {
    std::string temp = ps->LastData[ps->Depth];
    trimLeadingAndTrailing(temp);
    parameter->SetStep(temp);
    }
  else if (parameter && (name == "reference"))
    {
    // No-op: merely make "reference" element accepted. Its contents are handled in startElement
    }
  else if (name != "executable")
    {
    std::string error("ModuleDescriptionParser Error: Unrecognized element <" + name + std::string("> was found."));
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
    }

  if (!ps->OpenTags.empty())
    {
    ps->OpenTags.pop();
    ps->Depth--;
    }
}

//----------------------------------------------------------------------------
bool
ModuleDescriptionParser::processHiddenAttribute(const char* value, ModuleParameter* parameter, ParserState* ps)
{
  if ((strcmp(value, "true") == 0) ||
      (strcmp(value, "false") == 0))
    {
    parameter->SetHidden(value);
    return true;
    }
  else
    {
    std::string error("ModuleDescriptionParser Error: \"" + std::string(value) + "\" is not a valid argument for the attribute \"hidden\". Only \"true\" and \"false\" are accepted.");
    if (ps->ErrorDescription.size() == 0)
      {
      ps->ErrorDescription = error;
      ps->ErrorLine = XML_GetCurrentLineNumber(ps->Parser);
      ps->Error = true;
      }
    return false;
    }
}

//----------------------------------------------------------------------------
void
charData(void *userData, const char *s, int len)
{
  ParserState *ps = reinterpret_cast<ParserState *>(userData);
  if (len)
    {
    std::string str(s,len);
    ps->LastData[ps->Depth] += str;
    }
}

//----------------------------------------------------------------------------
int
ModuleDescriptionParser::Parse( const std::string& xml, ModuleDescription& description)
{
// Check the first line if the XML
  if (strncmp(xml.c_str(),"<?xml ", 6) != 0)
    {
    std::string required("<?xml version=\"1.0\" encoding=\"utf-8\"?>");
    std::cerr << "ModuleDescriptionParser: first line must be " << std::endl;
    std::cerr << required << std::endl;
    return 1;
    }

  ParserState parserState;
  parserState.CurrentDescription = description;

  XML_Parser parser = XML_ParserCreate(NULL);
  int done;

  parserState.Parser = parser;
  parserState.CurrentParameter = 0;
  parserState.CurrentGroup = 0;

  XML_SetUserData(parser, static_cast<void *>(&parserState));
  XML_SetElementHandler(parser, startElement, endElement);
  XML_SetCharacterDataHandler(parser, charData);

  // Parse the XML
  done = true;
  int status = 0;
  if (XML_Parse(parser, xml.c_str(), static_cast< int >(xml.size()), done) == 0)
    {
    std::cerr << XML_ErrorString(XML_GetErrorCode(parser))
              << " at line "
              << XML_GetCurrentLineNumber(parser)
              << std::endl;
    status = 1;
    }
  else
    {
    if (parserState.Error)
      {
      std::cerr << parserState.ErrorDescription
                << " at line "
                <<  parserState.ErrorLine
                << " while parsing "
                << xml
                << std::endl;
      status = 1;
      }
    }

  XML_ParserFree(parser);

  description = parserState.CurrentDescription;
  return status;
}
