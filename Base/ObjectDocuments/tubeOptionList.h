/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: tubeOptionList.h,v $
  Language:  C++
  Date:      $Date: 2004/05/06 13:45:52 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __tubeOptionList_h
#define __tubeOptionList_h

#include <itkArray.h>

#include <map>
#include <string>
#include <vector>

namespace tube
{

class OptionList
{

public:

  typedef std::multimap<std::string, std::string> OptionMap;
  typedef std::vector<std::string>                StringVector;
  typedef std::string                             LabelType;

  static const std::string              INPUT_DELIMITER;
  static const std::string              BOOLEAN_VALUE_TRUE;
  static const std::string              BOOLEAN_VALUE_FALSE;
  static LabelType                      NUMBER_SEPARATOR;
  static LabelType                      DUAL_ELEMENT_DELIMITER;


  OptionList( int argc, char * argv[] );
  ~OptionList( void ) {}

  class RequiredOptionMissing
    {
    public:
      RequiredOptionMissing(const std::string & tag) : OptionTag( tag ) {}
      std::string OptionTag;

    }; // End class RequiredOptionMissing

  int GetOption(const std::string & option_tag, StringVector & values) const;

  int DumpOption(const std::string & option_tag, bool withTag = true, bool withNewLine = false) const;

  int GetMultiDoubleOption(const std::string & tag, std::vector<double> & args, bool required) const;

  int GetMultiDoubleOption(const std::string & tag, itk::Array<double> & args, bool required) const;

  double GetDoubleOption(const std::string & tag, double default_value, bool required) const;

  bool GetBooleanOption(const std::string & tag, bool default_value, bool required) const;

  int GetMultiIntOption(const std::string & tag, std::vector<int> & args, bool required ) const;

  int GetIntOption(const std::string & tag, int default_value, bool required) const;

  int GetMultiUCharOption(const std::string & tag, std::vector< unsigned char > & args, bool required) const;

  unsigned char GetUCharOption(const std::string & tag, unsigned char default_value, bool required) const;

  int GetMultiUIntOption(const std::string & tag, std::vector< unsigned int > & args, bool required) const;

  unsigned int GetUIntOption(const std::string & tag, unsigned int default_value, bool required) const;

  int GetStringOption(const std::string & tag, std::string & ret, bool required) const;

  int GetMultiStringOption(const std::string & tag, std::vector< std::string > & ret,
                           bool required) const;
protected:

private:

  OptionMap m_Map;

}; // End class OptionList

} // End namespace tube

#endif // End !defined(__tubeOptionList_h)
