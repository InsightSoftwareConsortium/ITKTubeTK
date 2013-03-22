/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: OptionList.h,v $
  Language:  C++
  Date:      $Date: 2004/05/06 13:45:52 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __OPTIONLIST_H_
#define __OPTIONLIST_H_

#include <itkWin32Header.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "itkArray.h"

namespace tube
{

class OptionList
{

  public:

    typedef std::multimap<std::string, std::string> OptionMap ;
    typedef std::vector<std::string>                StringVector ;
    typedef const char *                            LabelType;

    static const char *                   INPUT_DELIMITER;
    static const char *                   BOOLEAN_VALUE_TRUE;
    static const char *                   BOOLEAN_VALUE_FALSE;
    static LabelType                      NUMBER_SEPARATOR;
    static LabelType                      DUAL_ELEMENT_DELIMITER;


    OptionList(int argc, char* argv[]) ;
    ~OptionList() {}

    class RequiredOptionMissing
    {
      public:

        RequiredOptionMissing(std::string tag) : OptionTag( tag ) {}
        std::string OptionTag ;
    };

    int GetOption(std::string option_tag, StringVector* values) ;

    int DumpOption(std::string option_tag, bool withTag = true, bool withNewLine = false) ;

    int GetMultiDoubleOption(std::string tag, std::vector<double>* args, bool required) ;

    int GetMultiDoubleOption(std::string tag, itk::Array<double>* args, bool required) ;

    double GetDoubleOption(std::string tag, double default_value, bool required);

    bool GetBooleanOption(std::string tag, bool default_value, bool required);

    int GetMultiIntOption(std::string tag, std::vector<int>* args, bool required );

    int GetIntOption(std::string tag, int default_value, bool required) ;

    int GetMultiUCharOption(std::string tag, std::vector< unsigned char >* args, bool required);

    unsigned char GetUCharOption(std::string tag, unsigned char default_value, bool required) ;

    int GetMultiUIntOption(std::string tag, std::vector< unsigned int >* args, bool required);

    unsigned int GetUIntOption(std::string tag, unsigned int default_value, bool required) ;

    int GetStringOption(std::string tag, std::string* ret, bool required);

    int GetMultiStringOption(std::string tag, std::vector< std::string >* ret,
                             bool required);
  protected:

  private:

    OptionMap m_Map;
};

} // End namespace tube

#endif
