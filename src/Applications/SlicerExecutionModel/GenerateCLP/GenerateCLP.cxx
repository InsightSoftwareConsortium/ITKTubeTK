/*=========================================================================

  Portions (c) Copyright 2006 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   GenerateCLP
  Module:    $URL: http://svn.na-mic.org:8000/svn/NAMICSandBox/trunk/CommandLineAPI/GenerateCLP.cxx $
  Date:      $Date: 2006-04-21 15:51:08 -0400 (Fri, 21 Apr 2006) $
  Version:   $Revision: 957 $

=========================================================================*/

/* Generate command line processing code from an xml description
   Usage: GenerateCLP input_xml_file output_include_file

   This program generates source code that processes command line
   arguments. The arguments are described in an xml file. The output
   include file defines a macro PARSE_ARGS. This macro generates the
   declarations for each command line argument and generates code to
   parse the command line. In addition to user specified coammnd line
   arguments, code to echo the xml file and code to print the command
   line arguments is also generated.

   Typical usage is:
   GenerateCLP foo.xml fooCLP.h

   foo.cxx contains:
   #include fooCLP.h
   int main (int argc, char *argv[])
   {
     PARSE_ARGS;
         .
         . Code to implement foo
         .
   }

   This first version of the code uses the TCLAP Templated C++ Command
   Line Parser (http://tclap.sourceforge.net). TCLAP is attractive
   because it is implemented entirely in header files. Other command line
   parsing may be supported in the future.

  GenerateCLP uses the expat XML parser (http://expat.sourceforge.net)
  to process the XML file.

  The generated C++ code relies on the kwsys library
  (http://www.cmake.org) to provide a portable implementaion of string
  streams.
*/

#include "GenerateCLPConfig.h" // For GENERATECLP_USE_MD5, GenerateCLP_USE_SERIALIZER

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "expat.h"
#include <string>
#include <vector>

#include <itksys/SystemTools.hxx>
#ifdef GENERATECLP_USE_MD5
# include <itksys/MD5.h>
#endif

#include "GenerateCLP.h"
#include "ModuleDescriptionUtilities.h"
#include "ModuleDescriptionParser.h"
#include "ModuleDescription.h"
#include "ModuleParameterGroup.h"
#include "ModuleParameter.h"

#ifdef GenerateCLP_USE_SERIALIZER
# include "itkJsonCppArchiver.h"
# include "itkSEMModuleDescriptionSerializer.h"
# include "ModuleDescriptionParser.h"
#endif

#if defined( _WIN32 )
#pragma warning ( disable : 4996 )
#endif

namespace
{

/* Comma separated arguments need a temporary variable to store the
 * string
 */
bool NeedsTemp(const ModuleParameter &parameter)
{
  std::string type = parameter.GetCPPType();
  std::string multi = parameter.GetMultiple();
  return (((type == "std::vector<int>" ||
           type == "std::vector<float>" ||
           type == "std::vector<double>" ||
           type == "std::vector<std::string>") ||
           type == "std::vector<std::vector<float> >"));
}
/* Some types need quotes in the initialization. */
bool NeedsQuotes(const ModuleParameter &parameter)
{
  std::string type = parameter.GetCPPType();
  std::string multi = parameter.GetMultiple();
  return (((type == "std::vector<int>" ||
           type == "std::vector<float>" ||
           type == "std::vector<double>" ||
           type == "std::vector<std::string>" ||
           type == "std::string") &&
           multi != "true") ||
           type == "std::vector<std::vector<float> >"
);
}
bool IsEnumeration(const ModuleParameter &parameter)
{
  std::string type = parameter.GetTag();
  return (type == "string-enumeration" ||
          type == "integer-enumeration" ||
          type == "float-enumeration" ||
          type == "double-enumeration" );
}

bool IsVectorOfVectors(const ModuleParameter &parameter)
{
  std::string type = parameter.GetCPPType();
  return (type == "std::vector<std::vector<float> >");
}

bool HasValue(const ModuleParameter &parameter)
{
  return (parameter.GetValue().size() > 0 && parameter.GetMultiple() != "true");
}

/* Generate the preamble to the code. This includes the required
 * include files and code to process comma separated arguments.
 */
void GeneratePre(std::ostream &, ModuleDescription &, int, char *[]);

#ifdef GenerateCLP_USE_SERIALIZER
/** Write the JSON description of the module in header. */
void GenerateJson( std::ostream & sout,
  itk::SEMModuleDescriptionSerializer::Pointer & serializer );
#endif // GenerateCLP_USE_SERIALIZER

#ifdef GenerateCLP_USE_JSONCPP
/** Create code to deserialize parameters from a file. */
void GenerateDeSerialization( std::ostream & sout,
  const ModuleDescription & module );

/** Create code to serialize the parameters to a file. */
void GenerateSerialization( std::ostream & sout,
  const ModuleDescription & module );
#endif // GenerateCLP_USE_JSONCPP

/* Generate the last statements. This defines the PARSE_ARGS macro */
void GeneratePost(std::ostream &);

/* Generate a function to split a string into a vector of strings. */
void GenerateSplitString(std::ostream &);

/* Generate a function to split a string into a vector of filenames
 * (which can contain commas in the name). */
void GenerateSplitFilenames(std::ostream &);

/* Generate the code that echos the XML file that describes the
 * command line arguments.
 */
void GenerateExports(std::ostream &);
void GeneratePluginDataSymbols(std::ostream &, std::vector<std::string>&, std::string);
void GeneratePluginEntryPoints(std::ostream &, std::vector<std::string> &);
void GeneratePluginProcedures(std::ostream &, std::vector<std::string> &);
void GenerateLOGO(std::ostream &, std::vector<std::string> &);
void GenerateXML(std::ostream &);

/* Generate the code that uses TCLAP to declare the parameters.
*/
void GenerateDeclare(std::ostream &, ModuleDescription &);

/* Generate the code that uses TCLAP to parse the command line
 * arguments.
 */
void GenerateTCLAP(std::ostream &, ModuleDescription &);

/** Generate the code that parses the command line arguments and puts them into
 * TCLAP data structures. */
void GenerateTCLAPParse(std::ostream & sout, ModuleDescription & module);

/* Generate the code that assigns the TCLAP parsed variables to their C++
 * values.  If onlyIfSet is true, the will only be set values have been passed
 * on the command line. */
void GenerateTCLAPAssignment(std::ostream & sout, const ModuleDescription & module, bool onlyIfSet=false);

/* Generate code to echo the command line arguments and their values. */
void GenerateEchoArgs(std::ostream &, ModuleDescription &);

/** Generate code to decode the address of a process information
 * structure */
void GenerateProcessInformationAddressDecoding(std::ostream &sout);

#ifdef GENERATECLP_USE_MD5
///** Compute the md5sum of a file */
bool ComputeFileMD5(const char* source, char* md5out);

///** Compute the md5sum of a string.  */
std::string ComputeStringMD5(const char* input);
#endif

} // end of anonymous namespace

int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  ModuleDescription module;
  ModuleDescriptionParser parser;

  // Read the XML file
  std::ifstream fin(InputXML.c_str(),std::ios::in|std::ios::binary);
  if (fin.fail())
    {
    std::cerr << argv[0] << ": Cannot open " << InputXML << " for input" << std::endl;
    perror(argv[0]);
    return EXIT_FAILURE;
    }

  // Get the length of the file
  fin.seekg (0, std::ios::end);
  const size_t len = fin.tellg();
  fin.seekg (0, std::ios::beg);
  char * XML = new char[len+1];
  fin.read (XML, len);
  XML[len] = '\0';

  // Parse the module description
  std::cout << "GenerateCLP";
  for (int i = 1; i < argc; i++)
    {
    std::cout << " " << argv[i];
    }
  std::cout << std::endl;

  if (parser.Parse(XML, module))
    {
    std::cerr << "GenerateCLP: One or more errors detected. Code generation aborted." << std::endl;
    return EXIT_FAILURE;
    }

#ifdef GenerateCLP_USE_SERIALIZER
  itk::SEMModuleDescription::Pointer itkModule =
    itk::SEMModuleDescription::New();
  itkModule->SetModuleDescription( module );
  itk::SEMModuleDescriptionSerializer::Pointer serializer =
    itk::SEMModuleDescriptionSerializer::New();
  serializer->SetTargetObject( itkModule );
  serializer->Serialize();
#endif // GenerateCLP_USE_SERIALIZER

  std::stringstream parametersGroupsMsg;

  // Print each command line arg
  parametersGroupsMsg << "GenerateCLP: Found " << module.GetParameterGroups().size() << " parameters groups" << std::endl;
  std::vector<ModuleParameterGroup>::const_iterator git;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    parametersGroupsMsg << "GenerateCLP: Group \"" << (*git).GetLabel() << "\" has " << (*git).GetParameters().size() << " parameters" << std::endl;
    }

  //
  // If GENERATE_USE_MD5 is defined, generated output will be appended into a stringstream.
  //
  // Following the generation, the MD5 of both the output buffer and
  // the existing file (if any) will be computed.
  //
  // Three cases can be numbered:
  //   1) Output file [OutputCxx] do NOT exists               => Write output buffer to file.
  //   2) Output file [OutputCxx] exists and MD5 match        => Discard Output buffer.
  //   3) Output file [OutputCxx] exists and MD5 do NOT match => Write output buffer to file.
  //
  // Computing the MD5 allows to update the file only if required. Since the
  // modified time associated the file will be updated only if it really changed,
  // extra compilation cycle could be avoided.
  //

  // Do the hard stuff
#ifdef GENERATECLP_USE_MD5
  std::stringstream sout;
#else
  std::cout << parametersGroupsMsg.str();
  std::cout.flush();

  std::ofstream sout(OutputCxx.c_str(), std::ios::out);
  if (sout.fail())
    {
    std::cerr << argv[0] << ": Cannot open " << OutputCxx << " for output" << std::endl;
    perror(argv[0]);
    return EXIT_FAILURE;
    }
#endif
  if (logoFiles.size() > 0 && !itksys::SystemTools::FileExists(logoFiles[0].c_str()))
    {
    std::cerr << argv[0] << ": Cannot open " << logoFiles[0] << " as a logo file" << std::endl;
    return EXIT_FAILURE;
    }

  GeneratePre(sout, module, argc, argv);
  GenerateExports(sout);
#ifdef GenerateCLP_USE_SERIALIZER
  // TODO make this create json-dl
  GenerateJson( sout, serializer );
#endif // GenerateCLP_USE_SERIALIZER
#ifdef GenerateCLP_USE_JSONCPP
  GenerateDeSerialization( sout, module );
  GenerateSerialization( sout, module );
#endif // GenerateCLP_USE_JSONCPP
  GeneratePluginEntryPoints(sout, logoFiles);
  GeneratePluginDataSymbols(sout, logoFiles, InputXML);
  GenerateSplitString(sout);
  GenerateSplitFilenames(sout);
  GeneratePluginProcedures(sout, logoFiles);
  GenerateLOGO(sout, logoFiles);
  GenerateXML(sout);

#ifdef GenerateCLP_USE_JSONCPP
  ModuleParameterGroup serializationGroup;
  serializationGroup.SetLabel("Parameter Serialization");

  ModuleParameter parametersSerialize;
  parametersSerialize.SetTag("file");
  parametersSerialize.SetCPPType("std::string");
  parametersSerialize.SetName("parametersSerialize");
  parametersSerialize.SetLongFlag("serialize");
  parametersSerialize.SetDescription("Store the module's parameters to a file.");
  parametersSerialize.SetValue("");
  parametersSerialize.SetChannel("output");
  serializationGroup.AddParameter(parametersSerialize);

  module.AddParameterGroup(serializationGroup);
#endif // GenerateCLP_USE_JSONCPP
  GenerateDeclare(sout, module);
  GenerateTCLAP(sout, module);
  GenerateTCLAPAssignment(sout, module, true);
  GenerateEchoArgs(sout, module);
  GenerateProcessInformationAddressDecoding(sout);
  GeneratePost(sout);
#ifndef GENERATECLP_USE_MD5
  sout.close();
#endif

#ifdef GENERATECLP_USE_MD5

  bool writeOutputBuffer = true;
  bool outputFileExists = itksys::SystemTools::FileExists(OutputCxx.c_str());
  if (outputFileExists)
    {
    // Compute MD5 of output buffer
    std::string outputBufferMD5 = itksys::SystemTools::LowerCase(ComputeStringMD5(sout.str().c_str()));

    // Compute MD5 of existing file
    char computedExistingOutputFileMD5[32];
    bool success = ComputeFileMD5(OutputCxx.c_str(), computedExistingOutputFileMD5);

    if (success)
      {
      std::string existingOutputFileMD5 =
          itksys::SystemTools::LowerCase(std::string(computedExistingOutputFileMD5, 32));

      if (outputBufferMD5.compare(existingOutputFileMD5) == 0)
        {
        std::cout << "GenerateCLP: File "
                  << itksys::SystemTools::GetFilenameName(OutputCxx) << " up-to-date." << std::endl;
        writeOutputBuffer = false;
        }
      }
    else
      {
      std::cerr << argv[0] << ": Failed to compute MD5 of file:" << OutputCxx << std::endl;
      }
    }

  if (writeOutputBuffer)
    {
    std::cout << parametersGroupsMsg.str();
    std::cout.flush();

    std::ofstream sOutputFile(OutputCxx.c_str(),std::ios::out);
    if (sOutputFile.fail())
      {
      std::cerr << argv[0] << ": Cannot open " << OutputCxx << " for output" << std::endl;
      perror(argv[0]);
      return EXIT_FAILURE;
      }
    sOutputFile << sout.str();
    sOutputFile.close();
    }
#endif

  return (EXIT_SUCCESS);
}

namespace
{

class MacroWriter
{
public:
  MacroWriter(std::ostream &sout, const char * macroName)
    : Stream(sout), Indent("    "), EOL(" \\\n")
    {
    sout << "#define " << macroName;
    }

  ~MacroWriter()
    {
    this->Stream << std::endl;
    }

  template <typename T> MacroWriter& operator+(const T& data)
    {
    this->Stream << data;
    return *this;
    }

  template <typename T> MacroWriter& operator|(const T& data)
    {
    this->Stream << this->EOL << this->Indent << data;
    return *this;
    }

private:
  std::ostream &Stream;
  const std::string Indent;
  const std::string EOL;
};

void GeneratePre(std::ostream &sout, ModuleDescription &, int argc, char *argv[])
{
  sout << "// This file was automatically generated by:" << std::endl;
  sout << "// ";
  for (int i = 0; i < argc; i++)
    {
    sout << " " << argv[i];
    }
  sout << std::endl;
  sout << "//" << std::endl;
  sout << "#include <cstdio>" << std::endl;
  sout << "#include <cstdlib>" << std::endl;
  sout << "#include <iostream>" << std::endl;
  sout << "#include <string>" << std::endl;
  sout << "#include <vector>" << std::endl;
  sout << "#include <map>" << std::endl;
  sout << std::endl;
  sout << "#include <sstream>" << std::endl;
#ifdef GenerateCLP_USE_JSONCPP
  sout << "#include <fstream>\n";
#endif // GenerateCLP_USE_JSONCPP
#ifdef GenerateCLP_USE_SERIALIZER
  sout << "#include <stdexcept>\n";
#endif // GenerateCLP_USE_SERIALIZER
  sout << std::endl;
  sout << "#include \"tclap/CmdLine.h\"" << std::endl;
#ifdef GenerateCLP_USE_JSONCPP
  sout << "#include \"JsonSerializationUtilities.h\"\n";
#endif // GenerateCLP_USE_JSONCPP
  sout << "#include \"ModuleProcessInformation.h\"" << std::endl;
  sout << std::endl;
}

#ifdef GenerateCLP_USE_SERIALIZER
void ReplaceCharacter(std::string& line, char character, const std::string replacement)
{
  size_t position = line.find( character );
    while( position != std::string::npos )
      {
      line = line.replace( position, 1, replacement );
      position = line.find( character, position + 2 );
      }
}

void GenerateJson( std::ostream & sout,
  itk::SEMModuleDescriptionSerializer::Pointer & serializer )
{
  sout << "\nextern \"C\" {\n";
  sout << "Module_EXPORT char JSONModuleDescription[] =" << std::endl;
  std::ostringstream oStream;
  itk::JsonCppArchiver::Pointer archiver =
    dynamic_cast< itk::JsonCppArchiver * >( serializer->GetArchiver() );
  archiver->WriteToStdStream( oStream );
  std::istringstream iStream( oStream.str() );
  while( !iStream.eof() )
    {
    std::string line;
    std::getline( iStream, line );
    ReplaceCharacter(line, '\\', "\\\\");
    ReplaceCharacter(line, '"', "\\\"");
    sout << "\"" << line << "\\n\"\n";
    }
  sout << ";\n";
  sout << "}" << std::endl;
}
#endif // GenerateCLP_USE_SERIALIZER

#ifdef GenerateCLP_USE_JSONCPP
void GenerateDeSerialization( std::ostream & sout,
  const ModuleDescription & module )
{
  MacroWriter w(sout, "GENERATE_DESERIALIZATION");

  w | "if (argc >= 2)"
    | "  {"
    | "  for (int arg_idx = 0; arg_idx < argc; ++arg_idx)"
    | "    {"
    | "    if (strcmp(argv[arg_idx],\"--deserialize\") == 0)"
    | "      {"
    | "      std::ifstream fStream( argv[arg_idx+1] );"
    | "      if( !fStream.is_open() )"
    | "        {"
    | "        std::cerr << \"Could not open file: \" << argv[arg_idx+1] << \" for writing.\" << std::endl;"
    | "        return EXIT_FAILURE;"
    | "        }"
    | "      Json::Value root;"
    | "      Json::Reader reader;"
    | "      reader.parse( fStream, root );"
    | "      const Json::Value & parameters = root[\"Parameters\"];";

  typedef std::vector<ModuleParameterGroup> ModuleParameterGroupsType;
  const ModuleParameterGroupsType & parameterGroups = module.GetParameterGroups();
  for( ModuleParameterGroupsType::const_iterator groupIt = parameterGroups.begin();
       groupIt != parameterGroups.end();
       ++groupIt )
    {
    const std::string & groupLabel = groupIt->GetLabel();
    typedef std::vector< ModuleParameter > ModuleParametersType;
    const ModuleParametersType & parameters = groupIt->GetParameters();
    for( ModuleParametersType::const_iterator paramIt = parameters.begin();
         paramIt != parameters.end();
         ++paramIt )
      {
      // Use LongFlag if present, Name otherwise
      const std::string & parameterName = paramIt->GetName();

      if (paramIt->IsReturnParameter())
        {
        continue;
        }

      w | "      if (!parameters[\"" + groupLabel + "\"][\"" + parameterName + "\"].isNull())"
        | "        {";
      std::string flag = "";
      if (!paramIt->GetLongFlag().empty())
        {
        flag = "--" + paramIt->GetLongFlag();
        }
      else if (!paramIt->GetFlag().empty())
        {
        flag = "-" + paramIt->GetFlag();
        }

      if (flag != "")
        {
        if (paramIt->GetCPPType() == "bool")
          {
          w | "        if (parameters[\"" + groupLabel + "\"][\"" + parameterName + "\"].asBool())"
            | "          {"
            | "          deserializedVectorFlaggedArgs.push_back(\"" + flag + "\");"
            | "          }";
          }
        else if (IsVectorOfVectors(*paramIt) || paramIt->GetMultiple() == "true")
          {
          w | "        Json::Value param = parameters[\"" + groupLabel + "\"][\"" + parameterName + "\"];"
            | "        if(param.size() > 0)"
            | "          {"
            | "          std::vector<std::string> values;"
            | "          for (unsigned int i = 0; i < param.size(); ++i)"
            | "            {";
          if (IsVectorOfVectors(*paramIt))
            {
            w | "            std::string value = \"\";"
              | "            for (unsigned int j = 0; j < param[i].size(); ++j)"
              | "              {"
              | "              value += param[i][j].asString();"
              | "              if (j < param[i].size() - 1)"
              | "                {"
              | "                value += \", \";"
              | "                }"
              | "              }"
              | "            values.push_back(value);";
            }
          else
            {
            w | "            values.push_back(param[i].asString());";
            }

          w | "            }"
            | "          if (!values.empty())"
            | "            {"
            | "            deserializedMultipleArgsMap[\"" + flag + "\"] = values;"
            | "            }"
            | "          }";
          }
        else
          {
          w | "        Json::Value param = parameters[\"" + groupLabel + "\"][\"" + parameterName + "\"];"
            | "        if (!param.isArray())"
            | "          {"
            | "          deserializedVectorFlaggedArgs.push_back(\"" + flag + "\");"
            | "          deserializedVectorFlaggedArgs.push_back(param.asString());"
            | "          }"
            | "        else if(param.size() > 0)"
            | "          {"
            | "          std::string value = \"\";"
            | "          for (unsigned int i = 0; i < param.size(); ++i)"
            | "            {"
            | "            value += param[i].asString();"
            | "            if (i < param.size() - 1)"
            | "              {"
            | "              value += \", \";"
            | "              }"
            | "            }"
            | "          deserializedVectorFlaggedArgs.push_back(\"" + flag + "\");"
            | "          deserializedVectorFlaggedArgs.push_back(value);"
            | "          }";
          }
        }
      else
        {
        if (paramIt->GetMultiple() == "true")
          {
          w | "        Json::Value param = parameters[\"" + groupLabel + "\"][\"" + parameterName + "\"]"
            | "        if(param.size() > 0)"
            | "          {"
            | "          for (unsigned int i = 0; i < param.size(); ++i)"
            | "            {"
            | "            deserializedVectorPositionalArgs.push_back(param[i].asString());"
            | "            }"
            | "          }";
          }
        else
          {
          w | "        deserializedVectorPositionalArgs.push_back(parameters[\"" + groupLabel + "\"][\"" + parameterName + "\"].asString());";
          }
        }
      w | "        }";
      }
    }
  w | "      }";
  w | "    }";
  w | "  }";
}

void GenerateSerialization( std::ostream & sout,
  const ModuleDescription & module )
{
  MacroWriter w(sout, "GENERATE_SERIALIZATION");
  w | "if( parametersSerializeArg.isSet() )"
    | "  {"
    | "  Json::Value root;"
    | "  Json::Value & parameters = root[\"Parameters\"];";
  typedef std::vector<ModuleParameterGroup> ModuleParameterGroupsType;
  const ModuleParameterGroupsType & parameterGroups = module.GetParameterGroups();
  for( ModuleParameterGroupsType::const_iterator groupIt = parameterGroups.begin();
       groupIt != parameterGroups.end();
       ++groupIt )
    {
    const std::string & groupLabel = groupIt->GetLabel();
    w | "    {"
      | "    Json::Value & parameterGroup = parameters[\""
         + groupLabel + "\"];";

    typedef std::vector< ModuleParameter > ModuleParametersType;
    const ModuleParametersType & parameters = groupIt->GetParameters();
    for( ModuleParametersType::const_iterator paramIt = parameters.begin();
         paramIt != parameters.end();
         ++paramIt )
      {
      const std::string & parameterName = paramIt->GetName();
      w | "    parameterGroup[\"" + parameterName
           + "\"] = JsonSerialize( " + parameterName + " );";
      }
      w | "    }";
    }
  w | "  std::ofstream fStream( parametersSerializeArg.getValue().c_str() );"
    | "  if( !fStream.is_open() )"
    | "    {"
    | "    std::cerr << \"Could not open file: \" << parametersSerializeArg.getValue() << \" for writing.\" << std::endl;"
    | "    return EXIT_FAILURE;"
    | "    }"
    | "  Json::StyledStreamWriter writer;"
    | "  writer.write( fStream, root );"
    | "  fStream.close();"
    | "  }";
}
#endif // GenerateCLP_USE_JSONCPP

void GenerateSplitString(std::ostream &sout)
{
  sout << "void" << std::endl;
  sout << "splitString (const std::string &text," << std::endl;
  sout << "             const std::string &separators," << std::endl;
  sout << "             std::vector<std::string> &words)" << std::endl;
  sout << "{" << std::endl;
  sout << "  const std::string::size_type n = text.length();" << std::endl;
  sout << "  std::string::size_type start = text.find_first_not_of(separators);" << std::endl;
  sout << "  while (start < n)" << std::endl;
  sout << "    {" << std::endl;
  sout << "    std::string::size_type stop = text.find_first_of(separators, start);" << std::endl;
  sout << "    if (stop > n) stop = n;" << std::endl;
  sout << "    words.push_back(text.substr(start, stop - start));" << std::endl;
  sout << "    start = text.find_first_not_of(separators, stop+1);" << std::endl;
  sout << "    }" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;
}

void GenerateSplitFilenames(std::ostream &sout)
{
  sout << "void" << std::endl;
  sout << "splitFilenames (const std::string &text," << std::endl;
  sout << "                std::vector<std::string> &words)" << std::endl;
  sout << "{" << std::endl;
  sout << "  const std::string::size_type n = text.length();" << std::endl;
  sout << "  bool quoted;" << std::endl;
  sout << "  std::string comma(\",\");" << std::endl;
  sout << "  std::string quote(\"\\\"\");" << std::endl;
  sout << "  std::string::size_type start = text.find_first_not_of(comma);" << std::endl;
  sout << "  while (start < n)" << std::endl;
  sout << "    {" << std::endl;
  sout << "    quoted = false;" << std::endl;
  sout << "    std::string::size_type startq = text.find_first_of(quote, start);" << std::endl;
  sout << "    std::string::size_type stopq = text.find_first_of(quote, startq+1);" << std::endl;
  sout << "    std::string::size_type stop = text.find_first_of(comma, start);" << std::endl;
  sout << "    if (stop > n) stop = n;" << std::endl;
  sout << "    if (startq != std::string::npos && stopq != std::string::npos)"
       << std::endl;
  sout << "      {" << std::endl;
  sout << "      while (startq < stop && stop < stopq && stop != n)" << std::endl;
  sout << "         {" << std::endl;
  sout << "         quoted = true;" << std::endl;
  sout << "         stop = text.find_first_of(comma, stop+1);" << std::endl;
  sout << "         if (stop > n) stop = n;" << std::endl;
  sout << "         }" << std::endl;
  sout << "      }" << std::endl;
  sout << "    if (!quoted)" << std::endl;
  sout << "      {" << std::endl;
  sout << "      words.push_back(text.substr(start, stop - start));" << std::endl;
  sout << "      }" << std::endl;
  sout << "    else" << std::endl;
  sout << "      {" << std::endl;
  sout << "      words.push_back(text.substr(start+1, stop - start-2));" << std::endl;
  sout << "      }" << std::endl;
  sout << "    start = text.find_first_not_of(comma, stop+1);" << std::endl;
  sout << "    }" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;
}

void GenerateExports(std::ostream &sout)
{
  sout << "#ifdef _WIN32" << std::endl;
  sout << "#define Module_EXPORT __declspec(dllexport)" << std::endl;
  sout << "#elif defined(MODULE_HIDDEN_VISIBILITY)" << std::endl;
  sout << "#define Module_EXPORT __attribute__((visibility(\"default\")))" << std::endl;
  sout << "#else" << std::endl;
  sout << "#define Module_EXPORT" << std::endl;
  sout << "#endif" << std::endl;
  sout << std::endl;
}

void GeneratePluginDataSymbols(std::ostream &sout, std::vector<std::string>& logos, std::string XMLFile)
{
  sout << "extern \"C\" {" << std::endl;
  sout << "Module_EXPORT char XMLModuleDescription[] = " << std::endl;

  std::string line;
  std::ifstream fin(XMLFile.c_str(),std::ios::in);
  while (!fin.eof())
    {
    std::getline( fin, line );

    // replace quotes with escaped quotes
    std::string cleanLine;
    for (std::string::size_type j = 0; j < line.length(); j++)
      {
      if (line[j] == '\"')
        {
        cleanLine.append("\\\"");
        }
      else
        {
        cleanLine.append(1,line[j]);
        }
      }
    sout << "\"" << cleanLine << "\\n\"" << std::endl;
    }
  sout << ";" << std::endl << std::endl;

  fin.close();

  if (logos.size() == 1)
    {
    std::string logo = logos[0];
    std::string fileName = itksys::SystemTools::GetFilenameWithoutExtension (logo);

    sout << "#define static Module_EXPORT" << std::endl;
    sout << "#define const" << std::endl;
    sout << "#define image_" << fileName << "_width ModuleLogoWidth"
         << std::endl;
    sout << "#define image_" << fileName << "_height ModuleLogoHeight"
         << std::endl;
    sout << "#define image_" << fileName << "_pixel_size ModuleLogoPixelSize"
         << std::endl;
    sout << "#define image_" << fileName << "_length ModuleLogoLength"
         << std::endl;
    sout << "#define image_" << fileName << " ModuleLogoImage"
         << std::endl;
    sout << "#include \"" << logos[0] << "\"" << std::endl;
    sout << "#undef static" << std::endl;
    sout << "#undef const" << std::endl;
    sout << "#undef image_" << fileName << "_width" << std::endl;
    sout << "#undef image_" << fileName << "_height" << std::endl;
    sout << "#undef image_" << fileName << "_pixel_size" << std::endl;
    sout << "#undef image_" << fileName << "_length" << std::endl;
    sout << "#undef image_" << fileName << std::endl;
    }
  sout << "}" << std::endl;
  sout << std::endl;
}

void GeneratePluginEntryPoints(std::ostream &sout, std::vector<std::string> &logos)
{
  sout << "#if defined(main) && !defined(REGISTER_TEST)" << std::endl;
  sout << "// If main defined as a preprocessor symbol, redefine it to the expected entry point." << std::endl;
  sout << "#undef main" << std::endl;
  sout << "#define main ModuleEntryPoint" << std::endl;
  sout << std::endl;

  sout << "extern \"C\" {" << std::endl;
  sout << "  Module_EXPORT char *GetXMLModuleDescription();" << std::endl;
  sout << "  Module_EXPORT int ModuleEntryPoint(int, char*[]);" << std::endl;
  if (logos.size() == 1)
    {
    sout << "  Module_EXPORT unsigned char *GetModuleLogo(int *width, int *height, int *pixel_size, unsigned long *bufferLength);" << std::endl;
    }
  sout << "}" << std::endl;
  sout << "#endif" << std::endl;
  sout << std::endl;
}

void GeneratePluginProcedures(std::ostream &sout, std::vector<std::string> &logos)
{
  if (logos.size() == 1)
    {
    sout << "unsigned char *GetModuleLogo(int *width," << std::endl;
    sout << "                             int *height," << std::endl;
    sout << "                             int *pixel_size," << std::endl;
    sout << "                             unsigned long *length)" << std::endl;
    sout << "{" << std::endl;

    sout << "  *width = ModuleLogoWidth;" << std::endl;
    sout << "  *height = ModuleLogoHeight;" << std::endl;
    sout << "  *pixel_size = ModuleLogoPixelSize;" << std::endl;
    sout << "  *length = ModuleLogoLength;" << std::endl;
    sout << "  return const_cast<unsigned char *>(ModuleLogoImage);" << std::endl;
    sout << "}" << std::endl;
    sout << std::endl;
    }

  sout << "char *GetXMLModuleDescription()" << std::endl;
  sout << "{" << std::endl;
  sout << "   return XMLModuleDescription;" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;
}

void GenerateLOGO(std::ostream &sout, std::vector<std::string> &logos)
{
  if (logos.size() > 0)
    {
    std::string EOL(" \\");

    sout << "#define GENERATE_LOGO \\" << std::endl;
    // Generate special section to produce logo description
    sout << "  if (argc >= 2 && (strcmp(argv[1],\"--logo\") == 0))" << EOL << std::endl;
    sout << "    {" << EOL << std::endl;
    sout << "    int width, height, pixel_size;    " << EOL << std::endl;
    sout << "    unsigned long length; " << EOL << std::endl;
    sout << "    unsigned char *logo = GetModuleLogo(&width, &height, &pixel_size, &length); " << EOL << std::endl;
    sout << "    std::cout << \"LOGO\" << std::endl; " << EOL << std::endl;
    sout << "    std::cout << width << std::endl; " << EOL << std::endl;
    sout << "    std::cout << height << std::endl; " << EOL << std::endl;
    sout << "    std::cout << pixel_size << std::endl; " << EOL << std::endl;
    sout << "    std::cout << length << std::endl; " << EOL << std::endl;
    sout << "    std::cout << logo << std::endl; " << EOL << std::endl;
    sout << "    return EXIT_SUCCESS; " << EOL << std::endl;
    sout << "    }" << std::endl;
    }
  else
    {
    sout << "#define GENERATE_LOGO" << std::endl;
    }
}

void GenerateXML(std::ostream &sout)
{
  std::string EOL(" \\");

  sout << "#define GENERATE_XML \\" << std::endl;
  // Generate special section to produce xml description
  sout << "  if (argc >= 2 && (strcmp(argv[1],\"--xml\") == 0))" << EOL << std::endl;
  sout << "    {" << EOL << std::endl;
  sout << "    std::cout << GetXMLModuleDescription();" << EOL << std::endl;
  sout << "    return EXIT_SUCCESS;" << EOL << std::endl;
  sout << "    }" << std::endl;
}

void GenerateEchoArgs(std::ostream &sout, ModuleDescription &module)
{
  std::string EOL(" \\");
  sout << "#define GENERATE_ECHOARGS \\" << std::endl;

  sout << "if (echoSwitch)" << EOL << std::endl;
  sout << "{" << EOL << std::endl;
  sout << "std::cout << \"Command Line Arguments\" << std::endl;" << EOL << std::endl;
  std::vector<ModuleParameterGroup>::const_iterator git;
  std::vector<ModuleParameter>::const_iterator pit;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (NeedsTemp(*pit) && pit->GetMultiple() != "true")
        {
        sout << "std::cout << "
             << "\"    "
             << pit->GetName()
             << ": \";"
             << EOL << std::endl;
        sout << "for (unsigned int _i =0; _i < "
             << pit->GetName()
             << ".size(); _i++)"
             << EOL << std::endl;
        sout << "{"
             << EOL << std::endl;
        sout << "std::cout << "
             << pit->GetName()
             << "[_i]"
             << " << \", \";"
             << EOL << std::endl;
        sout << "}"
             << EOL << std::endl;
        sout << "std::cout <<std::endl;"
             << EOL << std::endl;

        }
      else if (NeedsTemp(*pit) && pit->GetMultiple() == "true")
        {
        sout << "for (unsigned int _i= 0; _i < " << pit->GetName() << "Temp.size(); _i++)" << EOL << std::endl;
        sout << "{" << EOL << std::endl;
        sout << "std::cout << \"" << pit->GetName() << "[\" << _i << \"]: \";" << EOL << std::endl;
        sout << "std::vector<std::string> words;" << EOL << std::endl;
        sout << "words.clear();" << EOL << std::endl;
        if ((*pit).GetTag() == "file")
          {
          sout << "splitFilenames(" << pit->GetName() << "Temp[_i], words);" << EOL << std::endl;
          }
        else
          {
          sout << "      std::string sep(\",\");" << EOL << std::endl;
          sout << "splitString(" << pit->GetName() << "Temp[_i], sep, words);" << EOL << std::endl;
          }
        sout << "for (unsigned int _j= 0; _j < words.size(); _j++)" << EOL << std::endl;
        sout << "{" << EOL << std::endl;
        sout << "std::cout <<  words[_j] << \" \";" << EOL << std::endl;
        sout << "}" << EOL << std::endl;
        sout << "std::cout << std::endl;" << EOL << std::endl;
        sout << "}" << EOL << std::endl;
        }
      else if (pit->GetMultiple() == "true")
        {
        sout << "std::cout << "
             << "\"    "
             << pit->GetName()
             << ": \";"
             << EOL << std::endl;
        sout << "for (unsigned int _i =0; _i < "
             << pit->GetName()
             << ".size(); _i++)"
             << EOL << std::endl;
        sout << "{"
             << EOL << std::endl;
        sout << "std::cout << "
             << pit->GetName()
             << "[_i]"
             << " << \", \";"
             << EOL << std::endl;
        sout << "}"
             << EOL << std::endl;
        sout << "std::cout <<std::endl;"
             << EOL << std::endl;
        }
      else
        {
        sout << "std::cout << "
             << "\"    "
             << pit->GetName()
             << ": \" << "
             << pit->GetName()
             << " << std::endl;"
             << EOL << std::endl;
        }
      }
    }
  sout << "}" << std::endl;
}

void GenerateDeclare(std::ostream &sout, ModuleDescription &module)
{
  std::string EOL(" \\");
  sout << "#define GENERATE_DECLARE \\" << std::endl;

  ModuleParameterGroup autoParameters;

  // Add a switch argument to echo command line arguments
  ModuleParameter echoSwitch;
  echoSwitch.SetTag("boolean");
  echoSwitch.SetCPPType("bool");
  echoSwitch.SetName("echoSwitch");
  echoSwitch.SetLongFlag("echo");
  echoSwitch.SetDescription("Echo the command line arguments");
  echoSwitch.SetValue("false");
  autoParameters.AddParameter(echoSwitch);

  // Add a switch argument to produce xml output
  ModuleParameter xmlSwitch;
  xmlSwitch.SetTag("boolean");
  xmlSwitch.SetCPPType("bool");
  xmlSwitch.SetName("xmlSwitch");
  xmlSwitch.SetLongFlag("xml");
  xmlSwitch.SetDescription("Produce xml description of command line arguments");
  xmlSwitch.SetValue("false");
  autoParameters.AddParameter(xmlSwitch);

  // Add an argument to accept an address for storing process
  // information
  ModuleParameter processInformationAddressArg;
  processInformationAddressArg.SetTag("string");
  processInformationAddressArg.SetCPPType("std::string");
  processInformationAddressArg.SetName("processInformationAddressString");
  processInformationAddressArg.SetLongFlag("processinformationaddress");
  processInformationAddressArg.SetDescription("Address of a structure to store process information (progress, abort, etc.).");
  processInformationAddressArg.SetValue("0");
  autoParameters.AddParameter(processInformationAddressArg);

  // Add an argument to accept a filename for simple return types
  ModuleParameter returnParameterArg;
  returnParameterArg.SetTag("file");
  returnParameterArg.SetCPPType("std::string");
  returnParameterArg.SetName("returnParameterFile");
  returnParameterArg.SetLongFlag("returnparameterfile");
  returnParameterArg.SetDescription("Filename in which to write simple return parameters (int, float, int-vector, etc.) as opposed to bulk return parameters (image, geometry, transform, measurement, table).");
  returnParameterArg.SetValue("");
  autoParameters.AddParameter(returnParameterArg);


  // Add the parameter group to the module
  module.AddParameterGroup(autoParameters);

#ifdef GenerateCLP_USE_JSONCPP
  // Add the deserialize argument to the serialization group
  ModuleParameter parametersDeSerialize;
  parametersDeSerialize.SetTag("file");
  parametersDeSerialize.SetCPPType("std::string");
  parametersDeSerialize.SetName("parametersDeSerialize");
  parametersDeSerialize.SetLongFlag("deserialize");
  parametersDeSerialize.SetDescription("Restore the module's parameters that were previously archived.");
  parametersDeSerialize.SetValue("");
  parametersDeSerialize.SetChannel("input");

  // Find the serialization group and add to it.
  std::vector<ModuleParameterGroup>::iterator groupsIt;
  for (groupsIt = module.GetParameterGroups().begin();
    groupsIt != module.GetParameterGroups().end(); ++groupsIt)
    {
    if (groupsIt->GetLabel() == "Parameter Serialization")
      {
      groupsIt->AddParameter(parametersDeSerialize);
      break;
      }
    }
#endif

  sout << "    /* These two vectors are used to store the JSON deserialized value */" << EOL << std::endl;
  sout << "    /* that are then compiled with the command line. */" << EOL << std::endl;
  sout << "    std::vector< std::string > deserializedVectorFlaggedArgs;" << EOL << std::endl;
  sout << "    std::vector< std::string > deserializedVectorPositionalArgs;" << EOL << std::endl;
  sout << "    /* This map is used to store the JSON deserialized value of multiple args*/" << EOL << std::endl;
  sout << "    /* where the key is the argument flag and the value the values of each arg. */" << EOL << std::endl;
  sout << "    std::map< std::string, std::vector<std::string> > deserializedMultipleArgsMap;" << EOL << std::endl;
  sout << EOL << std::endl;
  sout << "    /* This vector is used to look up if a flag requires an argument after it. */" << EOL << std::endl;
  sout << "    /* This is used to differentiate between: */" << EOL << std::endl;
  sout << "    /* ./myExec --boolFlag /my/first/arg */" << EOL << std::endl;
  sout << "    /* ./myExec --flag flagArg */" << EOL << std::endl;
  sout << "    std::vector< std::string > nonbooleanFlags;" << EOL << std::endl;
  std::vector<ModuleParameterGroup>::const_iterator git;
  std::vector<ModuleParameter>::const_iterator pit;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (pit->GetCPPType() != "bool")
        {
        if(!pit->GetFlag().empty())
          {
          sout << "    nonbooleanFlags.push_back(\"-" << pit->GetFlag() <<"\");" << EOL << std::endl;
          }
        if (!pit->GetLongFlag().empty())
          {
          sout << "    nonbooleanFlags.push_back(\"--" << pit->GetLongFlag() <<"\");" << EOL << std::endl;
          }
        }
      }
    }
  sout << "    /* This map use is twofold: */" << EOL << std::endl;
  sout << "    /*  - to find whether a flag is multiple */" << EOL << std::endl;
  sout << "    /*  - to know if we need to reset the multiple arg value because it was */" << EOL << std::endl;
  sout << "    /*    in the JSON and it's also in the command line. */" << EOL << std::endl;
  sout << "    std::map<std::string, bool> multipleFlags;" << EOL << std::endl;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (pit->GetMultiple() == "true" && !pit->GetFlag().empty())
        {
        sout << "    multipleFlags[\"-" << pit->GetFlag() << "\"] = false;" << EOL << std::endl;
        }
      if (pit->GetMultiple() == "true" && !pit->GetLongFlag().empty())
        {
        sout << "    multipleFlags[\"--" << pit->GetLongFlag() << "\"] = false;" << EOL << std::endl;
        }
      }
    }

  // First pass generates argument declarations
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (!IsEnumeration(*pit))
        {
        sout << "    ";
        if (pit->GetMultiple() == "true")
          {
          sout << "std::vector<std::string>";
          }
        else if (NeedsQuotes(*pit))
          {
          sout << "std::string";
          }
        else
          {
          sout << pit->GetCPPType();
          }
        sout << " ";
        sout << pit->GetName();
        if (NeedsTemp(*pit))
          {
          sout << "Temp";
          }

        const std::string cppType = pit->GetCPPType();
        if (!HasValue(*pit) &&
            cppType != "bool")
          {
          // Initialized to avoid compiler warnings.
          if (cppType.compare("int") == 0)
            {
            sout << " = 0";
            }
          else if (cppType.compare("float") == 0)
            {
            sout << " = 0.0f";
            }
          else if (cppType.compare("double") == 0)
            {
            sout << " = 0.0";
            }
          sout << ";"
               << EOL << std::endl;
          }
        else
          {
          std::string ValueString = pit->GetValue();
          if ((*pit).GetCPPType() == "bool")
            {
            ValueString = "false";
            }
          replaceSubWithSub(ValueString, "\"", "\\\"");
          sout << " = ";
          if (NeedsQuotes(*pit))
            {
            sout << "\"";
            }
          sout << ValueString;
          if (NeedsQuotes(*pit))
            {
            sout << "\"";
            }
          sout << ";";
          sout << EOL << std::endl;
          }
        if (NeedsTemp(*pit))
          {
          sout << "    "
               << pit->GetCPPType()
               << " "
               << pit->GetName()
               << ";"
               << EOL << std::endl;
          }
        }
      else // IsEnumeration(*pit)
        {
        sout << "    "
             << pit->GetCPPType()
             << " ";
        sout << pit->GetName();
        if (!HasValue(*pit))
          {
          sout << ";"
               << EOL << std::endl;
          }
        else
          {
          sout << " = ";
          if (NeedsQuotes(*pit))
            {
            sout << "\"";
            }
          std::string ValueString = pit->GetValue();
          if ((*pit).GetCPPType() == "bool")
            {
            ValueString = "false";
            }
          replaceSubWithSub(ValueString, "\"", "\\\"");
          sout << ValueString;
          if (NeedsQuotes(*pit))
            {
            sout << "\"";
            }
          sout << ";"
               << EOL << std::endl;
          }
        sout << "    "
             << "std::vector<" << pit->GetCPPType() << "> "
             <<  pit->GetName() << "Allowed;"
             << EOL << std::endl;
        for (unsigned int e = 0; e < (pit->GetElements()).size(); e++)
          {
          sout << "    "
               << pit->GetName() << "Allowed.push_back(";
          if (NeedsQuotes(*pit))
            {
            sout << "\"";
            }
          std::string element = pit->GetElements()[e];
          replaceSubWithSub(element, "\"", "\\\"");

          sout << element;
          if (NeedsQuotes(*pit))
            {
            sout << "\"";
            }
          sout << "); "
               << EOL << std::endl;
          }
          sout << "    "
            "TCLAP::ValuesConstraint<" << pit->GetCPPType() << "> "
               << pit->GetName() << "AllowedVals ("
               << pit->GetName() << "Allowed); "
               << EOL << std::endl;
        }
      }
    }

  sout << std::endl;
}

void GenerateTCLAP(std::ostream &sout, ModuleDescription &module)
{
  GenerateTCLAPParse(sout, module);
  GenerateTCLAPAssignment(sout, module);
  sout << "#define GENERATE_TCLAP GENERATE_TCLAP_PARSE;GENERATE_TCLAP_ASSIGNMENT" << std::endl;
}

void GenerateTCLAPParse(std::ostream &sout, ModuleDescription &module)
{
  std::string EOL(" \\");
  sout << "#define GENERATE_TCLAP_PARSE \\" << std::endl;

  sout << "    std::string fullDescription(\"Description: \");" << EOL << std::endl;
  sout << "    fullDescription += \"" << module.GetDescription() << "\";" << EOL << std::endl;
  sout << "    if (!std::string(\"" << module.GetContributor() << "\").empty())" << EOL << std::endl;
  sout << "      {" << EOL << std::endl;
  sout << "      fullDescription += \"\\nAuthor(s): " << module.GetContributor() << "\";" << EOL << std::endl;
  sout << "      }" << EOL << std::endl;
  sout << "    if (!std::string(\"" << module.GetAcknowledgements() << "\").empty())" << EOL << std::endl;
  sout << "      {" << EOL << std::endl;
  sout << "      fullDescription += \"\\nAcknowledgements: " << module.GetAcknowledgements() << "\";" << EOL << std::endl;
  sout << "      }" << EOL << std::endl;
  sout << "    TCLAP::CmdLine commandLine (fullDescription," << EOL << std::endl;
  sout << "       " << "' '," << EOL << std::endl;
  sout << "      " << "\"" << module.GetVersion() << "\"";
  sout << " );" << EOL << std::endl << EOL << std::endl;
  sout << "      std::ostringstream msg;" << EOL << std::endl;

  // Second pass generates TCLAP declarations
  std::vector<ModuleParameterGroup>::const_iterator git;
  std::vector<ModuleParameter>::const_iterator pit;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      // Don't generate TCLAP data structures for parameters that are
      // simple return types
      if ((*pit).IsReturnParameter())
        {
        continue;
        }

      sout << "    msg.str(\"\");";
      sout << "msg << "
           << "\""
           << pit->GetDescription();
      if (pit->GetValue().empty())
        {
        sout << "\";";
        }
      else
        {
        sout << " (value: \" "
             << "<< ";
        if (HasValue(*pit))
          {
          sout << pit->GetName();
          }
        else
          {
          sout << "\"None\"";
          }

        if (NeedsTemp(*pit) && HasValue(*pit))
          {
          sout << "Temp";
          }
        sout << " << "
             << "\")"
             << "\";"
             << EOL << std::endl;
        }

      if (pit->GetCPPType() == "bool")
        {
        sout << "    TCLAP::SwitchArg "
             << pit->GetName()
             << "Arg" << "(\""
             << pit->GetFlag()
             << "\", \""
             << pit->GetLongFlag()
             << "\", msg.str(), "
             << "commandLine"
             << ", "
             << pit->GetName()
             << ");"
             << EOL << std::endl << EOL << std::endl;
        }
      else
        {
        if (pit->GetFlag().empty() && pit->GetLongFlag().empty())
          {
          if (pit->GetMultiple() == "true")
            {
            sout << "    TCLAP::UnlabeledMultiArg<";
            sout << pit->GetArgType();            }
          else
            {
            sout << "    TCLAP::UnlabeledValueArg<";
            sout << pit->GetCPPType();
            }
          sout << "> "
               << pit->GetName()
               << "Arg" << "(\""
               << pit->GetName()
               << "\", msg.str(), ";
          if (pit->GetMultiple() != "true")
            {
            sout << true
                 << ", "
                 << pit->GetName();
            }
          else
            {
            sout << false;
            }
          sout << ", "
               << "\""
               << pit->GetCPPType()
               << "\""
               << ", "
               << "commandLine);"
               << EOL << std::endl << EOL << std::endl;
          }
        else if (IsEnumeration(*pit))
          {
          sout << "    TCLAP::ValueArg<";
          sout << pit->GetCPPType();
          sout << "> "
               << pit->GetName()
               << "Arg" << "(\""
               << pit->GetFlag()
               << "\", \""
               << pit->GetLongFlag()
               << "\", msg.str(), "
               << false
               << ", "
               << pit->GetName();
          sout << ", "
               << "&" << pit->GetName() << "AllowedVals"
               << ", "
               << "commandLine);"
               << EOL << std::endl << EOL << std::endl;
          }
        else
          {
          if (pit->GetMultiple() == "true")
            {
            sout << "    TCLAP::MultiArg<";
            if (NeedsTemp(*pit))
              {
              sout << "std::string";
              }
            else
              {
              sout << pit->GetCPPType();
              }
            }
          else
            {
            sout << "    TCLAP::ValueArg<";
            if (NeedsTemp(*pit))
              {
              sout << "std::string";
              }
            else
              {
              sout << pit->GetCPPType();
              }
            }
          sout << " > "
               << pit->GetName()
               << "Arg" << "(\""
               << pit->GetFlag()
               << "\", \""
               << pit->GetLongFlag()
               << "\", msg.str(), "
               << false;
          if (pit->GetMultiple() != "true")
            {
            sout << ", "
                 << pit->GetName();
            if (NeedsTemp(*pit))
              {
              sout << "Temp";
              }
            }
          sout << ", "
               << "\""
               << pit->GetCPPType()
               << "\""
               << ", "
               << "commandLine);"
               << EOL << std::endl << EOL << std::endl;
          }
        }
      }
    }

  sout << "try" << EOL << std::endl;
  sout << "  {" << EOL << std::endl;

  // Remap any aliases in the flags or long flags
  sout << "    /* Build a map of flag aliases to the true flag */" << EOL << std::endl;
  sout << "    std::map<std::string,std::string> flagAliasMap;" << EOL << std::endl;
  sout << "    std::map<std::string,std::string> deprecatedFlagAliasMap;" << EOL << std::endl;
  sout << "    std::map<std::string,std::string> longFlagAliasMap;" << EOL << std::endl;
  sout << "    std::map<std::string,std::string> deprecatedLongFlagAliasMap;" << EOL << std::endl;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (pit->GetFlagAliases().size() != 0)
        {
        std::vector<std::string>::const_iterator ait;
        for (ait = pit->GetFlagAliases().begin();
             ait != pit->GetFlagAliases().end(); ++ait)
          {
          sout << "    flagAliasMap[\"-" << (*ait) << "\"] = \"-"
               << pit->GetFlag() << "\";" << EOL << std::endl;
          }
        }
      }
    }
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (pit->GetDeprecatedFlagAliases().size() != 0)
        {
        std::vector<std::string>::const_iterator ait;
        for (ait = pit->GetDeprecatedFlagAliases().begin();
             ait != pit->GetDeprecatedFlagAliases().end(); ++ait)
          {
          sout << "    deprecatedFlagAliasMap[\"-" << (*ait) << "\"] = \"-"
               << pit->GetFlag() << "\";" << EOL << std::endl;
          }
        }
      }
    }
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (pit->GetLongFlagAliases().size() != 0)
        {
        std::vector<std::string>::const_iterator ait;
        for (ait = pit->GetLongFlagAliases().begin();
             ait != pit->GetLongFlagAliases().end(); ++ait)
          {
          sout << "    longFlagAliasMap[\"--" << (*ait) << "\"] = \"--"
               << pit->GetLongFlag() << "\";" << EOL << std::endl;
          }
        }
      }
    }
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if (pit->GetDeprecatedLongFlagAliases().size() != 0)
        {
        std::vector<std::string>::const_iterator ait;
        for (ait = pit->GetDeprecatedLongFlagAliases().begin();
             ait != pit->GetDeprecatedLongFlagAliases().end(); ++ait)
          {
          sout << "    deprecatedLongFlagAliasMap[\"--" << (*ait) << "\"] = \"--"
               << pit->GetLongFlag() << "\";" << EOL << std::endl;
          }
        }
      }
    }
  sout << EOL << std::endl;

  sout << "    /* Go through argc and consolidate the JSON with the parameters from the command line */" << EOL << std::endl;
  sout << "    /* In case of conflict, take the command line. */" << EOL << std::endl;
  sout << "    std::map<std::string,std::string>::iterator ait;" << EOL << std::endl;
  sout << "    std::map<std::string,std::string>::iterator dait;" << EOL << std::endl;
  sout << "    std::vector<std::string>::iterator dvOptionnalArgsIt = deserializedVectorFlaggedArgs.begin();" << EOL << std::endl;
  sout << "    std::map<std::string, std::vector<std::string> >::iterator dvMultipleArgsIt;" << EOL << std::endl;
  sout << "    size_t noFlagCounter = 0;" << EOL << std::endl;
  sout << "    size_t ac = 1;" << EOL << std::endl;
  sout << EOL << std::endl;
  sout << "    while (ac < static_cast<size_t>(argc))" << EOL << std::endl;
  sout << "      {" << EOL << std::endl;
  sout << "      /* short flag case && long flag case */" << EOL << std::endl;
  // Simple short flags and long flags
  sout << "       if ((strlen(argv[ac]) == 2 && argv[ac][0]=='-')" << EOL << std::endl;
  sout << "           || (strlen(argv[ac]) > 2 && argv[ac][0]=='-' && argv[ac][1]=='-'))" << EOL << std::endl;
  sout << "         {" << EOL << std::endl;
  sout << "         std::string flag = argv[ac];" << EOL << std::endl;
  //   Flag re-mapping
  sout << "         /* Remap the flag if necessary */" << EOL << std::endl;
  sout << "         if (strlen(argv[ac]) == 2)" << EOL << std::endl;
  sout << "           {" << EOL << std::endl;
  //     short flag re-mapping
  sout << "           /* Short flag remapping */" << EOL << std::endl;
  sout << "           ait = flagAliasMap.find(flag);" << EOL << std::endl;
  sout << "           dait = deprecatedFlagAliasMap.find(flag);" << EOL << std::endl;
  sout << "           if (ait != flagAliasMap.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             /* remap the flag */" << EOL << std::endl;
  sout << "             flag = (*ait).second;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           else if (dait != deprecatedFlagAliasMap.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             std::cout << \"Flag \\\"\" << flag << \"\\\" is deprecated. Please use flag \\\"\" << (*dait).second << \"\\\" instead. \" << std::endl;" << EOL << std::endl;
  sout << "             /* remap the flag */" << EOL << std::endl;
  sout << "             flag = (*dait).second;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         else" << EOL << std::endl;
  sout << "           {" << EOL << std::endl;
  sout << "           /* Long flag remapping */" << EOL << std::endl;
  sout << "           ait = longFlagAliasMap.find(flag);" << EOL << std::endl;
  sout << "           dait = deprecatedLongFlagAliasMap.find(flag);" << EOL << std::endl;
  sout << "           if (ait != longFlagAliasMap.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             /* remap the flag */" << EOL << std::endl;
  sout << "             flag = (*ait).second;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           else if (dait != deprecatedLongFlagAliasMap.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             std::cout << \"Long flag \\\"\" << flag << \"\\\" is deprecated. Please use long flag \\\"\" << (*dait).second << \"\\\" instead. \" << std::endl;" << EOL << std::endl;
  sout << "             /* remap the flag */" << EOL << std::endl;
  sout << "             flag = (*dait).second;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         bool isMultiple = multipleFlags.find(flag) != multipleFlags.end();" << EOL << std::endl;
  sout << "         bool isBoolean = std::find(nonbooleanFlags.begin(), nonbooleanFlags.end(), flag) == nonbooleanFlags.end();" << EOL << std::endl;
  sout << "         dvOptionnalArgsIt = std::find(deserializedVectorFlaggedArgs.begin(), deserializedVectorFlaggedArgs.end(), flag);" << EOL << std::endl;
  sout << "         bool isPresentVFA = dvOptionnalArgsIt != deserializedVectorFlaggedArgs.end();" << EOL << std::endl;
  //   Compiling: 3 booleans => 8 cases
  //   But isBoolean && isMultiple is impossible, so a total of => 6 cases
  sout << "         if (isBoolean)" << EOL << std::endl;
  //     Boolean cases:
  sout << "           {" << EOL << std::endl;
  sout << "           /*Ignore if boolean and already present*/" << EOL << std::endl;
  sout << "           /*Otherwise add it*/" << EOL << std::endl;
  sout << "           if (!isPresentVFA)" << EOL << std::endl;
  //       Already present -> Do nothing
  //       Not present -> Add the flag
  sout << "             {" << EOL << std::endl;
  sout << "             deserializedVectorFlaggedArgs.push_back(flag);" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           ++ac;" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         else if (isMultiple)" << EOL << std::endl;
  sout << "           {" << EOL << std::endl;
  //     Multiple cases:
  sout << "           dvMultipleArgsIt = deserializedMultipleArgsMap.find(flag);" << EOL << std::endl;
  sout << "           bool isPresentMA = dvMultipleArgsIt != deserializedMultipleArgsMap.end();" << EOL << std::endl;
  sout << "           /*Ignore if boolean and already present*/" << EOL << std::endl;
  sout << "           /*Reset/Add the value if first deserialize or not present*/" << EOL << std::endl;
  sout << "           if (!isPresentMA || !multipleFlags[flag])" << EOL << std::endl;
  //       We need to reset/add the value if:
  //         - Not present
  //         - Was deserialized and it's now found in the command line
  sout << "             {" << EOL << std::endl;
  sout << "             deserializedMultipleArgsMap[flag] = std::vector<std::string>();" << EOL << std::endl;
  sout << "             multipleFlags[flag] = true;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  //       In any case, add the value and move on
  sout << "           ++ac;" << EOL << std::endl;
  sout << "           std::string value = \"\";" << EOL << std::endl;
  sout << "           if (ac < static_cast<size_t>(argc))" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             value = argv[ac];" << EOL << std::endl;
  sout << "             ++ac;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           deserializedMultipleArgsMap[flag].push_back(value);" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         else" << EOL << std::endl;
  //     Other cases:
  sout << "           {" << EOL << std::endl;
  sout << "           /*Add the flag and if needed*/" << EOL << std::endl;
  sout << "           if (!isPresentVFA)" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  //       Not here -> Add flag
  sout << "             deserializedVectorFlaggedArgs.push_back(flag);" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           ++ac;" << EOL << std::endl;
  sout << "           std::string value = \"\";" << EOL << std::endl;
  sout << "           if (ac < static_cast<size_t>(argc))" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             value = argv[ac];" << EOL << std::endl;
  sout << "             ++ac;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  //       In any case add/set the value and move on
  sout << "           if (!isPresentVFA)" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             deserializedVectorFlaggedArgs.push_back(value); "<< EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           else" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             *(++dvOptionnalArgsIt) = value;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         }" << EOL << std::endl;
  sout << "       /* short flag case where multiple flags are given at once */" << EOL << std::endl;
  sout << "       else if (strlen(argv[ac]) > 2 && argv[ac][0]=='-' && argv[ac][1]!='-')" << EOL << std::endl;
  sout << "         {" << EOL << std::endl;
  // Composed short flags
  sout << "         std::string rflag(argv[ac], 1, std::string::npos);" << EOL << std::endl;
  sout << "         for (std::string::size_type fi=0; fi < rflag.size(); ++fi)" << EOL << std::endl;
  sout << "           {" << EOL << std::endl;
  sout << "           std::string tf(rflag, fi, 1);" << EOL << std::endl;
  sout << "           std::string newFlag =\"-\";" << EOL << std::endl;
  sout << "           newFlag += tf;" << EOL << std::endl;
  //   Remapping
  sout << "           ait = flagAliasMap.find(newFlag);" << EOL << std::endl;
  sout << "           dait = deprecatedFlagAliasMap.find(newFlag);" << EOL << std::endl;
  sout << "           if (ait != flagAliasMap.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             /* remap the flag */" << EOL << std::endl;
  sout << "             newFlag = (*ait).second;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           else if (dait != deprecatedFlagAliasMap.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             std::cout << \"Flag \\\"\" << newFlag << \"\\\" is deprecated. Please use flag \\\"\" << (*dait).second << \"\\\" instead. \" << std::endl;" << EOL << std::endl;
  sout << "             /* remap the flag */" << EOL << std::endl;
  sout << "             newFlag = (*dait).second;" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  //   Compiling - i.e. add the boolean flag if not already there.
  sout << "           dvOptionnalArgsIt = std::find(deserializedVectorFlaggedArgs.begin(), deserializedVectorFlaggedArgs.end(), newFlag);" << EOL << std::endl;
  sout << "           /*These flags are always boolean, just add it if it's not there already */" << EOL << std::endl;
  sout << "           if (dvOptionnalArgsIt == deserializedVectorFlaggedArgs.end())" << EOL << std::endl;
  sout << "             {" << EOL << std::endl;
  sout << "             deserializedVectorFlaggedArgs.push_back(newFlag);" << EOL << std::endl;
  sout << "             }" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         ++ac;" << EOL << std::endl;
  sout << "         }" << EOL << std::endl;
  sout << "       else" << EOL << std::endl;
  sout << "         {" << EOL << std::endl;
  // No flags args - Replace if needed, otherwise append to the list.
  sout << "         /* Replace if needed, otherwise append.*/" << EOL << std::endl;
  sout << "         if (noFlagCounter < deserializedVectorPositionalArgs.size())" << EOL << std::endl;
  sout << "           {" << EOL << std::endl;
  sout << "           deserializedVectorPositionalArgs[noFlagCounter] = argv[ac];" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         else" << EOL << std::endl;
  sout << "           {" << EOL << std::endl;
  sout << "           deserializedVectorPositionalArgs.push_back(argv[ac]);" << EOL << std::endl;
  sout << "           }" << EOL << std::endl;
  sout << "         ++ac;" << EOL << std::endl;
  sout << "         ++noFlagCounter;" << EOL << std::endl;
  sout << "         }" << EOL << std::endl;
  sout << "       }" << EOL << std::endl;
  sout << EOL << std::endl;

  sout << "    /* Put the now compiled arguments in the argvVector */" << EOL << std::endl;
  sout << "    std::vector<std::string> argvVector;" << EOL << std::endl;
  sout << "    argvVector.push_back(argv[0]);" << EOL << std::endl;
  sout << "    argvVector.insert(argvVector.end(), deserializedVectorFlaggedArgs.begin(), deserializedVectorFlaggedArgs.end());" << EOL << std::endl;
  sout << "    std::map<std::string, std::vector<std::string> >::iterator mavit;" << EOL << std::endl;
  sout << "    for (mavit = deserializedMultipleArgsMap.begin(); mavit != deserializedMultipleArgsMap.end(); ++mavit)" << EOL << std::endl;
  sout << "      {" << EOL << std::endl;
  sout << "      for (size_t i = 0; i < mavit->second.size(); ++i)" << EOL << std::endl;
  sout << "        {" << EOL << std::endl;
  sout << "        argvVector.push_back(mavit->first);" << EOL << std::endl;
  sout << "        argvVector.push_back(mavit->second.at(i));" << EOL << std::endl;
  sout << "        }" << EOL << std::endl;
  sout << "      }" << EOL << std::endl;
  sout << "    argvVector.insert(argvVector.end(), deserializedVectorPositionalArgs.begin(), deserializedVectorPositionalArgs.end());" << EOL << std::endl;
  sout << EOL << std::endl;
  sout << "   /* Remap args to a structure that CmdLine::parse() can understand*/" << EOL << std::endl;
  sout << "   std::vector<char*> vargs;" << EOL << std::endl;
  sout << "   for (ac = 0; ac < argvVector.size(); ++ac)" << EOL << std::endl;
  sout << "     { " << EOL << std::endl;
  sout << "     vargs.push_back(const_cast<char *>(argvVector[ac].c_str()));" << EOL << std::endl;
  sout << "     }" << EOL << std::endl;

  // Generate the code to parse the command line
  sout << "    commandLine.parse ( static_cast<int>(vargs.size()), (char**) &(vargs[0]) );" << EOL << std::endl;

  // Wrapup the block and generate the catch block
  sout << "  }" << EOL << std::endl;
  sout << "catch ( TCLAP::ArgException & e )" << EOL << std::endl;
  sout << "  {" << EOL << std::endl;
  sout << "  std::cerr << \"error: \" << e.error() << \" for arg \" << e.argId() << std::endl;" << EOL << std::endl;
  sout << "  return ( EXIT_FAILURE );" << EOL << std::endl;
  sout << "  }" << std::endl;
}

void GenerateTCLAPAssignment(std::ostream & sout, const ModuleDescription & module, bool onlyIfSet)
{
  std::string EOL(" \\");
  sout << "#define GENERATE_TCLAP_ASSIGNMENT";
  if( onlyIfSet )
    {
    sout << "_IFSET \\" << std::endl;
    }
  else
    {
    sout << " \\" << std::endl;
    }

  std::vector<ModuleParameterGroup>::const_iterator git;
  std::vector<ModuleParameter>::const_iterator pit;
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if ((*pit).IsReturnParameter())
        {
        continue;
        }

      if( onlyIfSet )
        {
        sout << "    if( " << pit->GetName() << "Arg.isSet() )" << EOL << '\n'
             << "      {" << EOL << '\n';
        }
      sout << "      "
           << pit->GetName();
      if (NeedsTemp(*pit))
        {
        sout << "Temp";
        }
      sout << " = "
           << pit->GetName()
           << "Arg.getValue();"
           << EOL << std::endl;
      if( onlyIfSet )
        {
        sout << "      }" << EOL << '\n';
        }
      }
    }

  // Finally, for any arrays, split the strings into words
  for (git = module.GetParameterGroups().begin();
       git != module.GetParameterGroups().end();
       ++git)
    {
    for (pit = git->GetParameters().begin();
         pit != git->GetParameters().end();
         ++pit)
      {
      if ((*pit).IsReturnParameter())
        {
        continue;
        }

      if (NeedsTemp(*pit) && pit->GetMultiple() != "true")
        {
        if( onlyIfSet )
          {
          sout << "    if( " << pit->GetName() << "Arg.isSet() )" << EOL << '\n';
          }
        sout << "      { " << "/* Assignment for " << pit->GetName() << " */" << EOL << std::endl;
        sout << "      " << pit->GetName() << ".clear();" << EOL << std::endl;
        sout << "      std::vector<std::string> words;"
             << EOL << std::endl;
        sout << "      std::string sep(\",\");"
             << EOL << std::endl;
        if ((*pit).GetTag() == "file")
          {
          sout << "      splitFilenames("
               << pit->GetName()
               << "Temp"
               << ", "
               << "words);"
               << EOL << std::endl;
          }
        else
          {
          sout << "      splitString("
               << pit->GetName()
               << "Temp"
               << ", "
               << "sep, "
               << "words);"
               << EOL << std::endl;
          }
        sout << "      for (unsigned int _j = 0; _j < words.size(); _j++)"
             << EOL << std::endl;
        sout << "        {"
             << EOL << std::endl;
        sout << "        "
             << pit->GetName() << ".push_back("
             << pit->GetStringToType()
             << "(words[_j].c_str()));"
             << EOL << std::endl;
        sout << "        }"
             << EOL << std::endl;
        sout << "      }"
             << EOL << std::endl;
        }
      else if (NeedsTemp(*pit) && pit->GetMultiple() == "true")
        {
        sout << "      { " << "/* Assignment for " << pit->GetName() << " */" << EOL << std::endl;
        sout << "      for (unsigned int _i = 0; _i < ";
        sout << pit->GetName() << "Temp.size(); _i++)" << EOL << std::endl;
        sout << "        {" << EOL << std::endl;
        sout << "        std::vector<std::string> words;" << EOL << std::endl;
        sout << "        words.clear();" << EOL << std::endl;
        if ((*pit).GetTag() == "file")
          {
          sout << "        splitFilenames(" << pit->GetName() << "Temp[_i], words);" << EOL << std::endl;
          }
        else
          {
          sout << "      std::string sep(\",\");" << EOL << std::endl;
          sout << "        splitString(" << pit->GetName() << "Temp[_i], sep, words);" << EOL << std::endl;
          }
        if (IsVectorOfVectors(*pit))
          {
          sout << "        std::vector<" << pit->GetArgType() << "> elements;" << EOL << std::endl;
          sout << "        for (unsigned int _j= 0; _j < words.size(); _j++)" << EOL << std::endl;
          sout << "          {" << EOL << std::endl;
          sout << "          elements.push_back("
               << pit->GetStringToType()
               << "(words[_j].c_str()));"
               << EOL << std::endl;
          sout << "          }" << EOL << std::endl;
          sout << "        " << pit->GetName() << ".push_back(elements);"
               << EOL << std::endl;
          }
        else
          {
          sout << "        for (unsigned int _j= 0; _j < words.size(); _j++)" << EOL << std::endl;
          sout << "          {" << EOL << std::endl;
          sout << "            " << pit->GetName() << ".push_back("
               << pit->GetStringToType()
               << "(words[_j].c_str()));"
               << EOL << std::endl;
          sout << "          }" << EOL << std::endl;
          }
        sout << "        }" << EOL << std::endl;
        sout << "      }"
             << EOL << std::endl;
        }
      }
    }
  sout << "" << std::endl;
}

void GeneratePost(std::ostream &sout)
{
  sout << "#define PARSE_ARGS "
       << "GENERATE_LOGO;"
       << "GENERATE_XML;"
       << "GENERATE_DECLARE;"
#ifdef GenerateCLP_USE_JSONCPP
       << "GENERATE_DESERIALIZATION;"
#endif // GenerateCLP_USE_JSONCPP
       << "GENERATE_TCLAP;"
       << "GENERATE_ECHOARGS;"
       << "GENERATE_ProcessInformationAddressDecoding;"
#ifdef GenerateCLP_USE_JSONCPP
       << "GENERATE_SERIALIZATION;"
#endif // GenerateCLP_USE_JSONCPP
       << std::endl;
}

void GenerateProcessInformationAddressDecoding(std::ostream &sout)
{
  std::string EOL(" \\");
  sout << "#define GENERATE_ProcessInformationAddressDecoding \\" << std::endl;

  sout << "ModuleProcessInformation *CLPProcessInformation = 0;" << EOL << std::endl;
  sout << "if (processInformationAddressString != \"\")" << EOL << std::endl;
  sout << "{" << EOL << std::endl;
  sout << "sscanf(processInformationAddressString.c_str(), \"%p\", &CLPProcessInformation);" << EOL << std::endl;
  sout << "}" << std::endl;
}

#ifdef GENERATECLP_USE_MD5
// Adapted from CMake/Source/cmSystemTools.cxx
bool ComputeFileMD5(const char* source, char* md5out)
{
  if(!itksys::SystemTools::FileExists(source))
    {
    return false;
    }

  // Open files
#if defined(_WIN32) || defined(__CYGWIN__)
  std::ifstream fin(source, std::ios::binary | std::ios::in);
#else
  std::ifstream fin(source);
#endif
  if(!fin)
    {
    return false;
    }

  itksysMD5* md5 = itksysMD5_New();
  itksysMD5_Initialize(md5);

  // Should be efficient enough on most system:
  const int bufferSize = 4096;
  char buffer[bufferSize];
  unsigned char const* buffer_uc =
    reinterpret_cast<unsigned char const*>(buffer);
  // This copy loop is very sensitive on certain platforms with
  // slightly broken stream libraries (like HPUX).  Normally, it is
  // incorrect to not check the error condition on the fin.read()
  // before using the data, but the fin.gcount() will be zero if an
  // error occurred.  Therefore, the loop should be safe everywhere.
  while(fin)
    {
    fin.read(buffer, bufferSize);
    if(int gcount = static_cast<int>(fin.gcount()))
      {
      itksysMD5_Append(md5, buffer_uc, gcount);
      }
    }
  itksysMD5_FinalizeHex(md5, md5out);
  itksysMD5_Delete(md5);

  fin.close();
  return true;
}

// Adapted from CMake/Source/cmSystemTools.cxx
std::string ComputeStringMD5(const char* input)
{
  char md5out[32];
  itksysMD5* md5 = itksysMD5_New();
  itksysMD5_Initialize(md5);
  itksysMD5_Append(md5, reinterpret_cast<unsigned char const*>(input), -1);
  itksysMD5_FinalizeHex(md5, md5out);
  itksysMD5_Delete(md5);
  return std::string(md5out, 32);
}
#endif

} // end of anonymous namespace
