/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

// Compare a test JSON file to a baseline JSON file.  Implementation is
// adapted from CalaTK.

#include <json/json.h>

#include <metaCommand.h>

#include <vnl/vnl_math.h>

/** Return 0 if the test file and baseline file have the same JSON content and 1
 * otherwise. Values within the JSON object can be ignored by setting them to
 * "Regression test NA"
 */
int RegressionTestJson( const char *testJSONFileName,
  const char *baselineJSONFileName,
  bool reportErrors = true,
  bool verbose = true,
  double toleranceValue = 0.00001 );

int main( int argc, char * argv[] )
{
  int bestBaselineStatus = 2001;

  // Process some command-line arguments intended for BatchMake
  MetaCommand command;

  // Option for setting the tolerable difference in floating point values.
  command.SetOption( "toleranceValue", "d", false,
    "Acceptable differences in numeric values in the files" );
  command.AddOptionField( "toleranceValue", "value", MetaCommand::FLOAT,
    true );

  // Option for setting the filename of the test JSON file.
  command.SetOption( "testFile", "t", true,
    "Filename of the JSON file to be tested against the baseline JSON files." );
  command.AddOptionField( "testFile", "filename", MetaCommand::STRING,
    true );

  // Option for setting the filename of multiple baseline JSON documents.
  command.SetOption( "baselineFiles", "B", false,
    "List of baseline JSON files <N> <JSON1> <JSON2>...<JSONN>" );
  command.AddOptionField( "baselineFiles", "filename", MetaCommand::LIST,
    true );

  // Option for setting the filename of a single baseline JSON document.
  command.SetOption( "baselineFile", "b", false,
    "Baseline JSON file filename" );
  command.AddOptionField( "baselineFile", "filename", MetaCommand::STRING,
    true );

  command.Parse( argc, argv );

  double        toleranceValue               = 0.00001;
  std::string   testFilename;
  std::string   baselineFilename;

  // If a value tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceValue" ) )
    {
    toleranceValue = command.GetValueAsFloat( "toleranceValue",
      "value" );
    }

  // Get the filename of the JSON document to be tested
  if( command.GetOptionWasSet( "testFile" ) )
    {
    testFilename =
      command.GetValueAsString( "testFile", "filename" );
    }

  std::list< std::string > baselineFilenames;
  baselineFilenames.clear();

  bool singleBaselineFile = true;

  if( !command.GetOptionWasSet( "baselineFile" )
    && !command.GetOptionWasSet( "baselineFiles" ) )
    {
    std::cerr <<
      "You must provide a -BaselineFile or -BaselineFiles option"
      << std::endl;
    return EXIT_FAILURE;
    }

  // Get the filename of the baseline file
  if( command.GetOptionWasSet( "baselineFile" ) )
    {
    singleBaselineFile = true;
    baselineFilename = command.GetValueAsString( "baselineFile",
      "filename" );
    }

  // Get the filename of the base line JSON documents
  if( command.GetOptionWasSet( "baselineFiles" ) )
    {
    singleBaselineFile = false;
    baselineFilenames = command.GetValueAsList( "baselineFiles" );
    }

  try
    {
    if( singleBaselineFile )
      {
      bestBaselineStatus =
        RegressionTestJson( 
            testFilename.c_str(), baselineFilename.c_str(),
            true, true, toleranceValue );
      }
    else
      {
      typedef std::list< std::string >::const_iterator  nameIterator;
      nameIterator baselineFileItr = baselineFilenames.begin();
      while( baselineFileItr != baselineFilenames.end() )
        {
        const int currentStatus =
          RegressionTestJson( 
              testFilename.c_str(), baselineFileItr->c_str(),
              true, true, toleranceValue );
        if( currentStatus < bestBaselineStatus )
          {
          bestBaselineStatus = currentStatus;
          }
        if( bestBaselineStatus == 0 )
          {
          break;
          }
        ++baselineFileItr;
        }
      }
    }
  catch( const std::exception& e )
    {
    std::cerr << "Compare JSON files caught an exception:\n";
    std::cerr << e.what() << "\n";
    bestBaselineStatus = -1;
    }
  catch( ... )
    {
    std::cerr << "Compare JSON files caught an unknown exception!!!\n";
    bestBaselineStatus = -1;
    }
  if( bestBaselineStatus == 0 )
    {
    std::cout << "Files match." << std::endl;
    return EXIT_SUCCESS;
    }
  else
    {
    std::cout << "Files differ." << std::endl;
    return EXIT_FAILURE;
    }
} // end main


// Recursively compare test and baseline.
int compareJSON( const Json::Value & test,
  const Json::Value & baseline,
  bool reportErrors,
  bool verbose,
  const double & toleranceValue )
{
  int same = EXIT_SUCCESS;
  for( Json::Value::const_iterator baselineIt = baseline.begin(),
         testIt = test.begin();
       testIt != test.end();
       ++baselineIt, ++testIt )
    {
    // Certain points in the hierarchy can be avoided by setting them to
    // "Regression test NA"
    if( ( *baselineIt ).isString() &&
        ( *baselineIt ).asString() == "Regression test NA" )
      {
      continue;
      }
    if( !( ( *baselineIt ).isNull() ) &&
        ( ( *baselineIt ).isArray() || ( *baselineIt ).isObject() ) )
      {
      if( compareJSON( *testIt, *baselineIt, reportErrors, verbose,
        toleranceValue ) )
        {
        same = EXIT_FAILURE;
        break;
        }
      }
    else
      {
      if( verbose &&
          !( *baselineIt ).isNull() && !( *testIt ).isNull() )
        {
        std::cout << "Comparing: " << *testIt
                  << " to " << *baselineIt << std::endl;
        }
      if( ( *baselineIt ).isNull() )
        {
        if( !( *testIt ).isNull() )
          {
          if( reportErrors )
            {
            std::cerr << "The test value was non-null "
                      << "when the baseline value was null." << std::endl;
            }
          same = EXIT_FAILURE;
          break;
          }
        }
      else if( ( *baselineIt ).isBool() )
        {
        if( ( *baselineIt ).asBool() != ( *testIt ).asBool() )
          {
          if( reportErrors )
            {
            std::cerr << "The test value: " << ( *testIt ).asBool()
                      << " does not equal the baseline value: "
                      << ( *baselineIt ).asBool() << std::endl;
            }
          same = EXIT_FAILURE;
          break;
          }
        }
      else if( ( *baselineIt ).isIntegral() )
        {
        if( ( *baselineIt ).asInt() != ( *testIt ).asInt() )
          {
          if( reportErrors )
            {
            std::cerr << "The test value: " << ( *testIt ).asInt()
                      << " does not equal the baseline value: "
                      << ( *baselineIt ).asInt() << std::endl;
            }
          same = EXIT_FAILURE;
          break;
          }
        }
      else if( ( *baselineIt ).isDouble() )
        {
        if( vnl_math_abs( ( *baselineIt ).asDouble() - ( *testIt ).asDouble() )
              > toleranceValue )
          {
          if( reportErrors )
            {
            std::cerr << "The test value: " << ( *testIt ).asDouble()
                      << " does not equal the baseline value: "
                      << ( *baselineIt ).asDouble() << std::endl;
            }
          same = EXIT_FAILURE;
          break;
          }
        }
      else if( ( *baselineIt ).isString() )
        {
        if( ( *baselineIt ).asString() != ( *testIt ).asString() )
          {
          if( reportErrors )
            {
            std::cerr << "The test value: " << ( *testIt ).asString()
                      << " does not equal the baseline value: "
                      << ( *baselineIt ).asString() << std::endl;
            }
          same = EXIT_FAILURE;
          break;
          }
        }
      }
    }
  return same;
}


int RegressionTestJson( const char *testJSONFileName,
  const char *baselineJSONFileName,
  bool reportErrors,
  bool verbose,
  double toleranceValue )
{
  Json::Value testRoot;
  Json::Value baselineRoot;

  Json::Reader reader;
  std::ifstream testFile( testJSONFileName );
  if( !testFile.is_open() )
    {
    std::cerr << "Could not open test file: " << testJSONFileName << std::endl;
    return EXIT_FAILURE;
    }
  if( !reader.parse( testFile, testRoot ) )
    {
    std::cerr << "Could not parse test file: " << testJSONFileName << std::endl;
    testFile.close();
    return EXIT_FAILURE;
    }
  testFile.close();
  std::ifstream baselineFile( baselineJSONFileName );
  if( !baselineFile.is_open() )
    {
    std::cerr << "Could not open baseline file: " << baselineJSONFileName
      << std::endl;
    return EXIT_FAILURE;
    }
  if( !reader.parse( baselineFile, baselineRoot ) )
    {
    std::cerr << "Could not parse baseline file: " << baselineJSONFileName
      << std::endl;
    baselineFile.close();
    return EXIT_FAILURE;
    }
  baselineFile.close();

  return compareJSON( testRoot,
    baselineRoot, reportErrors,
    verbose, toleranceValue );
}
