/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageCompareCommand.cxx, v $
  Language:  C++
  Date:      $Date: 2008-11-09 18:18:52 $
  Version:   $Revision: 1.12 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkNumericTraits.h"

#include "metaCommand.h"

#include <iostream>


int RegressionTestFile ( const char *testFilename,
  const char *baselineFilename, bool reportErrors,
  bool createDifferenceFile, double valueTolerance, 
  int numberOfDifferenceTolerance );

int main( int argc, char **argv )
{
  int bestBaselineStatus = 2001;

  // Process some command-line arguments intended for BatchMake
  MetaCommand command;

  // Option for setting the tolerable difference in intensity values
  // between the two images.
  command.SetOption( "toleranceValue", "d", false, 
    "Acceptable differences in numeric values in the files" );
  command.AddOptionField( "toleranceValue", "value", MetaCommand::FLOAT,
    true );

  // Option for setting the number of pixel that can be tolerated to 
  // have different intensities.
  command.SetOption( "toleranceNumberOfDifferences", "n", false, 
    "Number of differences that are acceptable to have intensity differences" );
  command.AddOptionField( "toleranceNumberOfDifferences", "value",
    MetaCommand::INT, true );

  // Option for setting the filename of the test image.
  command.SetOption( "testFile", "t", true, 
    "Filename of the text file to be tested against the baseline images" );
  command.AddOptionField( "testFile", "filename", MetaCommand::STRING,
    true );

  // Option for setting the filename of multiple baseline images.
  command.SetOption( "baselineFiles", "B", false, 
    "List of baseline text files <N> <image1> <image2>...<imageN>" );
  command.AddOptionField( "baselineFiles", "filename", MetaCommand::LIST,
    true );

  // Option for setting the filename of a single baseline image.
  command.SetOption( "baselineFile", "b", false, 
    "Baseline text file filename" );
  command.AddOptionField( "baselineFile", "filename", MetaCommand::STRING,
    true );

  command.Parse( argc, argv );

  double        toleranceValue               = 0.0;
  unsigned long toleranceNumberOfDifferences = 0;
  std::string   testFilename;
  std::string   baselineFilename;

  // If a value tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceValue" ) )
    {
    toleranceValue = command.GetValueAsFloat( "toleranceValue",
      "value" );
    }
 
  // If a number of differences tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceNumberOfDifferences" ) )
    {
    toleranceNumberOfDifferences = command.GetValueAsInt(
      "toleranceNumberOfDifferences", "value" );
    }
     
  // Get the filename of the image to be tested
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
     
  // Get the filename of the base line file
  if( command.GetOptionWasSet( "baselineFile" ) )
    {
    singleBaselineFile = true;
    baselineFilename = command.GetValueAsString( "baselineFile",
      "filename" );
    }

  // Get the filename of the base line image
  if( command.GetOptionWasSet( "baselineFiles" ) )
    {
    singleBaselineFile = false;
    baselineFilenames = command.GetValueAsList( "baselineFiles" );
    }
  
  std::string bestBaselineFilename; 

  try
    {
    if( singleBaselineFile ) 
      {
      bestBaselineStatus = 
        RegressionTestFile( 
            testFilename.c_str(), baselineFilename.c_str(), 
            false, false, toleranceValue, toleranceNumberOfDifferences );
      bestBaselineFilename = baselineFilename;
      }
    else
      {
      typedef std::list< std::string >::const_iterator  nameIterator;
      nameIterator baselineFileItr = baselineFilenames.begin();
      while( baselineFileItr != baselineFilenames.end() )
        {
        const int currentStatus =
          RegressionTestFile( 
              testFilename.c_str(), baselineFileItr->c_str(), 
              false, false, toleranceValue, toleranceNumberOfDifferences );
        if( currentStatus < bestBaselineStatus )
          {
          bestBaselineStatus = currentStatus;
          bestBaselineFilename = *baselineFileItr;
          }
        if( bestBaselineStatus == 0 )
          {
          break;
          }
        ++baselineFileItr;
        }
      }
    // generate images of our closest match
    if( bestBaselineStatus == 0 )
      {
      RegressionTestFile( 
        testFilename.c_str(), bestBaselineFilename.c_str(), true, false, 
        toleranceValue, toleranceNumberOfDifferences );
      }
    else
      {
      RegressionTestFile( 
        testFilename.c_str(), 
        bestBaselineFilename.c_str(), true, true, 
        toleranceValue, toleranceNumberOfDifferences );
      }
    
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Compare text files caught an ITK exception:\n";
    std::cerr << e << "\n";
    bestBaselineStatus = -1;
    }
  catch( const std::exception& e )
    {
    std::cerr << "Compare text files caught an exception:\n";
    std::cerr << e.what() << "\n";
    bestBaselineStatus = -1;
    }
  catch( ... )
    {
    std::cerr << "Compare text files caught an unknown exception!!!\n";
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
}

// Regression Testing Code
int RegressionTestFile ( const char *testFilename,
  const char *baselineFilename, bool reportErrors,
  bool createDifferenceFile, double valueTolerance, 
  int numberOfDifferenceTolerance )
{
  std::ifstream baselineFin( baselineFilename );
  if( ! baselineFin.is_open() )
    {
    std::cerr << "Cannot load baseline = " << baselineFilename << std::endl;
    return EXIT_FAILURE;
    }

  // Read the file generated by the test
  std::ifstream testFin( testFilename );
  if( ! testFin.is_open() )
    {
    std::cerr << "Cannot load test = " << testFilename << std::endl;
    return EXIT_FAILURE;
    }

  std::ofstream diffFout;
  if( createDifferenceFile )
    {
    std::string diffFilename( testFilename );
    diffFilename += std::string( ".diff" );
    diffFout.open( diffFilename.c_str() );
    if( ! diffFout.is_open() )
      {
      std::cerr << "Cannot create difference file = " << diffFilename 
        << std::endl;
      return EXIT_FAILURE;
      }
    }

  unsigned int differences = 0;
  std::string baselineStr;
  std::string testStr;
  while( testFin >> testStr && baselineFin >> baselineStr )
    {
    if( baselineStr.compare( testStr ) != 0 )
      {
      std::istringstream testInStream( testStr );
      std::istringstream baselineInStream( baselineStr );
      double testVal = 0.0;
      double baselineVal = 0.0;
      if( testInStream >> testVal && baselineInStream >> baselineVal )
        {
        if( fabs( testVal - baselineVal ) > valueTolerance )
          {
          if( reportErrors )
            {
            std::cout << testVal << " != " << baselineVal << std::endl;
            }
          if( createDifferenceFile )
            {
            diffFout << testVal << " != " << baselineVal << std::endl;
            }
          ++differences;
          }
        }
      else
        {
        if( reportErrors )
          {
          std::cout << testStr << " != " << baselineStr << std::endl;
          }
        if( createDifferenceFile )
          {
          diffFout << testStr << " != " << baselineStr << std::endl;
          }
        ++differences;
        }
      }
    }

  if( numberOfDifferenceTolerance >= 0 &&
    differences > (unsigned int)(numberOfDifferenceTolerance) )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
