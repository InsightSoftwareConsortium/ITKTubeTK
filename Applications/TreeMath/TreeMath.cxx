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
#include <iostream>
#include <sstream>
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "TreeMathCLP.h"
#include <metaCommand.h>

//TubeTK imports
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "tubeTreeFilters.h"

template< unsigned int VDimension >
int DoIt( MetaCommand & command )
{
  if( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D data is currently supported." );
    return EXIT_FAILURE;
    }

  typename itk::GroupSpatialObject< VDimension >::Pointer inputTubes;

  MetaCommand::OptionVector parsed = command.GetParsedOptions();

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;
  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();
  try
    {
    tubeFileReader->SetFileName( command.GetValueAsString( "infile" ).
      c_str() );
    tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE File: " + std::string(
      err.GetDescription() ) );
    return EXIT_FAILURE;
    }

  inputTubes = tubeFileReader->GetGroup();
  inputTubes->ComputeObjectToWorldTransform();

  MetaCommand::OptionVector::const_iterator it = parsed.begin();
  while( it != parsed.end() )
    {
    if( it->name == "Write" )
      {
      tubeStandardOutputMacro( << "\n>> Writing TRE File" );

      typedef itk::SpatialObjectWriter< VDimension > TubeWriterType;
      typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
      try
        {
        tubeWriter->SetFileName( command.GetValueAsString( *it,
          "filename" ).c_str() );
        tubeWriter->SetInput( inputTubes );
        tubeWriter->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        tube::ErrorMessage( "Error writing TRE file: " + std::string(
          err.GetDescription() ) );
        return EXIT_FAILURE;
        }
      }
    else if( it->name == "FillGapsInTubeTree" )
      {
      tubeStandardOutputMacro( << "\n>> Filling gaps in input tree" );
      tube::TreeFilters< VDimension >::FillGap( inputTubes, command.GetValueAsString( *it,
        "InterpolationMethod" ).c_str()[0] );
      }
    ++it;
    }
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  MetaCommand command;

  command.SetName( "TreeMath" );
  command.SetVersion( "1.0" );
  command.SetAuthor( "Sumedha Singla" );
  command.SetDescription( "Perform several filters on a tube tree." );

  command.AddField( "infile", "infile filename",
  MetaCommand::STRING, MetaCommand::DATA_IN );

  command.SetOption( "Write", "w", false,
    "Writes current tubes to the designated file." );
  command.AddOptionField( "Write", "filename", MetaCommand::STRING, true,
    "", "Output filename", MetaCommand::DATA_OUT );

  command.SetOption( "FillGapsInTubeTree", "f", false,
    "Connects the parent and child tube if they have a gap inbetween,"
    " by interpolating the path inbetween." );
  command.AddOptionField( "FillGapsInTubeTree", "InterpolationMethod",
    MetaCommand::STRING, true, "",
    "[S]traight Line, [L]Linear Interp., [C]urve Fitting, [M]inimal Path",
    MetaCommand::DATA_IN );

  if( !command.Parse( argc, argv ) )
    {
    return EXIT_FAILURE;
    }

  return DoIt< 3 >( command );
}
