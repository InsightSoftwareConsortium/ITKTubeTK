/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// The following three should be used in every CLI application
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"

// Includes specific to this CLI application
#include <iostream>

// Must include CLP before including Helper functions
#include "SendAndReceivePositionMessageUsingOpenIGTLinkTestCLP.h"

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  pid_t pid;

  std::cout << "Starting sender ... " << std::endl;
  pid = fork();
  
  if (pid == 0) // parent fork run server sending executable
    {
    std::stringstream command;
    command << serverexecutable;
    command << " ";
    command << port;
    command << " ";
    command << 2;  // duration
    system( command.str().c_str() );
    command.clear();
    }
  else // child fork run client receiving executable
    {
    std::stringstream command;
    command << clientexecutable;
    command << " ";
    command << port;
    command << " ";
    command << hostname;
    system( command.str().c_str() );
    command.clear();
    }
  
  return EXIT_SUCCESS;
}
