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

#ifdef WIN32
#include <stdio.h>
#include <process.h>
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

#ifdef WIN32
  char *args[4];
  args[0] = &serverexecutable[0];
  args[1] = &port[0];
  args[2] = "10";
  args[3] = NULL;

  char *args1[4];
  args1[0] = &clientexecutable[0];
  args1[1] = &port[0];
  args1[2] = &hostname[0];
  args1[3] = NULL;


  _spawnv(P_NOWAIT, serverexecutable.c_str(), args);
  
  _spawnv(P_NOWAIT, clientexecutable.c_str(), args1);

#else
  pid_t pid = fork();
  
  if (pid == 0) // parent fork run server sending executable
    {
    std::stringstream command;
    command << serverexecutable;
    command << " ";
    command << port;
    command << " ";
    command << 10;  // duration
    system( command.str().c_str() );
    command.clear();
    }
  else // child fork run client receiving executable
    {
    sleep( 3 );
    std::stringstream command;
    command << clientexecutable;
    command << " ";
    command << port;
    command << " ";
    command << hostname;
    system( command.str().c_str() );
    command.clear();
    }
#endif

  return EXIT_SUCCESS;
}

