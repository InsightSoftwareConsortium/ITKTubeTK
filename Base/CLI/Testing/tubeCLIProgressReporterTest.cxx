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

#include "tubeCLIProgressReporter.h"

void cliProgressReporterTestCallbackFunction( void * data )
{
  if( data )
    {
    std::cout << "cliProgressReporterTestCallbackFunction with pointer: "
              << data << std::endl;
    }
  else
    {
    std::cout << "cliProgressReporterTestCallbackFunction with pointer: NULL"
              << std::endl;
    }
}

int tubeCLIProgressReporterTest( int argc, char * itkNotUsed( argv )[] )
{
  if( argc != 1 )
    {
    std::cerr << "Usage: " << std::endl;
    return EXIT_FAILURE;
    }

  ModuleProcessInformation CLPProcessInfo;
  CLPProcessInfo.Initialize();
  CLPProcessInfo.SetProgressCallback( 
    cliProgressReporterTestCallbackFunction, NULL );

  std::string processName( "tubeCLIProgressReporterTest" );

  tube::CLIProgressReporter reporter( processName.c_str(),
                                      &CLPProcessInfo,
                                      true );
  reporter.Start();
  reporter.Report( 0.0 );
  reporter.Report( 0.5 );
  CLPProcessInfo.Abort = 1;
  reporter.Report( 1.0 );
  reporter.End();
  if( reporter.GetProcess() != processName )
    {
    std::cerr << "Process name does not match." << std::endl;
    std::cerr << reporter.GetProcess() << " != " << processName
              << std::endl;
    }
  if( reporter.GetProcessInformation() != &CLPProcessInfo )
    {
    std::cerr << "Process info does not match." << std::endl;
    }

  return EXIT_SUCCESS;
}
