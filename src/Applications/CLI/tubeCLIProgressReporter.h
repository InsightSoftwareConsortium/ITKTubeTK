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

#ifndef __tubeCLIProgressReporter_h
#define __tubeCLIProgressReporter_h

#include "tubeMacro.h"

#include <itkTimeProbe.h>

#include <ModuleProcessInformation.h>

namespace tube
{

/** \class CLIProgressReporter
 * \brief Simple mechanism for monitoring the pipeline events of a
 * filter and reporting these events to std::cout. Formats reports
 * with xml.
 */
class CLIProgressReporter
{

public:

  CLIProgressReporter( const char * process,
                       ModuleProcessInformation * inf,
                       bool useStdCout = false )
    {
    m_Process = process;
    m_ProcessInformation = inf;
    m_UseStdCout = useStdCout;
    }

  virtual ~CLIProgressReporter( void )
    {
    }

  virtual bool Report( double fraction )
    {
    if( m_ProcessInformation )
      {
      std::strncpy( m_ProcessInformation->ProgressMessage,
              this->m_Process.c_str(), 1023 );
      m_ProcessInformation->Progress = fraction;
      // if( m_Fraction != 1.0 )
      //   {
      //   m_ProcessInformation->StageProgress =
      //   this->GetProcess()->GetProgress();
      //   }

      m_TimeProbe.Stop();
      m_ProcessInformation->ElapsedTime = m_TimeProbe.GetMean()
                                          * m_TimeProbe.GetNumberOfStops();
      m_TimeProbe.Start();

      if( m_ProcessInformation->Abort )
        {
        m_ProcessInformation->Progress = 0;
        m_ProcessInformation->StageProgress = 0;
        return false;
        }

      if( m_ProcessInformation->ProgressCallbackFunction
          && m_ProcessInformation->ProgressCallbackClientData )
        {
        ( *( m_ProcessInformation->ProgressCallbackFunction ) )(
          m_ProcessInformation->ProgressCallbackClientData );
        }
      }

    if( !m_ProcessInformation || m_UseStdCout )
      {
      std::cout << "<filter-progress>"
                << fraction
                << "</filter-progress>"
                << std::endl;
      // if( m_Fraction != 1.0 )
      //   {
      //   std::cout << "<filter-stage-progress>"
      //             << this->GetProcess()->GetProgress()
      //             << "</filter-stage-progress>"
      //             << std::endl;
      //   }
      // std::cout << std::flush;
      }

    return true;

    }

  /** Callback method to show the StartEvent */
  virtual void Start( void )
    {
    m_TimeProbe.Start();
    if( m_ProcessInformation )
      {
      m_ProcessInformation->Progress = 0;
      m_ProcessInformation->StageProgress = 0;
      std::strncpy( m_ProcessInformation->ProgressMessage,
              m_Process.c_str(), 1023 );

      if( m_ProcessInformation->ProgressCallbackFunction
          && m_ProcessInformation->ProgressCallbackClientData )
        {
        ( *( m_ProcessInformation->ProgressCallbackFunction ) )(
          m_ProcessInformation->ProgressCallbackClientData );
        }
      }

    if( !m_ProcessInformation || m_UseStdCout )
      {
      std::cout << "<filter-start>"
                << std::endl;
      std::cout << "<filter-name>"
                << m_Process
                << "</filter-name>"
                << std::endl;
      std::cout << "<filter-comment>"
                << " \"" << m_Process << "\" "
                << "</filter-comment>"
                << std::endl;
      std::cout << "</filter-start>"
                << std::endl;
      std::cout << std::flush;
      }
    }

  /** Callback method to show the EndEvent */
  virtual void End( void )
    {
    m_TimeProbe.Stop();
    if( m_ProcessInformation )
      {
      m_ProcessInformation->Progress = 1;
      m_ProcessInformation->StageProgress = 1;

      m_ProcessInformation->ElapsedTime = m_TimeProbe.GetMean()
                                          * m_TimeProbe.GetNumberOfStops();

      if( m_ProcessInformation->ProgressCallbackFunction
          && m_ProcessInformation->ProgressCallbackClientData )
        {
        ( *( m_ProcessInformation->ProgressCallbackFunction ) )(
          m_ProcessInformation->ProgressCallbackClientData );
        }
      }

    if( !m_ProcessInformation || m_UseStdCout )
      {
      std::cout << "<filter-end>"
                << std::endl;
      std::cout << "<filter-name>"
                << m_Process
                << "</filter-name>"
                << std::endl;
      std::cout << "<filter-time>"
                << m_TimeProbe.GetMean()
                   * m_TimeProbe.GetNumberOfStops()
                << "</filter-time>"
                << std::endl;
      std::cout << "</filter-end>";
      std::cout << std::flush;
      }
    }

  virtual void SetUseStdCout( bool useStdCout )
    {
    m_UseStdCout = useStdCout;
    }

  virtual bool GetUseStdCout( void )
    {
    return m_UseStdCout;
    }

  virtual std::string GetProcess( void )
    {
    return m_Process;
    }

  virtual ModuleProcessInformation * GetProcessInformation( void )
    {
    return m_ProcessInformation;
    }

protected:

  bool                       m_UseStdCout;
  itk::TimeProbe             m_TimeProbe;
  std::string                m_Process;
  ModuleProcessInformation * m_ProcessInformation;

}; // End class CLIProgressReporter

} // End namespace itk

#endif // End !defined( __tubeCLIProgressReporter_h )
