/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeCLIFilterWatcher_h
#define __tubeCLIFilterWatcher_h

#include <itkSimpleFilterWatcher.h>

#include <ModuleProcessInformation.h>

namespace tube
{

/** \class CLIFilterWatcher
 * \brief Simple mechanism for monitoring the pipeline events of a
 * filter and reporting these events to std::cout. Formats reports
 * with xml.
 */
class CLIFilterWatcher : public itk::SimpleFilterWatcher
{
public:
  CLIFilterWatcher( itk::ProcessObject* o,
                      const char *comment="",
                      ModuleProcessInformation *inf=0,
                      double fraction = 1.0,
                      double start = 0.0,
                      bool useStdCout = false )
    : SimpleFilterWatcher( o, comment )
    {
    m_ProcessInformation = inf;
    m_Fraction = fraction;
    m_Start = start;
    m_StartCalled = false;
    m_UseStdCout = useStdCout;
    }

protected:

  /** Callback method to show the ProgressEvent */
  virtual void ShowProgress( void )
    {
    if( !m_StartCalled )
      {
      this->StartFilter();
      }

    if( this->GetProcess() )
      {
      this->SetSteps( this->GetSteps()+1 );
      if( !this->GetQuiet() )
        {
        if( m_ProcessInformation )
          {
          std::strncpy( m_ProcessInformation->ProgressMessage,
                  this->GetComment().c_str(), 1023 );
          m_ProcessInformation->Progress =
            ( this->GetProcess()->GetProgress() * m_Fraction + m_Start );
          if( m_Fraction != 1.0 )
            {
            m_ProcessInformation->StageProgress =
              this->GetProcess()->GetProgress();
            }

          this->GetTimeProbe().Stop();
          m_ProcessInformation->ElapsedTime = this->GetTimeProbe().GetMean()
            * this->GetTimeProbe().GetNumberOfStops();
          this->GetTimeProbe().Start();

          if( m_ProcessInformation->Abort )
            {
            this->GetProcess()->AbortGenerateDataOn();
            m_ProcessInformation->Progress = 0;
            m_ProcessInformation->StageProgress = 0;
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
                    << ( this->GetProcess()->GetProgress() * m_Fraction )
                       + m_Start
                    << "</filter-progress>"
                    << std::endl;
          if( m_Fraction != 1.0 )
            {
            std::cout << "<filter-stage-progress>"
                      << this->GetProcess()->GetProgress()
                      << "</filter-stage-progress>"
                      << std::endl;
            }
          std::cout << std::flush;
          }
        }
      }
    }

  /** Callback method to show the StartEvent */
  virtual void StartFilter( void )
    {
    this->SetSteps( 0 );
    this->SetIterations( 0 );
    this->GetTimeProbe().Start();
    this->m_StartCalled = true;
    if( !this->GetQuiet() )
      {
      if( m_ProcessInformation )
        {
        m_ProcessInformation->Progress = 0;
        m_ProcessInformation->StageProgress = 0;
        std::strncpy( m_ProcessInformation->ProgressMessage,
                this->GetComment().c_str(), 1023 );

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
                  << ( this->GetProcess()
                      ? this->GetProcess()->GetNameOfClass() : "None" )
                  << "</filter-name>"
                  << std::endl;
        std::cout << "<filter-comment>"
                  << " \"" << this->GetComment() << "\" "
                  << "</filter-comment>"
                  << std::endl;
        std::cout << "</filter-start>"
                  << std::endl;
        std::cout << std::flush;
        }
      }
    }

  /** Callback method to show the EndEvent */
  virtual void EndFilter( void )
    {
    this->GetTimeProbe().Stop();
    if( !this->GetQuiet() )
      {
      if( m_ProcessInformation )
        {
        m_ProcessInformation->Progress = 1;
        m_ProcessInformation->StageProgress = 1;

        m_ProcessInformation->ElapsedTime = this->GetTimeProbe().GetMean()
          * this->GetTimeProbe().GetNumberOfStops();

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
                  << ( this->GetProcess()
                      ? this->GetProcess()->GetNameOfClass() : "None" )
                  << "</filter-name>"
                  << std::endl;
        std::cout << "<filter-time>"
                  << this->GetTimeProbe().GetMean()
                  << "</filter-time>"
                  << std::endl;
        std::cout << "</filter-end>";
        std::cout << std::flush;
        }
      }
    }

  ModuleProcessInformation * m_ProcessInformation;

  double m_Fraction;
  double m_Start;
  bool   m_StartCalled;
  bool   m_UseStdCout;

}; // End class CLIFilterWatcher

} // End namespace itk

#endif // End !defined( __tubeCLIFilterWatcher_h )
