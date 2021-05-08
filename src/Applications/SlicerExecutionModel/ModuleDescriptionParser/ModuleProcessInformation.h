/*=========================================================================

  Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   Module Description Parser
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __ModuleProcessInformation_h
#define __ModuleProcessInformation_h

#include <cstring>
#include <ostream>
#include <string.h>

#if defined(_WIN32)
#pragma warning ( disable : 4996 )
#endif

extern "C" {
  struct ModuleProcessInformation
  {
    /** Inputs from calling application to the module **/
    unsigned char Abort;

    /** Outputs from the module to the calling application **/
    float Progress;      /// Overall progress
    float StageProgress; /// Progress of a single stage in an algorithm
    char  ProgressMessage[1024];
    void (*ProgressCallbackFunction)(void *);
    void *ProgressCallbackClientData;

    double ElapsedTime;

    ModuleProcessInformation()
      {
        ProgressCallbackFunction = nullptr;
        ProgressCallbackClientData = nullptr;
        Initialize();
      }

    void Initialize()
      {
        Abort = 0;
        Progress = 0;
        StageProgress = 0;
        strcpy(ProgressMessage, "");
        ElapsedTime = 0.0;
      }

    void SetProgressCallback( void (*fun)(void *), void *who )
      {
        ProgressCallbackFunction = fun;
        ProgressCallbackClientData = who;
      }
  };
}

std::ostream& operator<<(std::ostream &os, const ModuleProcessInformation &p);


#endif
