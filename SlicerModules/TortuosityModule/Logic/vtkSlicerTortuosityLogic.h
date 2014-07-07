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

// .NAME vtkSlicerTortuosityLogic -
// \todo
// .SECTION Description
// \todo

#ifndef __vtkSlicerTortuosityLogic_h
#define __vtkSlicerTortuosityLogic_h

#include <vtkSlicerModuleLogic.h>
#include <vtkSlicerTortuosityModuleLogicExport.h>

class VTK_SLICER_TORTUOSITY_MODULE_LOGIC_EXPORT vtkSlicerTortuosityLogic
 : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerTortuosityLogic *New( void );
  vtkTypeRevisionMacro(vtkSlicerTortuosityLogic,vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkSlicerTortuosityLogic( void );
  ~vtkSlicerTortuosityLogic( void );
  vtkSlicerTortuosityLogic(const vtkSlicerTortuosityLogic&);
  void operator=(const vtkSlicerTortuosityLogic&);

}; // End class vtkSlicerTortuosityLogic

#endif // End !defined(__vtkSlicerTortuosityLogic_h)
