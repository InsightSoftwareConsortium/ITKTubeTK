/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 1U24CA194354-01

==============================================================================*/

#ifndef __ModuleDescriptionTestingUtilities_h
#define __ModuleDescriptionTestingUtilities_h

/// This module provides functions to facilitate writing tests.
///
/// Example:
///
/// \code{.cpp}
/// int current = 40 + 2;
/// int expected = 43;
/// if (!CheckInt(__LINE__, "40 + 2", current, expected))
///   {
///   return false;
///   }
/// \endcode
///
/// Usually these test methods are used by single-line convenience macros
/// defined in ModuleDescriptionTestingMacros.h.
///

// Adapted from vtkMRMLCoreTestingUtilities.h

// STD includes
#include <string>

namespace ModuleDescriptionTestingUtilities
{

bool CheckInt(int line, const std::string& description,
              int current, int expected);

bool CheckNotNull(int line, const std::string& description,
                  const void* pointer);

bool CheckNull(int line, const std::string& description,
               const void* pointer);

bool CheckPointer(int line, const std::string& description,
                  void* current, void* expected, bool errorIfDifferent = true);

bool CheckString(int line, const std::string& description,
                 const char* current, const char* expected, bool errorIfDifferent = true );

}

#endif
