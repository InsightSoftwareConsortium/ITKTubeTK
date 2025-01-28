/*=========================================================================

Library:   TubeTKLib

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

#include "tubeObject.h"

#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ > 0)
#  include <cxxabi.h>
#endif

namespace tube
{

// Constructor.
Object::Object(void) {}

// Destructor.
Object::~Object(void) {}

// Print out information about this object.
void
Object::Print(std::ostream & os, Indent indent) const
{
  this->PrintHeader(os, indent);
  this->PrintSelf(os, indent.GetNextIndent());
  this->PrintFooter(os, indent);
}

// Header for when printing out information about this object.
void
Object::PrintHeader(std::ostream & os, Indent indent) const
{
  os << indent << this->GetNameOfClass() << " ( " << this << " )" << std::endl;
}

// Print out information about the member variables of this object.
void
Object::PrintSelf(std::ostream & os, Indent indent) const
{
  os << indent << "RTTI typeinfo: ";

#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ > 0)
  const char * mangled = typeid(*this).name();
  int          status;
  char *       demangled = abi::__cxa_demangle(mangled, NULL, NULL, &status);

  if (status == 0)
  {
    os << demangled;
    free(demangled);
  }
  else
  {
    os << mangled;
  }
#else
  os << typeid(*this).name();
#endif

  os << std::endl;
}

// Footer for when printing out information about this object.
void
Object::PrintFooter(std::ostream & tubeNotUsed(os), Indent tubeNotUsed(indent)) const
{}

// Print out information about the specified object.
std::ostream &
operator<<(std::ostream & os, const Object & object)
{
  object.Print(os);

  return os;
}

} // End namespace tube
