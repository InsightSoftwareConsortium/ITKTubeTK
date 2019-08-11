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

/* This file is intended to remove compilation warning when building against
an ITK build tree.
Without this file, the TubeTK target should be link as an INTERFACE
library against TubeTKFiltering. But when building it as an ITK external
module, it will be linked against ITKCommon using the LINK_PRIVATE or
LINK_PUBLIC signatures because ITK doesn't use modern Cmake yet.
This is creating the following conflict depending on the policy :

Target \"TubeTKITK\" has an INTERFACE_LINK_LIBRARIES property.  This should
be preferred as the source of the link interface for this library but
because CMP0022 is not set CMake is ignoring the property and using the
link implementation as the link interface instead.

INTERFACE_LINK_LIBRARIES: TubeTKFiltering;ITKCommon
Link implementation: ITKCommon

We can avoid this warning by linking TubeTKITK against TubeTKFiltering using
the PUBLIC signature.
Then this file is needed to remove the following warning :

You have called ADD_LIBRARY for library TubeTKITK without any source files.
This typically indicates a problem with your CMakeLists.txt file
*/

namespace itk {

namespace tube {

class EmptyClass
  {
  };

}

}
