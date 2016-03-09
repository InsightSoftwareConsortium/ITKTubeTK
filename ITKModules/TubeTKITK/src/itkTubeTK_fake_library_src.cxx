/* This file is intended to remove compilation warning when building against
an ITK build tree.
Without this file, the TubeTKITK target should be link as an INTERFACE
library against TubeTKFiltering. But when building it as an ITK external
module, it will be linked against ITKCommon using the LINK_PRIVATE or
LINK_PUBLIC signatures because ITK doesn't use modern Cmake yet.
This is creating the following conflict depending on the policy :

  Target \"TubeTKITK\" has an INTERFACE_LINK_LIBRARIES property.  This should
  be preferred as the source of the link interface for this library but
  because CMP0022 is not set CMake is ignoring the property and using the
  link implementation as the link interface instead.

  INTERFACE_LINK_LIBRARIES:
    TubeTKFiltering;ITKCommon
  Link implementation:
    ITKCommon

We can avoid this warning by linking TubeTKITK against TubeTKFiltering using
the PUBLIC signature.
Then this file is needed to remove the following warning :

  You have called ADD_LIBRARY for library TubeTKITK without any source files.
  This typically indicates a problem with your CMakeLists.txt file
*/
