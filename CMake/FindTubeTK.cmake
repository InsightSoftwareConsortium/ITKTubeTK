##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

# Find a TubeTK installation or build tree.
# The following variables are set if TubeTK is found.  If TubeTK is not
# found, TubeTK_FOUND is set to false.
#  TubeTK_FOUND         - Set to true when TubeTK is found.
#  TubeTK_USE_FILE      - CMake file to use TubeTK.
#  TubeTK_MAJOR_VERSION - The TubeTK major version number.
#  TubeTK_MINOR_VERSION - The TubeTK minor version number
#                       (odd non-release).
#  TubeTK_BUILD_VERSION - The TubeTK patch level
#                       (meaningless for odd minor).
#  TubeTK_INCLUDE_DIRS  - Include directories for TubeTK
#  TubeTK_LIBRARY_DIRS  - Link directories for TubeTK libraries
#  TubeTK_LIBRARIES     - List of libraries
#
# The following cache entries must be set by the user to locate TubeTK:
#  TubeTK_DIR  - The directory containing TubeTKConfig.cmake.
#             This is either the root of the build tree,
#             or the lib/tubetk directory.  This is the
#             only cache entry.
#

# Assume not found.
set( TubeTK_FOUND 0 )

# Construct consitent error messages for use below.
set( TubeTK_DIR_MESSAGE
"TubeTK not found.  Set the TubeTK_DIR cmake cache entry to the directory containing TubeTKConfig.cmake.  This is either the root of the build tree or PREFIX/lib/tubetk for an installation." )

# Use the Config mode of the find_package command to find TubeTKConfig.
# If this succeeds (possibly because TubeTK_DIR is already set), the
# command will have already loaded TubeTKConfig.cmake and set TubeTK_FOUND.
if( NOT TubeTK_FOUND )
  find_package( TubeTK QUIET NO_MODULE )
endif( NOT TubeTK_FOUND )

if( NOT TubeTK_DIR )
  find_path( TubeTK_DIR
    NAMES UseTubeTK.cmake
    PATH_SUFFIXES TubeTK tubetk
    HINTS $ENV{TubeTK_DIR}
    PATHS
    ..

    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10] )

endif( NOT TubeTK_DIR )

if( NOT TubeTK_FOUND )
  # TubeTK not found, explain to the user how to specify its location.
  if( TubeTK_FIND_REQUIRED )
    message( FATAL_ERROR ${TubeTK_DIR_MESSAGE} )
  else( TubeTK_FIND_REQUIRED )
    if( NOT TubeTK_FIND_QUIETLY )
      message( STATUS ${TubeTK_DIR_MESSAGE} )
    endif( NOT TubeTK_FIND_QUIETLY )
  endif( TubeTK_FIND_REQUIRED )
endif( NOT TubeTK_FOUND )
