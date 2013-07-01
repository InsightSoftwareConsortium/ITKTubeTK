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

# - Macro to add modules in TubeTK.
# This macro takes a list of modules directory in input.
# Example:
# set( My_list_of_modules m1 m2 )
# include( TubeTKMacroAddModules.cmake )
# TubeTKAddModules( MODULES ${My_list_of_modules} )
#
# It adds the option TubeTK_BUILD_ALL_MODULES and a dependent option
# TubeTK_BUILD_yourmodulename for each given module. This option can then be
# queried by to know if the user want to build the module or not.
# Example (following the example above):
#
# if( TubeTK_BUILD_m1 )
#   # Then I can build m1
# endif( TubeTK_BUILD_m1 )
#
# All the modules chosen to be built are added to the TubeTK_MODULES property

include( CMakeParseArguments )

macro( TubeTKAddModules )
  set( options )
  set( oneValueArgs )
  set( multiValueArgs MODULES )
  CMAKE_PARSE_ARGUMENTS( MY_TUBETK
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN} )

  # Add Build all modules option
  option( TubeTK_BUILD_ALL_MODULES "Build all TubeTK modules or not" ON )
  mark_as_advanced( TubeTK_BUILD_ALL_MODULES )

  # Add option for each module
  if( MY_TUBETK_MODULES )
    include( CMakeDependentOption )

    set( tubetk_modules )
    foreach( module ${MY_TUBETK_MODULES} )
      CMAKE_DEPENDENT_OPTION(
        TubeTK_BUILD_${module} "Build ${module} or not. This does not take module dependencies into account." OFF
         "NOT TubeTK_BUILD_ALL_MODULES" ON )
      mark_as_advanced( TubeTK_BUILD_${module} )

      if( TubeTK_BUILD_${module} )
        list( APPEND tubetk_modules ${CMAKE_CURRENT_SOURCE_DIR}/${module} )
      endif( TubeTK_BUILD_${module} )
    endforeach( module ${MY_TUBETK_MODULES} )

    set_property( GLOBAL APPEND PROPERTY TubeTK_MODULES ${tubetk_modules} )
  endif( MY_TUBETK_MODULES )

endmacro( TubeTKAddModules )
