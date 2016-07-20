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
#
# Example:
#
#   set( My_list_of_modules m1 m2 )
#   include( TubeTKMacroAddModules )
#   TubeTKMacroAddModules( MODULES ${My_list_of_modules} )
#
#
# TubeTK_BUILD_ALL_MODULES:
#
# It adds the option TubeTK_BUILD_ALL_MODULES and a dependent option
# TubeTK_BUILD_yourmodulename for each given module. This option can then be
# queried by to know if the user want to build the module or not.
#
# By default, there are no "TubeTK_BUILD_yourmodulename" options. Setting
# TubeTK_BUILD_ALL_MODULES to OFF will make them available.
#
# Example (following the example above):
#
#   if( TubeTK_BUILD_m1 )
#     # Then I can build m1
#   endif( TubeTK_BUILD_m1 )
#
#
# TubeTK_BUILD_MODULE_DEFAULT:
#
# The default value of "TubeTK_BUILD_yourmodulename"  options will be ON, this
# can be changed by configuring TukeTK with option TubeTK_BUILD_MODULE_DEFAULT
# set to OO. Note that option 'TubeTK_BUILD_MODULE_DEFAULT' has to be specified
# when configuring for the first time (clean build).
#
#
# TubeTK_MODULES:
#
# All the modules chosen to be built are added to the TubeTK_MODULES property

include( CMakeParseArguments )

macro( TubeTKMacroAddModules )
  set( options )
  set( oneValueArgs )
  set( multiValueArgs MODULES )
  CMAKE_PARSE_ARGUMENTS( MY_TubeTK
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN} )

  set( _build_module_default ON )
  if( DEFINED TubeTK_BUILD_MODULE_DEFAULT )
    if( TubeTK_BUILD_ALL_MODULES )
      message( AUTHOR_WARNING
        "Specifying 'TubeTK_BUILD_MODULE_DEFAULT' is effective only if 'TubeTK_BUILD_ALL_MODULES' is disabled" )
    endif()
    set( _build_module_default ${TubeTK_BUILD_MODULE_DEFAULT} )
  endif()

  # Add Build all modules option
  option( TubeTK_BUILD_ALL_MODULES "Build all TubeTK modules or not" ON )
  mark_as_advanced( TubeTK_BUILD_ALL_MODULES )
  if( TubeTK_USE_SUPERBUILD )
    mark_as_superbuild( TubeTK_BUILD_ALL_MODULES )
  endif()

  # Add option for each module
  if( MY_TubeTK_MODULES )
    include( CMakeDependentOption )

    set( tubetk_modules )
    foreach( module ${MY_TubeTK_MODULES} )

      CMAKE_DEPENDENT_OPTION(
        TubeTK_BUILD_${module}
        "Build ${module}. Does not take module dependencies into account." OFF
        "NOT TubeTK_BUILD_ALL_MODULES" ${_build_module_default} )
      mark_as_advanced( TubeTK_BUILD_${module} )
      if( TubeTK_USE_SUPERBUILD )
        mark_as_superbuild( TubeTK_BUILD_${module} )
      endif()

      if( TubeTK_BUILD_${module} )
        list( APPEND tubetk_modules ${CMAKE_CURRENT_SOURCE_DIR}/${module} )
      endif( TubeTK_BUILD_${module} )
    endforeach( module ${MY_TubeTK_MODULES} )

    set_property( GLOBAL APPEND PROPERTY TubeTK_MODULES ${tubetk_modules} )
  endif( MY_TubeTK_MODULES )

endmacro( TubeTKMacroAddModules )
