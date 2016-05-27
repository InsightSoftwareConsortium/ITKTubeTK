# - Find python libraries
# returns:
# - ${library}_ERROR - TRUE if more than one argument is given to the macro
# - ${library}_FOUND - TRUE if the library could be found
#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
# Copyright 2012 Continuum Analytics, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
# nor the names of their contributors may be used to endorse or promote
# products derived from this software without specific prior written
# permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

include(CMakeParseArguments)

macro(TubeTKCheckPythonLibrary library)
  set(${library}_ERROR TRUE)
  set(${library}_FOUND FALSE)
  #Sanity check: Make sure only one library name is given
  if(${ARGC} EQUAL 1)
    set(list_var "${ARGN}")
    find_package(PythonInterp REQUIRED)

    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
      "exec(\"import sys\\ntry:\\n  import ${library}\\nexcept Exception as e:\\n  print(e)\\n  sys.exit(0)\\nsys.exit(1)\")"
      RESULT_VARIABLE ${library}_FOUND
      OUTPUT_VARIABLE _PYTHON_VALUES
      ERROR_VARIABLE _PYTHON_ERROR_VALUE
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(${library}_ERROR FALSE)
  endif()
endmacro()

macro(TubeTKCheckPythonLibraries)
  set(options REQUIRED)
  set(oneValueArgs ERROR_MESSAGE FOUND_VARNAME)
  set(multiValueArgs LIBRARIES)
  cmake_parse_arguments(CheckLibraries "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
   if(CheckLibraries_UNPARSED_ARGUMENTS)
    message(AUTHOR_WARNING "Unparsed arguments given [${CheckLibraries_UNPARSED_ARGUMENTS}]")
  endif()
  if(NOT DEFINED CheckLibraries_LIBRARIES)
    message(AUTHOR_WARNING "LIBRARIES are mandatory: [${CheckLibraries_LIBRARIES}]")
  endif()
  if(DEFINED CheckLibraries_FOUND_VARNAME)
    set(${CheckLibraries_FOUND_VARNAME} ON)
  endif()
  find_package(PythonInterp)
  if(NOT PythonInterp_FOUND)
    if(DEFINED CheckLibraries_FOUND_VARNAME)
      set(${CheckLibraries_FOUND_VARNAME} OFF)
    endif()
    if(CheckLibraries_REQUIRED)
      message(FATAL_ERROR
        "Python interpreter was not found. ${CheckLibraries_ERROR_MESSAGE}")
    endif()
  else()
    foreach(var ${CheckLibraries_LIBRARIES})
      TubeTKCheckPythonLibrary(${var})
      if( NOT ${var}_FOUND )
        if(DEFINED CheckLibraries_FOUND_VARNAME)
          set(${CheckLibraries_FOUND_VARNAME} OFF)
        endif()
        if(CheckLibraries_REQUIRED)
          message(FATAL_ERROR
            "${var} was not found. ${CheckLibraries_ERROR_MESSAGE}")
        endif()
      endif()
    endforeach()
  endif()
endmacro()
