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

#
# Today macro - Allows to retrieve current date in a cross-platform fashion
# (Unix-like and Windows).
#
# Adapted from the work of Benoit Rat
# See http://www.cmake.org/pipermail/cmake/2009-February/027014.html
#

function( TubeTKFunctionToday RESULT_VAR )
  set( today_date )
  if( WIN32 )
    execute_process( COMMAND cmd /c "date /T"
                     ERROR_VARIABLE getdate_error
                     RESULT_VARIABLE getdate_result
                     OUTPUT_VARIABLE today_date
                     OUTPUT_STRIP_TRAILING_WHITESPACE )

    string( REGEX REPLACE ".*(..)/(..)/(....)"
                          "\\3-\\1-\\2"
                          today_date
                          ${today_date} )
  elseif( UNIX )
    execute_process( COMMAND date "+%Y-%m-%d"
                     OUTPUT_VARIABLE today_date
                     OUTPUT_STRIP_TRAILING_WHITESPACE )

  else( WIN32 )
    message( SEND_ERROR "Error: Failed to obtain today date - ${RESULT_VAR} set to 0000-00-00" )
    set( today_date "0000-00-00" )
  endif( WIN32 )
  set( ${RESULT_VAR} ${today_date} PARENT_SCOPE )
endfunction( TubeTKFunctionToday )
