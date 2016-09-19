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
macro( _subdirectory_conditional_add name )
  if( NOT DEFINED TubeTK_BUILD_BASE_${name} )
    set( TubeTK_BUILD_BASE_${name} 1 )
    if( ITK_SOURCE_DIR )
      # TubeTK base libraries required by TubeTKITK module will have
      # to be explicitly enabled.
      set( TubeTK_BUILD_BASE_${name} 0 )
    endif()
  endif()
  if( TubeTK_BUILD_BASE_${name} )
    add_subdirectory( ${name} )
  endif()
endmacro()

if( ITK_SOURCE_DIR )
  include("../CMake/TubeTKMacroInstallPlugins.cmake")
endif()

_subdirectory_conditional_add( CLI )
_subdirectory_conditional_add( Common )
_subdirectory_conditional_add( CompareTools )

_subdirectory_conditional_add( MetaIO )

_subdirectory_conditional_add( Numerics )
_subdirectory_conditional_add( ObjectDocuments )

if( TubeTK_USE_PYTHON )
  _subdirectory_conditional_add( Python )
endif( TubeTK_USE_PYTHON )

_subdirectory_conditional_add( Filtering )
_subdirectory_conditional_add( Registration )
_subdirectory_conditional_add( Segmentation )
_subdirectory_conditional_add( IO )
_subdirectory_conditional_add( USTK )
