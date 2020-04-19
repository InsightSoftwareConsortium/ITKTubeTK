##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
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

set( TubeTK_SRCS "" )

include( Common.cmake )
include( Filtering.cmake )
include( IO.cmake )
include( MetaIO.cmake )
include( Numerics.cmake )
include( ObjectDocuments.cmake )
include( Registration.cmake )
include( Segmentation.cmake )

if( TubeTK_USE_ARRAYFIRE )
  include(CheckLanguage)
  check_language(CUDA)
  enable_language(CUDA)
endif()

itk_module_add_library( TubeTK ${TubeTK_SRCS} )
target_compile_features( TubeTK PUBLIC cxx_std_11 )
set_target_properties( TubeTK PROPERTIES POSITION_INDEPENDENT_CODE ON )

if( TubeTK_USE_ARRAYFIRE )
  set_target_properties( TubeTK PROPERTIES CUDA_SEPERABLE_COMPILATION ON )
  install( TARGETS TubeTK OBJECTS DESTINATION bin/ptx )
endif( TubeTK_USE_ARRAYFIRE )

if( ITK_WRAP_PYTHON )
  add_subdirectory( Python )
endif()
