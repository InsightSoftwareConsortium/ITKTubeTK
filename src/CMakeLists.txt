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
if( TubeTK_USE_VTK )
  find_package(VTK REQUIRED)
  include( VTK.cmake )
endif( TubeTK_USE_VTK )

itk_module_add_library( TubeTK ${TubeTK_SRCS} )
if( TubeTK_USE_VTK )
  target_link_libraries( TubeTK PUBLIC ${VTK_LIBRARIES} )
endif( TubeTK_USE_VTK )

if( ITK_WRAP_PYTHON )
  add_subdirectory( Python )
endif()

if( TubeTK_BUILD_APPLICATIONS )
  add_subdirectory( Applications )
endif()
