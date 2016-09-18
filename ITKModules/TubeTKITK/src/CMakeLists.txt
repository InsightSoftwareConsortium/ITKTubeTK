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
set( TubeTKITK_SRC
  itkTubeTK_fake_library_src.cxx
  )

# Include depends file to set TubeTKITK_LIBRARIES in current scope
include( ../tubetk-depends.cmake )

add_library( TubeTKITK STATIC ${TubeTKITK_SRC} )
target_link_libraries( TubeTKITK PUBLIC ${TubeTKITK_LIBRARIES} )

itk_module_link_dependencies()
itk_module_target( TubeTKITK )
