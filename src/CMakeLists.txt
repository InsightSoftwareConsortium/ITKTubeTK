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

set( TubeTK_SRCS
  itkTubeTK_fake_library_src.cxx
  )

add_library( TubeTK ${TubeTK_SRCS} )

itk_module_link_dependencies()

add_dependencies( TubeTK ${TubeTKLib_LIBRARIES} )

target_link_libraries( TubeTK INTERFACE
  ${TubeTKLib_LIBRARIES}
  ${ITK_LIBRARIES} 
  )

target_include_directories( TubeTK INTERFACE
  ${TubeTKLib_INCLUDE_DIRS}
  )

itk_module_target( TubeTK )
