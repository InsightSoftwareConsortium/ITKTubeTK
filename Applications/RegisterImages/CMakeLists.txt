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

set( MODULE_NAME RegisterImages )
project( ${MODULE_NAME} )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

option(RegisterImages_Support_2D_Images "Support of 2D Images in RegisterImages" OFF)

find_package( ITK REQUIRED )
if( TubeTK_BUILD_WITHIN_SLICER )
  set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 )
endif( TubeTK_BUILD_WITHIN_SLICER )
include( ${ITK_USE_FILE} )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES} ITKOptimizers ITKIOTransformBase
    TubeTKRegistration
  INSTALL_RUNTIME_DESTINATION bin )

if(RegisterImages_Support_2D_Images)
  target_compile_definitions(${MODULE_NAME}Lib PUBLIC SUPPORT_2D_IMAGES=1 )
endif()

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
