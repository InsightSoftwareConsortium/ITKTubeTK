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

set( MODULE_NAME RegisterUsingSlidingGeometries )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/docs/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
    ITKIOMeta
    ITKSpatialObjects
    TubeTK
    TubeTKLibCommon
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
