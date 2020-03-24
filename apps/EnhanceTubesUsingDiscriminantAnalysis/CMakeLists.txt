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

set( MODULE_NAME EnhanceTubesUsingDiscriminantAnalysis )

set( ${MODULE_NAME}_EXTERNAL_LIBS )
if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
  set( ${MODULE_NAME}_EXTERNAL_LIBS ${LIBSVM_LIBRARIES} )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_RANDOMFOREST )
  find_package( RandomForest REQUIRED )
  include_directories( ${RandomForest_INCLUDE_DIR} )
endif( TubeTK_USE_RANDOMFOREST )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/docs/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    TubeTKLibCommon
    TubeTKLibFiltering
    TubeTKLibMetaIO
    ${${MODULE_NAME}_EXTERNAL_LIBS}
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
