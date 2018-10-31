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

# See https://github.com/KitwareMedical/TubeTK/wiki/Dependencies

# Sanity checks
set(expected_nonempty_vars github_protocol git_protocol)
foreach(varname ${expected_nonempty_vars})
  if("${${varname}}" STREQUAL "")
    message(FATAL_ERROR "Variable '${varname}' is empty")
  endif()
endforeach()


###########################################################
# KitwareMedical Modules
###########################################################

# Cppcheck
set( Cppcheck_URL ${github_protocol}://github.com/KitwareMedical/cppcheck.git )
set( Cppcheck_HASH_OR_TAG bb0c8cbe7667615bca37ea2a02e2652642a5112e )

# JsonCpp
set( JsonCpp_URL ${git_protocol}://github.com/KitwareMedical/jsoncpp.git )
set( JsonCpp_HASH_OR_TAG 110d054227e9eb63faad48a1fb6a828ad0670e61 )

# RandomForest
set( RandomForest_URL
  ${git_protocol}://github.com/KitwareMedical/random-forest.git )
set( RandomForest_HASH_OR_TAG 6e1d0e271fea967487655555b8f26915aa1004d4 )

# LIBSVM
set( LIBSVM_URL ${git_protocol}://github.com/KitwareMedical/libsvm.git )
set( LIBSVM_HASH_OR_TAG 9bc3630f0f15fed7a5119c228c4d260574b4b6b2 )


###########################################################
# Kitware Modules
###########################################################

# Slicer Execution Model
set( SlicerExecutionModel_URL
  ${github_protocol}://github.com/Slicer/SlicerExecutionModel.git )
set( SlicerExecutionModel_HASH_OR_TAG
  fa5f22e6346f8115381c1a242c4c50afe0ee8a50 )

# KWStyle
set( KWStyle_URL ${git_protocol}://github.com/Kitware/KWStyle.git )
set( KWStyle_HASH_OR_TAG e03980ff514d5248a9f95ea355dcd9eff78c62d3 )


###########################################################
# ITK Modules
###########################################################

# TubeTKITK: Source already available in TubeTK project
set( TubeTKITK_URL ${TubeTK_SOURCE_DIR}/ITKModules/TubeTKITK )
set( TubeTKITK_HASH_OR_TAG "")

# MinimalPathExtraction
set( MinimalPathExtraction_URL
  ${git_protocol}://github.com/InsightSoftwareConsortium/ITKMinimalPathExtraction.git )
set( MinimalPathExtraction_HASH_OR_TAG
  21b01c83ffe9717e280020edb7f73996341945e9 )

set( TubeTK_ITK_MODULES
  TubeTKITK
  MinimalPathExtraction
  )
