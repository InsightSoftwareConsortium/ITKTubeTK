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

# Cppcheck
set( Cppcheck_URL ${github_protocol}://github.com/KitwareMedical/cppcheck.git )
set( Cppcheck_HASH_OR_TAG a9c9482d6e1b42457aedf8065e21523654f46124 )

# TubeTK Image Viewer
set( ImageViewer_URL
  ${github_protocol}://github.com/KitwareMedical/ImageViewer.git )
set( ImageViewer_HASH_OR_TAG 2123736abe01a4cee91c231011086e3770058088 )

# JsonCpp
# http://midas3.kitware.com/midas/download/bitstream/366544/JsonCpp_r276.tar.gz
set( JsonCpp_URL ${git_protocol}://github.com/KitwareMedical/jsoncpp.git )
set( JsonCpp_HASH_OR_TAG 110d054227e9eb63faad48a1fb6a828ad0670e61 )

# KWStyle
set( KWStyle_URL ${git_protocol}://github.com/Kitware/KWStyle.git )
set( KWStyle_HASH_OR_TAG d5b6d0208200d2898d2ed3fdf72b2bb2810ca4fa )

# LIBSVM
set( LIBSVM_URL ${git_protocol}://github.com/KitwareMedical/libsvm.git )
set( LIBSVM_HASH_OR_TAG 9bc3630f0f15fed7a5119c228c4d260574b4b6b2 )

# RandomForest
set( RandomForest_URL
  ${git_protocol}://github.com/KitwareMedical/random-forest.git )
set( RandomForest_HASH_OR_TAG 6e1d0e271fea967487655555b8f26915aa1004d4 )

###########################################################
# ITK Modules
###########################################################

# TubeTKITK: Source already available in TubeTK project
set( TubeTKITK_URL ${TubeTK_SOURCE_DIR}/ITKModules/TubeTKITK )
set( TubeTKITK_HASH_OR_TAG "")

# MinimalPathExtraction
set( MinimalPathExtraction_URL
  ${git_protocol}://github.com/InsightSoftwareConsortium/ITKMinimalPathExtraction.git )
set( MinimalPathExtraction_HASH_OR_TAG 5a2017ef5d5c25db518ecae49408598f906dd307 )

set( TubeTK_ITK_MODULES
  TubeTKITK
  MinimalPathExtraction
  )

###########################################################
###########################################################
# The following were copied from Slicer on 12/03/2016
###########################################################
###########################################################

# Common Toolkit
# Not updated because compilation error remains in CTK head:
# Creating library bin\CTKVisualizationVTKCore.lib and object
#  bin\CTKVisualizationVTKCore.exp
# vtkLightBoxRendererManager.cpp.obj : error LNK2019: unresolved external
#  symbol "void __cdecl vtkRenderingFreeType_AutoInit_Construct(void)"
#  (?vtkRenderingFreeType_AutoInit_Construct@@YAXXZ) referenced in function
#  "public: __cdecl vtkRenderingCore_AutoInit::vtkRenderingCore_AutoInit(void)"
#  (??0vtkRenderingCore_AutoInit@@QEAA@XZ)
# vtkLightBoxRendererManager.cpp.obj : error LNK2019: unresolved external
#  symbol "void __cdecl vtkRenderingFreeType_AutoInit_Destruct(void)"
#  (?vtkRenderingFreeType_AutoInit_Destruct@@YAXXZ) referenced in function
#  "public: __cdecl vtkRenderingCore_AutoInit::~vtkRenderingCore_AutoInit(void)"
#  (??1vtkRenderingCore_AutoInit@@QEAA@XZ)
# bin\CTKVisualizationVTKCore.dll : fatal error LNK1120: 2 unresolved externals
set( CTK_URL ${github_protocol}://github.com/aylward/CTK.git )
set( CTK_HASH_OR_TAG d1ebb429c952a329615f5fa9223a5897c6cc459c )

# Insight Segmentation and Registration Toolkit
set( ITK_URL ${github_protocol}://github.com/InsightSoftwareConsortium/ITK.git )
set( ITK_HASH_OR_TAG f535cf3cd87178b06eb682358a6d19fd7d5ed7da )

# Slicer Execution Model
set( SlicerExecutionModel_URL
  ${github_protocol}://github.com/Slicer/SlicerExecutionModel.git )
set( SlicerExecutionModel_HASH_OR_TAG 62d0121dbb0fb057ebbd7c9ab84520accacec8bc )

# Visualization Toolkit (3D Slicer fork)
set( VTK_URL ${github_protocol}://github.com/Slicer/VTK.git )
set( VTK_HASH_OR_TAG a024cefc2acf25350734e6f04d2562f9a6a3b124 )
