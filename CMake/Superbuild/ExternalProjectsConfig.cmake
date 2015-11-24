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


# Cppcheck version 1.61
set( Cppcheck_URL ${github_protocol}://github.com/danmar/cppcheck.git )
set( Cppcheck_HASH_OR_TAG 1.61 )

# TubeTK Image Viewer
set( ImageViewer_URL
  ${github_protocol}://github.com/KitwareMedical/ImageViewer.git )
set( ImageViewer_HASH_OR_TAG 2f728d05fe159828e6221e1420b9c3a00295315b )

# JsonCpp
# ${svn_protocol}://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp )
# http://midas3.kitware.com/midas/download/bitstream/366544/JsonCpp_r276.tar.gz
set( JsonCpp_URL ${git_protocol}://github.com/KitwareMedical/jsoncpp.git )
set( JsonCpp_HASH_OR_TAG 34019a49daa7b7ee0692b6eeaf16df7dbdc3422a )

# KWStyle
set( KWStyle_URL
  ${git_protocol}://github.com/Kitware/KWStyle.git )
set( KWStyle_HASH_OR_TAG 22f2e0bd4ec4c94ca85af0ce24747e67a17b4c52 )

# LIBSVM version 3.2 (minimum version 3.1)
set( LIBSVM_URL
  ${git_protocol}://github.com/KitwareMedical/libsvm.git )
set( LIBSVM_HASH_OR_TAG 55fbd26a0bb53efffe89eb65b3a046c3ecc9965f )


###########################################################
###########################################################
# The following were copied from Slicer on 11/22/2015
###########################################################
###########################################################

# Common Toolkit
set( CTK_URL ${github_protocol}://github.com/commontk/CTK.git )
set( CTK_HASH_OR_TAG 5f581323b6526caf5ff293d9a810a3c9f68343d9 )

# Insight Segmentation and Registration Toolkit
set( ITK_URL ${github_protocol}://github.com/Slicer/ITK.git )
set( ITK_HASH_OR_TAG 03055046e2693349937339975df4e16b598609c1 )

# Slicer Execution Model
set( SlicerExecutionModel_URL
  ${github_protocol}://github.com/Slicer/SlicerExecutionModel.git )
set( SlicerExecutionModel_HASH_OR_TAG 311eff9038be165c1341dab90196fbd7466e715a )

# Visualization Toolkit (3D Slicer fork)
set( VTK_URL ${github_protocol}://github.com/Slicer/VTK.git )
set( VTK_HASH_OR_TAG fe92273888219edca422f3a308761ddcd2882e2b)
