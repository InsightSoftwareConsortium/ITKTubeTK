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
set( Cppcheck_URL
  ${github_protocol}://github.com/danmar/cppcheck.git )
set( Cppcheck_HASH_OR_TAG 1.61 )

# TubeTK Image Viewer
set( ImageViewer_URL ${github_protocol}://github.com/KitwareMedical/ImageViewer.git )
set( ImageViewer_HASH_OR_TAG 2f728d05fe159828e6221e1420b9c3a00295315b )

# JsonCpp
# ${svn_protocol}://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp )
set( JsonCpp_URL http://midas3.kitware.com/midas/download/bitstream/366544/JsonCpp_r276.tar.gz )
set( JsonCpp_HASH_OR_TAG 192f0cf2b00798d4f4fb29c99a3aa83c )

# KWStyle
set( KWStyle_URL
  ${git_protocol}://github.com/Kitware/KWStyle.git )
set( KWStyle_HASH_OR_TAG c502d2c92b76e925f709f7747353ff0ecb7200d8 )

# LIBSVM version 3.17 (minimum version 3.1)
set( LIBSVM_URL
  http://www.csie.ntu.edu.tw/~cjlin/libsvm/oldfiles/libsvm-3.17.tar.gz )
set( LIBSVM_HASH_OR_TAG 67f8b597ce85c1f5288d7838e57ea28a )


###########################################################
###########################################################
# The following were copied from Slicer on 2015.08.29
###########################################################
###########################################################

# Common Toolkit
set( CTK_URL ${github_protocol}://github.com/commontk/CTK.git )
set( CTK_HASH_OR_TAG 35a187ff5ae15f1c590a2ac317cd102c5746a81b )

# Insight Segmentation and Registration Toolkit
set( ITK_URL ${github_protocol}://github.com/Slicer/ITK.git )
set( ITK_HASH_OR_TAG 0905a43149320d8c2993e83ecf6c8a9d6b2b4232 )

# Slicer Execution Model
set( SlicerExecutionModel_URL
  ${github_protocol}://github.com/Slicer/SlicerExecutionModel.git )
set( SlicerExecutionModel_HASH_OR_TAG 608c7b06f9402b35f0ec01550a3b9e313582f683 )

# Visualization Toolkit (3D Slicer fork)
set( VTK_URL ${github_protocol}://github.com/Slicer/VTK.git )
set( VTK_HASH_OR_TAG d23ea3680dcdd375d6fdab50de9dae7090f73ed6 )
