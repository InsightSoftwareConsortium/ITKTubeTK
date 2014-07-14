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

# Common Toolkit snapshot 2014-01-03
set( CTK_URL ${github_protocol}://github.com/commontk/CTK.git )
set( CTK_HASH_OR_TAG 282abc745a7db700a822822375d0b8c170ec1a56 )

# TubeTK Image Viewer snapshot 2014-01-26
set( ImageViewer_URL ${github_protocol}://github.com/KitwareMedical/ImageViewer.git )
set( ImageViewer_HASH_OR_TAG d6d692d1a2a444cabbfa41a392d1c30143fd34bc )

# Insight Segmentation and Registration Toolkit
#set( ITK_URL ${github_protocol}://github.com/Kitware/ITK.git )
#set( ITK_HASH_OR_TAG 9fb84edbf9f6a5bb0846acbf6ecd767ca5cab2e1 )
set( ITK_URL ${github_protocol}://github.com/Slicer/ITK.git )
set( ITK_HASH_OR_TAG 3307a27de52d921ab40573e77097bfc72d70f5d6 )

# JsonCpp snapshot 2014-04-15 r276
# ${svn_protocol}://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp )
set( JsonCpp_URL http://midas3.kitware.com/midas/download/bitstream/366544/JsonCpp_r276.tar.gz )
set( JsonCpp_HASH_OR_TAG 192f0cf2b00798d4f4fb29c99a3aa83c )

# KWStyle snapshot 2012-04-19 04:05:19
set( KWStyle_URL
  ${git_protocol}://github.com/Kitware/KWStyle.git )
set( KWStyle_HASH_OR_TAG 9711cdbd35af37a9abcdd8b1dbd8e2b5a4ac8779 )

# LIBSVM version 3.17 (minimum version 3.1)
set( LIBSVM_URL
  http://www.csie.ntu.edu.tw/~cjlin/libsvm/oldfiles/libsvm-3.17.tar.gz )
set( LIBSVM_HASH_OR_TAG 67f8b597ce85c1f5288d7838e57ea28a )

# Parameter Serializer snapshot 2014-05-31 )
set( ParameterSerializer_URL
  ${github_protocol}://github.com/Slicer/ParameterSerializer.git )
set( ParameterSerializer_HASH_OR_TAG d7e94d3ff9de9a8782db6997f8e1bc45626cfa3b )

# Slicer Execution Model snapshot 2014-05-13
set( SlicerExecutionModel_URL
  ${github_protocol}://github.com/Slicer/SlicerExecutionModel.git )
set( SlicerExecutionModel_HASH_OR_TAG ecb58ad73c5dcdebf8c4d46257c7e44a07a7cdd4 )

# Visualization Toolkit (3D Slicer fork) snapshot 2013-08-20 06:54:45
set( VTK_URL ${github_protocol}://github.com/Slicer/VTK.git )
set( VTK_HASH_OR_TAG d8540e7ee356bbc025cb9917a41b7c7fa0548d4b )
