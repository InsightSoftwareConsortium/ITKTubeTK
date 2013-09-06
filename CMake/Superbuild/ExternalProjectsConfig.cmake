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

# See https://github.com/TubeTK/TubeTK/wiki/Dependencies

# Cppcheck version 1.61 (minimum version 1.52)
set( Cppcheck_URL https://github.com/danmar/cppcheck/archive/1.61.tar.gz )
set( Cppcheck_URL_MD5_TAG 362b6f8cf907262a7e3869b99183d734 )

# Common Toolkit snapshot 2013-08-20 16:21:07
set( CTK_URL ${github_protocol}://github.com/commontk/CTK.git )
set( CTK_URL_MD5_TAG 76219ba324ce1a976256fb7c7da5705d7e5f214b )

# TubeTK Image Viewer snapshot 2013-08-21 06:23:25
set( ImageViewer_URL
  ${github_protocol}://github.com/TubeTK/TubeTK-ImageViewer.git )
set( ImageViewer_URL_MD5_TAG 44adad4d8cc00e33b7c478c4bd0f4eb901bd7814 )

# Insight Segmentation and Registration Toolkit version 4.4.2
set( ITK_URL https://github.com/Kitware/ITK/archive/v4.4.2.tar.gz )
set( ITK_URL_MD5_TAG aa9006b3dcf05d615988c67aef241f4b )

# JsonCpp snapshot 2013-08-08 23:08:28
set( JsonCpp_URL ${svn_protocol}://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp )
set( JsonCpp_URL_MD5_TAG 274 )

# KWStyle snapshot 2012-04-19 04:05:19 (minimum version 1.0.0)
set( KWStyle_URL ${git_protocol}://public.kitware.com/KWStyle.git )
set( KWStyle_URL_MD5_TAG 16c5ca21e8133e6db155795dfdcb7d4bfa944af7 )

# LIBSVM version 3.17 (minimum version 3.1)
set( LIBSVM_URL
  http://www.csie.ntu.edu.tw/~cjlin/libsvm/oldfiles/libsvm-3.17.tar.gz )
set( LIBSVM_URL_MD5_TAG 67f8b597ce85c1f5288d7838e57ea28a )

# TubeTK Parameter Serializer snapshot 2013-08-21 11:14:41
set( ParameterSerializer_URL
  ${github_protocol}://github.com/TubeTK/TubeTK-ParameterSerializer.git )
set( ParameterSerializer_URL_MD5_TAG 819e3f6df97ab385f8c7b09203c07aa2fa56ade0 )

# Slicer Execution Model (TubeTK fork) snapshot 2013-08-05 10:16:50
set( SlicerExecutionModel_URL
  ${github_protocol}://github.com/TubeTK/TubeTK-SlicerExecutionModel.git )
set( SlicerExecutionModel_URL_MD5_TAG 5ba7ca1b8ed954e83c82fce604c8501a089d24e8 )

# Visualization Toolkit (3D Slicer fork) snapshot 2013-08-20 06:54:45
set( VTK_URL ${github_protocol}://github.com/Slicer/VTK.git )
set( VTK_URL_MD5_TAG 6da2c7235b8de9d8663a62f3246198d6865ecfae )
