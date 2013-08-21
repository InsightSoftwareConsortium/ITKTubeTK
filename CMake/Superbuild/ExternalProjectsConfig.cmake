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

# Cppcheck version 1.61
set( Cppcheck_GIT_REPOSITORY
  ${git_protocol}://github.com/danmar/cppcheck.git )
set( Cppcheck_GIT_TAG 9e16abeb4275f01b2041972a2edbd6ad46d945c1 )

# Common Toolkit snapshot 2013-08-20 16:21:07
set( CTK_GIT_REPOSITORY ${git_protocol}://github.com/commontk/CTK.git )
set( CTK_GIT_TAG 76219ba324ce1a976256fb7c7da5705d7e5f214b )

# TubeTK Image Viewer snapshot 2013-08-21 06:23:25
set( ImageViewer_GIT_REPOSITORY
  ${git_protocol}://github.com/TubeTK/TubeTK-ImageViewer.git )
set( ImageViewer_GIT_TAG 44adad4d8cc00e33b7c478c4bd0f4eb901bd7814 )

# Insight Segmentation and Registration Toolkit snapshot 2013-07-12 12:08:29
set( ITK_GIT_REPOSITORY ${git_protocol}://github.com/Kitware/ITK.git )
set( ITK_GIT_TAG 35b90133a793ffd884820e499175db19366fe627 )

# JsonCpp snapshot 2013-08-08 23:08:28
set( JsonCpp_SVN_REPOSITORY
  ${svn_protocol}://svn.code.sf.net/p/jsoncpp/code/trunk/jsoncpp )
set( JsonCpp_SVN_REVISION 274 )

# KWStyle snapshot 2012-04-19 04:05:19
set( KWStyle_GIT_REPOSITORY
  ${git_protocol}://public.kitware.com/KWStyle.git )
set( KWStyle_GIT_TAG 16c5ca21e8133e6db155795dfdcb7d4bfa944af7 )

# LIBSVM version 3.17 (minimum version 3.1)
set( LIBSVM_GIT_REPOSITORY
  ${git_protocol}://github.com/TubeTK/TubeTK-LIBSVM.git )
set( LIBSVM_GIT_TAG 74907dab39c105f88dcc3f0e14c6f5562b5affb4 )

# TubeTK Parameter Serializer snapshot 2013-08-21 11:14:41
set( ParameterSerializer_GIT_REPOSITORY
  ${git_protocol}://github.com/TubeTK/TubeTK-ParameterSerializer.git )
set( ParameterSerializer_GIT_TAG 819e3f6df97ab385f8c7b09203c07aa2fa56ade0 )

# Slicer Execution Model (TubeTK fork) snapshot 2013-08-05 10:16:50
set( SlicerExecutionModel_GIT_REPOSITORY
  ${git_protocol}://github.com/TubeTK/TubeTK-SlicerExecutionModel.git )
set( SlicerExecutionModel_GIT_TAG 5ba7ca1b8ed954e83c82fce604c8501a089d24e8 )

# Visualization Toolkit (3D Slicer fork) snapshot 2013-08-20 06:54:45
set( VTK_GIT_REPOSITORY ${git_protocol}://github.com/Slicer/VTK.git )
set( VTK_GIT_TAG 6da2c7235b8de9d8663a62f3246198d6865ecfae )
