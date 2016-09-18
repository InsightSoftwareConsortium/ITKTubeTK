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

#
# TubeTK Maintainers:
#
# TubeTK base components required to build TubeTKITK should be added to
#   TubeTKITK_DEPENDS list.
#
# * Components should be listed in topological order.
#
# * Name should be specified without the 'TubeTK' prefix.
#
# * All components should be listed: This will ensure all required folder
#     will be added when building TubeTK as an ITK remote module.
#

set( TubeTKITK_DEPENDS
  Common
  Numerics
  Filtering
  Segmentation
  Registration )

set( TubeTKITK_LIBRARIES )
foreach( component ${TubeTKITK_DEPENDS} )
  list( APPEND TubeTKITK_LIBRARIES TubeTK${component} )
endforeach()
