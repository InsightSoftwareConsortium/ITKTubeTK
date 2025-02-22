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
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( TubeTK_SRCS "" )

include( Common.cmake )
include( Filtering.cmake )
include( IO.cmake )
include( MetaIO.cmake )
include( Numerics.cmake )
include( ObjectDocuments.cmake )
include( Registration.cmake )
include( Segmentation.cmake )

add_library(TubeTK STATIC ${TubeTK_SRCS})
itk_module_link_dependencies()
itk_module_target(TubeTK)
