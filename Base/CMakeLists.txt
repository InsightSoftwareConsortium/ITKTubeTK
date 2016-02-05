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

add_subdirectory( CLI )
add_subdirectory( Common )
add_subdirectory( Numerics )
add_subdirectory( IO )
add_subdirectory( ObjectDocuments )

if( TubeTK_USE_PYTHON )
  add_subdirectory( Python )
endif( TubeTK_USE_PYTHON )

add_subdirectory( Filtering )
add_subdirectory( Registration )
add_subdirectory( Segmentation )
add_subdirectory( USTK )
