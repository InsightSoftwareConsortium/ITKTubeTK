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
# Configure tubetk.pth.
# Add a symlink to this file in the directory that results from::
#
#   python -m site --user-site
#
# to make TubeTK Python packages and modules importable.
#
configure_file( tubetk.pth.in tubetk.pth )

add_subdirectory( tubetk )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
