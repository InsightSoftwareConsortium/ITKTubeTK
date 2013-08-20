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

find_program( KWSTYLE_EXECUTABLE NAMES KWStyle PATHS /usr/local/bin
  "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Kitware Inc.\\KWStyle 1.0.0]/bin" )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( KWStyle DEFAULT_MSG KWSTYLE_EXECUTABLE )

mark_as_advanced( KWSTYLE_EXECUTABLE )
