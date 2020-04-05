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
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( TubeTK_Common_H_Files
  Common/tubeIndent.h
  Common/tubeMacro.h
  Common/tubeMessage.h
  Common/tubeObject.h
  Common/tubeStringUtilities.h
  Common/tubeTestMain.h )

set( TubeTK_Common_HXX_Files )

set( TubeTK_Common_CXX_Files
  Common/tubeIndent.cxx
  Common/tubeObject.cxx )

list( APPEND TubeTK_SRCS
  ${TubeTK_Common_H_Files}
  ${TubeTK_Common_HXX_Files}
  ${TubeTK_Common_CXX_Files} )
