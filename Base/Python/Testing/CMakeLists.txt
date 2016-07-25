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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

if( ${TubeTK_USE_NUMPY_STACK} )
  Midas3FunctionAddTestWithEnv( NAME VesselTubeToNumPyTest
    COMMAND ${PYTHON_TESTING_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_tubetk.py
      VesselTubeToNumPyTest
        MIDAS{tube.tre.md5}
        MIDAS{tube.tre.npy.md5}
    ENVIRONMENT TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR}
    )
  if( ${TubeTK_USE_PYQTGRAPH} )
    Midas3FunctionAddTestWithEnv( NAME PyQtGraphTubesAsCirclesTest
      COMMAND ${PYTHON_TESTING_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_tubetk.py
        PyQtGraphTubesAsCirclesTest
          MIDAS{tube.tre.md5}
          ${TEMP}/PyQtGraphTubesAsCirclesTest.png
      ENVIRONMENT TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR}
      )
    Midas3FunctionAddTestWithEnv( NAME PyQtGraphTubesAsCirclesTest-Compare
      COMMAND $<TARGET_FILE:CompareImages>
        -t ${TEMP}/PyQtGraphTubesAsCirclesTest.png
        -b MIDAS{PyQtGraphTubesAsCirclesTest.png.md5},MIDAS{PyQtGraphTubesAsCirclesTest.1.png.md5}
        MIDAS_FETCH_ONLY{PyQtGraphTubesAsCirclesTest.png.md5}
        MIDAS_FETCH_ONLY{PyQtGraphTubesAsCirclesTest.1.png.md5}
      ENVIRONMENT TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR}
      )
    set_property( TEST
      PyQtGraphTubesAsCirclesTest-Compare
      APPEND PROPERTY DEPENDS
      PyQtGraphTubesAsCirclesTest )
  endif( ${TubeTK_USE_PYQTGRAPH} )
endif( ${TubeTK_USE_NUMPY_STACK} )
