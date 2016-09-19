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

Midas3FunctionAddTest( NAME Python.ConvertTubesToImages
  COMMAND ${PYTHON_TESTING_EXECUTABLE}
    "${NOTEBOOK_TEST_DRIVER}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../ConvertTubesToImages.ipynb"
    "${TubeTK_BINARY_DIR}"
    MIDAS_FETCH_ONLY{Branch.n010.mha.md5}
    MIDAS_FETCH_ONLY{Branch-truth.tre.md5} )
