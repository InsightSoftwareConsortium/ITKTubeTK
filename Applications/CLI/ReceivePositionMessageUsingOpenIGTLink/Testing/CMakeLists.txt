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

include_regular_expression( "^.*$" )

set( TEST_EXE
 SendAndReceivePositionMessageUsingOpenIGTLinkTest )

set( TEST_SOURCE 
 SendAndReceivePositionMessageUsingOpenIGTLinkTest.cxx
  )

set( TEST_PATH 
 ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR} )

#######
# Test1
#######

generateclp( TEST_SOURCE SendAndReceivePositionMessageUsingOpenIGTLinkTest.xml )

add_executable( ${TEST_EXE}
                ${TEST_SOURCE} )
target_link_libraries( ${TEST_EXE} 
     ReceivePositionMessageUsingOpenIGTLinkModule 
     SendPositionMessageUsingOpenIGTLinkModule )

slicer3_set_plugins_output_path( ${TEST_EXE} )
set( TARGETS 
     ${TEST_EXE} )
slicer3_install_plugins( ${TARGETS} )

add_test( ${TEST_EXE}1  
     ${TEST_PATH}/${TEST_EXE} 
     ${TEST_PATH}/SendPositionMessageUsingOpenIGTLink
     ${TEST_PATH}/ReceivePositionMessageUsingOpenIGTLink
     18000 
     localhost)
