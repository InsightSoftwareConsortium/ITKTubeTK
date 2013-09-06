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

cmake_minimum_required( VERSION 2.8.7 )

project( SVM )

set( SVM_MAJOR_VERSION 3 )
set( SVM_MINOR_VERSION 17 )
set( SVM_PATCH_VERSION 0 )
set( SVM_VERSION ${SVM_MAJOR_VERSION}.${SVM_MINOR_VERSION}.${SVM_PATCH_VERSION} )
set( SVM_LIBRARY_SO_VERSION 2 )

set( SVM_TARGETS_NAME svm-targets )
set( SVM_INSTALL_BIN_DIR bin )
set( SVM_INSTALL_LIB_DIR lib )
set( SVM_INSTALL_INCLUDE_DIR include )
set( SVM_INSTALL_SHARE_DIR share/libsvm-${SVM_VERSION} )

set( EXECUTABLE_OUTPUT_PATH "${SVM_BINARY_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH "${SVM_BINARY_DIR}/lib" )

if( NOT CMAKE_BUILD_TYPE )
  message( STATUS "CMAKE_BUILD_TYPE was not specified, defaulting to Release" )
  set( CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration type" FORCE )
endif( NOT CMAKE_BUILD_TYPE )

option( BUILD_SHARED_LIBS "Build LIBSVM as a shared library." ON )

add_executable( svm-train svm-train.c svm.h svm.cpp )
add_executable( svm-predict svm-predict.c svm.h svm.cpp )
add_executable( svm-scale svm-scale.c svm.h svm.cpp )

set( SVM_MATH_LIBRARY m )
target_link_libraries( svm-train ${SVM_MATH_LIBRARY} )
target_link_libraries( svm-predict ${SVM_MATH_LIBRARY} )
target_link_libraries( svm-scale ${SVM_MATH_LIBRARY} )

add_library( svm svm.cpp )
target_link_libraries( svm ${SVM_MATH_LIBRARY} )
set_target_properties( svm PROPERTIES SOVERSION ${SVM_LIBRARY_SO_VERSION} VERSION ${SVM_VERSION} )

install( TARGETS svm svm-train svm-predict svm-scale
  EXPORT ${SVM_TARGETS_NAME}
  RUNTIME DESTINATION ${SVM_INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${SVM_INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${SVM_INSTALL_LIB_DIR} )

install( EXPORT ${SVM_TARGETS_NAME} DESTINATION ${SVM_INSTALL_SHARE_DIR} )

install( FILES svm.h DESTINATION ${SVM_INSTALL_INCLUDE_DIR} )
