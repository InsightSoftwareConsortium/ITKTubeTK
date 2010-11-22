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

project(Registration)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if(TubeTK_USE_VTK)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif(TubeTK_USE_VTK)
 
#set(Registration_SRCS
#  TubeRegistrator.cpp
#  )

#add_library(${PROJECT_NAME} ${Registration_SRCS})
#set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(${PROJECT_NAME} ${ITK_LIBRARIES})

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

