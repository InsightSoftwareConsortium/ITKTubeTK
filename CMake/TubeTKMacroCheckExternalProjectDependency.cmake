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

##############################################################################
#
#  Program: 3D Slicer
#
#  Copyright (c) Kitware Inc.
#
#  See COPYRIGHT.txt
#  or http://www.slicer.org/copyright/copyright.txt for details.
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#  This file was originally developed by Jean-Christophe Fillion-Robin,
#  Kitware Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
#
##############################################################################

if( NOT EXISTS "${EXTERNAL_PROJECT_DIR}" )
  set( EXTERNAL_PROJECT_DIR ${${CMAKE_PROJECT_NAME}_SOURCE_DIR}/CMake/Superbuild )
endif( NOT EXISTS "${EXTERNAL_PROJECT_DIR}" )

macro( TubeTKMacroCheckExternalProjectDependency proj )

  if( NOT DEFINED __indent )
    set( __indent "" )
  else( NOT DEFINED __indent )
    set( __indent "${__indent}  " )
  endif( NOT DEFINED __indent )

  # Sanity checks
  if( NOT DEFINED ${proj}_DEPENDENCIES )
    message( FATAL_ERROR "${__indent}${proj}_DEPENDENCIES variable is NOT defined" )
  endif( NOT DEFINED ${proj}_DEPENDENCIES )

  # Display dependency of project being processed
  if( "${${proj}_DEPENDENCIES}" STREQUAL "" )
    message( STATUS "SuperBuild - ${__indent}${proj}[OK]" )
  else( "${${proj}_DEPENDENCIES}" STREQUAL "" )
    set( dependency_str " " )
    foreach( dep ${${proj}_DEPENDENCIES} )
      if( External_${dep}_FILE_INCLUDED )
        set( dependency_str "${dependency_str}${dep}[INCLUDED], " )
      else( External_${dep}_FILE_INCLUDED )
        set( dependency_str "${dependency_str}${dep}, " )
      endif( External_${dep}_FILE_INCLUDED )
    endforeach()
    message( STATUS "SuperBuild - ${__indent}${proj} => Requires${dependency_str}" )
  endif( "${${proj}_DEPENDENCIES}" STREQUAL "" )

  foreach( dep ${${proj}_DEPENDENCIES} )
    if( NOT External_${dep}_FILE_INCLUDED )
      set( external_project_found OFF )
      if( EXISTS "${EXTERNAL_PROJECT_DIR}/External_${dep}.cmake" )
        include( ${EXTERNAL_PROJECT_DIR}/External_${dep}.cmake )
        set( external_project_found ON )
      endif( EXISTS "${EXTERNAL_PROJECT_DIR}/External_${dep}.cmake" )

      if( EXISTS "${Slicer_ADDITIONAL_EXTERNAL_PROJECT_DIR}/External_${dep}.cmake"
        AND NOT ${external_project_found} )
        include( ${${proj}_ADDITIONAL_EXTERNAL_PROJECT_DIR}/External_${dep}.cmake )
        set( external_project_found ON )
      endif( EXISTS "${Slicer_ADDITIONAL_EXTERNAL_PROJECT_DIR}/External_${dep}.cmake"
        AND NOT ${external_project_found} )

      if( NOT ${external_project_found} )
        foreach( dir ${${proj}_ADDITIONAL_EXTERNAL_PROJECT_DIRS} )
          if( EXISTS "${dir}/External_${dep}.cmake" AND NOT ${external_project_found} )
            include( ${${proj}_ADDITIONAL_EXTERNAL_PROJECT_DIR}/External_${dep}.cmake )
            set( external_project_found ON )
          endif( EXISTS "${dir}/External_${dep}.cmake" AND NOT ${external_project_found} )
        endforeach()
      endif( NOT ${external_project_found} )

      if( NOT ${external_project_found} )
        message( FATAL_ERROR "Cannot find External_${dep}.cmake" )
      endif( NOT ${external_project_found} )
    endif( NOT External_${dep}_FILE_INCLUDED )
  endforeach()

  set( __${proj}_superbuild_message "SuperBuild - ${__indent}${proj}[OK]" )
  set( __${proj}_indent ${__indent} )

  # If project being process has dependencies, indicates it has also been added.
  if( NOT "${${proj}_DEPENDENCIES}" STREQUAL "" )
    message( STATUS ${__${proj}_superbuild_message} )
  endif( NOT "${${proj}_DEPENDENCIES}" STREQUAL "" )


  # Update indent variable
  string( LENGTH "${__indent}" __indent_length )
  math( EXPR __indent_length "${__indent_length}-2" )
  if( NOT ${__indent_length} LESS 0 )
    string( SUBSTRING "${__indent}" 0 ${__indent_length} __indent )
  endif( NOT ${__indent_length} LESS 0 )
endmacro( TubeTKMacroCheckExternalProjectDependency proj )
