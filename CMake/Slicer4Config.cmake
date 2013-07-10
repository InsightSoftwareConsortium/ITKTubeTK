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

#  We use Slicer4 variable names to simplify porting to Slicer4
set( Slicer4_INSTALL_BIN_DIR "bin" )
set( Slicer4_INSTALL_LIB_DIR "lib/TubeTK" )
set( Slicer4_INSTALL_INCLUDE_DIR "include/TubeTK" )
set( Slicer4_INSTALL_SHARE_DIR "share/TubeTK" )
set( Slicer4_INSTALL_ITKFACTORIES_DIR "${Slicer4_INSTALL_LIB_DIR}/ITKFactories" )

set( Slicer4_INSTALL_PLUGINS_BIN_DIR "${Slicer4_INSTALL_BIN_DIR}" )
set( Slicer4_INSTALL_PLUGINS_LIB_DIR "${Slicer4_INSTALL_LIB_DIR}" )
set( Slicer4_INSTALL_PLUGINS_INCLUDE_DIR "${Slicer4_INSTALL_INCLUDE_DIR}" )
set( Slicer4_INSTALL_PLUGINS_SHARE_DIR "${Slicer4_INSTALL_SHARE_DIR}" )
set( Slicer4_INSTALL_PLUGINS_CACHE_DIR "${Slicer4_INSTALL_LIB_DIR}/CACHE" )

set( Slicer4_HOME "${Slicer4_BINARY_DIR}" )

include( ${TubeTK_SOURCE_DIR}/CMake/Slicer4PluginsMacros.cmake )
include( ${SlicerExecutionModel_DIR}/CMake/SEMMacroBuildCLI.cmake )

# Needed by SlicerBaseCLI
set( Slicer_INSTALL_BIN_DIR ${Slicer4_INSTALL_BIN_DIR} )
set( Slicer_INSTALL_LIB_DIR ${Slicer4_INSTALL_LIB_DIR} )
