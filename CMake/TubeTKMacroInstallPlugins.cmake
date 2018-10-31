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

# Install one or more plugins to the default plugin location.
macro( TubeTKMacroInstallPlugins )
  set( _extra )
  if( ITK_SOURCE_DIR )
    set( _extra EXPORT ITKTargets )
    set( TubeTK_INSTALL_PLUGINS_BIN_DIR ${ITK_INSTALL_RUNTIME_DIR} )
    set( TubeTK_INSTALL_PLUGINS_LIB_DIR ${ITK_INSTALL_LIBRARY_DIR} )
  endif()
  install( TARGETS ${ARGN}
    ${_extra}
    RUNTIME DESTINATION ${TubeTK_INSTALL_PLUGINS_BIN_DIR} COMPONENT
	  RuntimeLibraries
    LIBRARY DESTINATION ${TubeTK_INSTALL_PLUGINS_LIB_DIR} COMPONENT
	  RuntimeLibraries
    ARCHIVE DESTINATION ${TubeTK_INSTALL_PLUGINS_LIB_DIR} COMPONENT
	  Development  )
endmacro( TubeTKMacroInstallPlugins )
