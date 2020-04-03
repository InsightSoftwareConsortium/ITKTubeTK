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

set( TubeTK_IO_H_Files
  IO/itktubePDFSegmenterParzenIO.h
  IO/itktubeRidgeSeedFilterIO.h
  IO/itktubeTubeExtractorIO.h
  IO/itktubeTubeXIO.h )

if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTK_IO_H_Files
    IO/itktubePDFSegmenterSVMIO.h )
endif()

if( TubeTK_USE_RANDOMFOREST )
  list( APPEND TubeTK_IO_H_Files
    IO/itktubePDFSegmenterRandomForestIO.h )
endif()

set( TubeTK_IO_HXX_Files
  IO/itktubePDFSegmenterParzenIO.hxx
  IO/itktubeRidgeSeedFilterIO.hxx
  IO/itktubeTubeExtractorIO.hxx
  IO/itktubeTubeXIO.hxx )

if( TubeTK_USE_LIBSVM )
  list( APPEND TubeTK_IO_HXX_Files
    IO/itktubePDFSegmenterSVMIO.hxx )
endif()

if( TubeTK_USE_RANDOMFOREST )
  list( APPEND TubeTK_IO_HXX_Files
    IO/itktubePDFSegmenterRandomForestIO.hxx )
endif()

set( TubeTK_IO_CXX_Files )

list( APPEND TubeTK_SRCS
  ${TubeTK_IO_H_Files}
  ${TubeTK_IO_HXX_Files}
  ${TubeTK_IO_CXX_Files} )
