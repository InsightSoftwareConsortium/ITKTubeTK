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
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( TubeTK_MetaIO_H_Files
  MetaIO/itktubeMetaLDA.h
  MetaIO/itktubeMetaNJetLDA.h
  MetaIO/itktubeMetaClassPDF.h
  MetaIO/itktubeMetaRidgeSeed.h
  MetaIO/itktubeMetaTubeExtractor.h )

set( TubeTK_MetaIO_HXX_Files )

set( TubeTK_MetaIO_CXX_Files
  MetaIO/itktubeMetaClassPDF.cxx
  MetaIO/itktubeMetaLDA.cxx
  MetaIO/itktubeMetaNJetLDA.cxx
  MetaIO/itktubeMetaRidgeSeed.cxx
  MetaIO/itktubeMetaTubeExtractor.cxx )

list( APPEND TubeTK_SRCS
  ${TubeTK_MetaIO_H_Files}
  ${TubeTK_MetaIO_HXX_Files}
  ${TubeTK_MetaIO_CXX_Files} )
