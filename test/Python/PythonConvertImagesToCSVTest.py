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
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################
#!/usr/bin/env python

import os
import sys

def GetRequiredEnvironmentVariable( varName ):
    if varName in os.environ:
        return os.environ[ varName ]
    else:
        print( '%s not found!' )%varName
        print( '  Set environment variable' )
        sys.exit( 1 )

def CheckIfPathExists( path, name ):
    if not os.path.exists( path ):
        print( 'Directory not found!' )
        print( '  %s = %s' )%( name, path )
        sys.exit( 1 )

def AppendSysPath( path ):
    BUILD_TYPE = GetRequiredEnvironmentVariable( 'BUILD_TYPE' )
    # Append path libs
    sys.path.append( os.path.join( os.path.join( path, \
                               'Wrapping/Generators/Python' ), BUILD_TYPE ) )
    # Folder containing *py files (and *a/*so files on Linux)
    sys.path.append( os.path.join( path, 'lib') )
    # Folder containing *lib files on Windows
    sys.path.append( os.path.join( os.path.join( path, \
                                   'lib' ), BUILD_TYPE) )
    # Windows needs this to load the DLL's
    os.environ[ 'PATH' ] += os.pathsep \
                         + os.path.join( os.path.join( path, 'bin' ),\
                         BUILD_TYPE )

# Path for ITK
ITK_BUILD_DIR = GetRequiredEnvironmentVariable( 'ITK_BUILD_DIR' )
CheckIfPathExists( ITK_BUILD_DIR, 'ITK_BUILD_DIR' )
AppendSysPath( ITK_BUILD_DIR )
# Path for TubeTK libs
TubeTK_BUILD_DIR = GetRequiredEnvironmentVariable( 'TubeTK_BUILD_DIR' )
CheckIfPathExists( TubeTK_BUILD_DIR, 'TubeTK_BUILD_DIR' )
AppendSysPath( TubeTK_BUILD_DIR )

import itk
from itk import TubeTKITK
import numpy

def main():

  if len(sys.argv) != 5:
    print("Usage: %s InputImage InputImageList OutputCSVFile stride"%sys.argv[0])
    return 1
  inputImage=sys.argv[1]
  inputImageFileNameList=sys.argv[2]
  outputCSVFile=sys.argv[3]
  stride=int(sys.argv[4])

  readerMask=itk.ImageFileReader.New(FileName=inputImage)
  readerMask.Update()

  imageFileNameList = inputImageFileNameList.split(',')
  reader1=itk.ImageFileReader.New(FileName=imageFileNameList[1])
  reader1.Update()

  convertFilter=TubeTKITK.ConvertImagesToCSV[reader1.GetOutput(), readerMask.GetOutput()].New()
  convertFilter.SetInputMask(readerMask.GetOutput())

  numImages=0
  for image in imageFileNameList:
    reader=itk.ImageFileReader.New(FileName=image)
    reader.Update()
    # Cast
    castFilter = itk.CastImageFilter[reader.GetOutput(),reader1.GetOutput()].New()
    castFilter.SetInput(reader)
    castFilter.Update()
    # Add image
    convertFilter.AddImage(castFilter.GetOutput())
    numImages += 1

  convertFilter.SetStride(stride)
  convertFilter.SetNumImages(numImages)

  convertFilter.Update()
  numberRows=convertFilter.GetNumberRows()
  matrix=convertFilter.GetOutput()

  narray=numpy.zeros((numberRows+1,matrix.cols()))

  for ii in range(0,numberRows):
      for jj in range(0,matrix.cols()):
        narray[ii,jj]=matrix.get(ii,jj)

  with open(outputCSVFile,'wb') as f :
    for jj in range (0,matrix.cols()-1) :
      f.write(os.path.basename(imageFileNameList[jj])+',')
    f.write('Class\n')
    numpy.savetxt(f, narray, fmt='%.6g', delimiter=",")

if __name__ == "__main__":
  sys.exit(main())
