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

  if len(sys.argv) != 6:
    print("Usage: %s InputImage ScaleImage PointsImage OutputListFile threshold"%sys.argv[0])
    return 1
  inImage=sys.argv[1]
  inScale=sys.argv[2]
  inPoint=sys.argv[3]
  outputList=sys.argv[4]
  threshold=float(sys.argv[5])

  inImageReader=itk.ImageFileReader.New(FileName=inImage)
  inScaleReader=itk.ImageFileReader.New(FileName=inScale)
  inPointReader=itk.ImageFileReader.New(FileName=inPoint)

  convertFilter=TubeTKITK.ConvertShrunkenSeedImageToList[inScaleReader.GetOutput(), inPointReader.GetOutput()].New()

  castFilter = itk.CastImageFilter[inImageReader.GetOutput(),inScaleReader.GetOutput()].New()
  castFilter.SetInput(inImageReader)

  convertFilter.SetInput(castFilter.GetOutput())
  convertFilter.SetScaleImage(inScaleReader.GetOutput())

  convertFilter.SetPointsImage(inPointReader.GetOutput())
  convertFilter.SetThreshold(threshold)

  convertFilter.Update()
  matrix=convertFilter.GetOutput()

  narray=numpy.zeros((matrix.rows(),matrix.cols()))

  for ii in range(0,matrix.rows()):
      for jj in range(0,matrix.cols()):
        narray[ii,jj]=matrix.get(ii,jj)

  with open(outputList,'wb') as f :
    numpy.savetxt(f, narray, fmt='%.6g', delimiter=",")

if __name__ == "__main__":
  sys.exit(main())
