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

def main():
  if len( sys.argv ) != 6:
    print( "Usage: %s InputImage OutputVesselMask OutputNotVesselMask gap\
            notVesselWidth"%sys.argv[0] )
    return 1
  inputImage = sys.argv[1]
  outputVesselMask = sys.argv[2]
  outputNotVesselMask = sys.argv[3]
  gap=float( sys.argv[4] )
  notVesselWidth=float( sys.argv[5] )

  reader=itk.ImageFileReader.New( FileName=inputImage )
  reader.Update()
  trainingMask=TubeTKITK.ComputeTrainingMask.New( reader )
  trainingMask.SetGap(gap)
  trainingMask.SetNotVesselWidth( notVesselWidth )
  trainingMask.Update()
  writer=itk.ImageFileWriter.New( trainingMask, FileName=outputVesselMask )
  writer.Update()

  writer=itk.ImageFileWriter.New( trainingMask.GetNotVesselMask(),\
                                  FileName=outputNotVesselMask )
  writer.Update()


if __name__ == "__main__":
  sys.exit( main() )
