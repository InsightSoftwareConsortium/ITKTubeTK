#!/usr/bin/env python

import os
import sys


# Path for TubeTK libs
TubeTK_BUILD_DIR=None
if 'TubeTK_BUILD_DIR' in os.environ:
    TubeTK_BUILD_DIR = os.environ['TubeTK_BUILD_DIR']
if not os.path.exists(TubeTK_BUILD_DIR):
    print('TubeTK_BUILD_DIR not found!')
    print('  Set environment variable')
    sys.exit(1)

try:
    import itk
except:
    ITK_BUILD_DIR = None
    if 'ITK_BUILD_DIR' in os.environ:
        ITK_BUILD_DIR = os.environ['ITK_BUILD_DIR']
    else:
        print('ITK_BUILD_DIR not found!')
        print('  Set environment variable')
        sys.exit( 1 )

    if not os.path.exists(ITK_BUILD_DIR):
        print('ITK_BUILD_DIR set by directory not found!')
        print('  ITK_BUIDL_DIR = ' + ITK_BUILD_DIR )
        sys.exit(1)
    # Append ITK libs
    sys.path.append(os.path.join(ITK_BUILD_DIR, 'Wrapping/Generators/Python'))
    sys.path.append(os.path.join(ITK_BUILD_DIR, 'lib'))

    # Append TubeTK libs
    sys.path.append(os.path.join(TubeTK_BUILD_DIR, 'TubeTK-build/lib/TubeTK'))
    import itk

from itk import TubeTKITK
import sys

def main():
  if len(sys.argv) != 6:
    print("Usage: %s InputImage OutputVesselMask OutputNotVesselMask gap notVesselWidth"%sys.argv[0])
    return 1
  inputImage=sys.argv[1]
  outputVesselMask=sys.argv[2]
  outputNotVesselMask=sys.argv[3]
  gap=float(sys.argv[4])
  notVesselWidth=float(sys.argv[5])

  reader=itk.ImageFileReader.New(FileName=inputImage)
  reader.Update()
  trainingMask=TubeTKITK.ComputeTrainingMask.New(reader)
  trainingMask.SetGap(gap)
  trainingMask.SetNotVesselWidth(notVesselWidth)
  trainingMask.Update()
  writer=itk.ImageFileWriter.New(trainingMask,FileName=outputVesselMask)
  writer.Update()

  writer=itk.ImageFileWriter.New(trainingMask.GetNotVesselMask(),FileName=outputNotVesselMask)
  writer.Update()


if __name__ == "__main__":
  sys.exit(main())
