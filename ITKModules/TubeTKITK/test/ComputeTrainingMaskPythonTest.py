#!/usr/bin/env python

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
