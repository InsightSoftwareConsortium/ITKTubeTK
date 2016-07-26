#!/usr/bin/python

###########################################################################
# saveSlabs.py :
#
# Save each slabs as a single image and compute expert training mask.
#
###########################################################################

import os, glob, sys, subprocess

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

import numpy as np

## Append ITK libs
sys.path.append('/home/lucas/Projects/ITK-Release/Wrapping/Generators/Python')
sys.path.append('/home/lucas/Projects/ITK-Release/Modules/ThirdParty/VNL/src/vxl/lib')
import itk

# Compute Training mask
def computeTrainingMask( expertInput, trainingMaskOutput ):
  subprocess.call( ["SegmentBinaryImageSkeleton",
                    expertInput,
                    expertInput+"tmp.png"] )
  subprocess.call( ["ImageMath", expertInput,
                    "-t", "0", "0.5", "0", "255",
                    "-W", "0", expertInput+"tmp2.png",
                    "-M", "1", "1", "255", "0",
                    "-a", "1", "-1", expertInput+"tmp2.png",
                    "-a", "0.5", "255", expertInput+"tmp.png",
                    "-W", "0", trainingMaskOutput] )
  subprocess.call( ["rm", "-rf",expertInput+"tmp.png"])
  subprocess.call( ["rm", "-rf",expertInput+"tmp2.png"])

  # WARNING: Couldn't write to PNG using the implemented CLI
  #subprocess.call( ["ComputeTrainingMask",
                    #expertInput,
                    #trainingMaskOutput,
                    #"--notVesselWidth","1"] )

# Save slab
def saveSlabs( fileList ):
  for file in fileList:
    print file
    reader.SetFileName(file)
    reader.Update()
    img = reader.GetOutput()
    buf = itk.PyBuffer[ImageType].GetArrayFromImage(img)
    # filenames definition
    filePrefix = os.path.basename(os.path.splitext(file)[0])
    fileDirectory = os.path.dirname(os.path.abspath(file))
    fileDirectory = fileDirectory + "/"
    # Iterate through each slab
    for i in range(0,buf.shape[0]):
      outputImage = fileDirectory + str(i) + "_" + filePrefix + ".png"
      plt.imsave( outputImage, buf[i,:,:],cmap='Greys_r' )
      # Compute the expert training mask
      if "expert" in outputImage:
        computeTrainingMask( outputImage, outputImage )

# Path variables
hardDrive_root = "/media/lucas/krs0014/"

trainOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/training/*"
testOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/*"

trainFiles = glob.glob( os.path.join( trainOutputDir, "*.mha" ) )
testFiles = glob.glob( os.path.join( testOutputDir, "*.mha" ) )

PixelType = itk.F
Dimension = 3
ImageType=itk.Image[PixelType,Dimension]
ReaderType = itk.ImageFileReader[ImageType]

reader = ReaderType.New()

saveSlabs(trainFiles)
saveSlabs(testFiles)
