#!/usr/bin/python

###########################################################################
# ConvertData.py :
#
# Shrink input images into slabs. For each input, the associated
# expert segmentation is first converted into an image and then shrunk into
# slabs the same way.
#
###########################################################################

import sys
import os
import glob
import subprocess

# Relative path variables
caffe_root = "./"  # this file should be run from {caffe_root} (otherwise change this line)
hardDrive_root = "/media/lucas/krs0014/"

# Shrink images
def shrink( inputImage, expertImage, outputImagePrefix ):
  subprocess.call( ["ShrinkImage", "-n 512,512,10",
                    "-o 0,0,20",
                    "-p", outputImagePrefix + "_points.mha",
                    inputImage,
                    outputImagePrefix + ".mha"] )
  subprocess.call( ["ShrinkImage", "-n 512,512,10",
                    expertImage,
                    outputImagePrefix + "_expert.mha"] )

# Convert TRE file to Image
def convert( templateFile, tubeFile, outputFile ):
  subprocess.call( ["ConvertTubesToImage", "-r",
                    templateFile,
                    tubeFile,
                    outputFile ] )

# Copy infile to outFile
def copy( inFile, outFile ):
  subprocess.call( ["cp", inFile, outFile] )


def main(argv):
  # Input data directories
  controls = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/Controls/*"
  tumors = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/LargeTumor/*"
  # Output data directories
  expertDir = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/expert/"
  controlsOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/controls/"
  tumorsOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/tumors/"

  # Sanity checks
  if not os.path.exists( hardDrive_root + "SegmentVesselsUsingNeuralNetworks/" ):
    subprocess.call( ["mkdir" , hardDrive_root + "SegmentVesselsUsingNeuralNetworks/"] )
  if not os.path.exists( expertDir ):
    subprocess.call( ["mkdir" , expertDir] )
  if not os.path.exists( controlsOutputDir ):
    subprocess.call( ["mkdir" , controlsOutputDir] )
  if not os.path.exists( tumorsOutputDir ):
    subprocess.call( ["mkdir" , tumorsOutputDir] )

  # Files to process
  controlFiles = glob.glob( os.path.join( controls, "*.mhd" ) )
  tumorFiles = glob.glob( os.path.join( tumors, "*.mhd" ) )

  # Output directories where to copy the data for training and testing
  i=0
  trainOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/training/"
  testOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/"

  # Process control files
  for file in controlFiles:
      # Create filenames
      filePrefix = os.path.basename(os.path.splitext(file)[0])
      fileDirectory = os.path.dirname(os.path.abspath(file))
      templateFile = file
      tubeFile = fileDirectory + "/TRE/" + filePrefix + ".tre"
      expertFile = expertDir + filePrefix + "_expert.mha"

      # Process
      convert( templateFile, tubeFile, expertFile )
      shrink( templateFile, expertFile, controlsOutputDir + filePrefix )

      # Mixing controls for testing and training.
      if i%2 == 0 :
        copy( controlsOutputDir + filePrefix + ".mha", trainOutputDir+"images/" + filePrefix + ".mha")
        copy( controlsOutputDir + filePrefix + "_points.mha", trainOutputDir+"points/" + filePrefix + "_points.mha")
        copy( controlsOutputDir + filePrefix + "_expert.mha", trainOutputDir+"expert/" + filePrefix + "_expert.mha")
      else :
        copy( controlsOutputDir + filePrefix + ".mha", testOutputDir+"images/" + filePrefix + ".mha")
        copy( controlsOutputDir + filePrefix + "_points.mha", testOutputDir+"points/" + filePrefix + "_points.mha")
        copy( controlsOutputDir + filePrefix + "_expert.mha", testOutputDir+"expert/" + filePrefix + "_expert.mha")
      i = i+1

  i=0
  # Process tumor files
  for file in tumorFiles:
      # Create filenames
      filePrefix = os.path.basename(os.path.splitext(file)[0])
      fileDirectory = os.path.dirname(os.path.abspath(file))
      templateFile = file
      tubeFile = fileDirectory + "/TRE/" + filePrefix + ".tre"
      expertFile = expertDir + filePrefix + "_expert.mha"

      # Process
      convert( templateFile, tubeFile, expertFile )
      shrink( templateFile, expertFile, tumorsOutputDir + filePrefix )

      # Mixing tumors for testing and training.
      if i%2 == 0 :
        copy( tumorsOutputDir + filePrefix + ".mha", trainOutputDir+"images/" + filePrefix + ".mha")
        copy( tumorsOutputDir + filePrefix + "_points.mha", trainOutputDir+"points/" + filePrefix + "_points.mha")
        copy( tumorsOutputDir + filePrefix + "_expert.mha", trainOutputDir+"expert/" + filePrefix + "_expert.mha")
      else :
        copy( tumorsOutputDir + filePrefix + ".mha", testOutputDir+"images/" + filePrefix + ".mha")
        copy( tumorsOutputDir + filePrefix + "_points.mha", testOutputDir+"points/" + filePrefix + "_points.mha")
        copy( tumorsOutputDir + filePrefix + "_expert.mha", testOutputDir+"expert/" + filePrefix + "_expert.mha")
      i = i+1

if __name__ == "__main__":
  main( sys.argv[0:] )
