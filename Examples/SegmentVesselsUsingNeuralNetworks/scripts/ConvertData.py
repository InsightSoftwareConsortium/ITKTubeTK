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
import shutil
import json

# Define paths
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']

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

  # create path if it doesnt exist
  out_path = os.path.dirname( outFile )
  if not os.path.exists( out_path ):
    os.makedirs( out_path )

  # copy file
  shutil.copyfile(inFile, outFile)


def main(argv):

  # Input data directories
  controls = os.path.join(caffe_root, "data/SegmentVesselsUsingNeuralNetworks/Controls/*")
  tumors = os.path.join(caffe_root, "data/SegmentVesselsUsingNeuralNetworks/LargeTumor/*")

  # Output data directories
  expertDir = os.path.join(caffe_root, "data/SegmentVesselsUsingNeuralNetworks/expert/")
  controlsOutputDir = os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/controls/")
  tumorsOutputDir = os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/tumors/")

  # Sanity checks
  if not os.path.exists( os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/") ):
    os.mkdir( os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/") )

  if not os.path.exists( expertDir ):
    os.mkdir( expertDir )

  if not os.path.exists( controlsOutputDir ):
    os.mkdir( controlsOutputDir )

  if not os.path.exists( tumorsOutputDir ):
    os.mkdir( tumorsOutputDir )

  # Files to process
  controlFiles = glob.glob( os.path.join( controls, "*.mhd" ) )
  tumorFiles = glob.glob( os.path.join( tumors, "*.mhd" ) )

  # Output directories where to copy the data for training and testing
  i=0

  trainOutputDir = os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/training/")
  testOutputDir = os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/testing/")

  # Process control files
  for file in controlFiles:
      # Create filenames
      filePrefix = os.path.basename(os.path.splitext(file)[0])
      fileDirectory = os.path.dirname(os.path.abspath(file))
      templateFile = file
      tubeFile = os.path.join(fileDirectory, "TRE", filePrefix + ".tre")
      expertFile = os.path.join(expertDir, filePrefix + "_expert.mha")

      # Process
      convert( templateFile, tubeFile, expertFile )
      shrink( templateFile, expertFile, os.path.join(controlsOutputDir, filePrefix) )

      # Mixing controls for testing and training.
      if i%2 == 0 :
        copy( os.path.join(controlsOutputDir, filePrefix + ".mha"), os.path.join(trainOutputDir, "images", filePrefix + ".mha") )
        copy( os.path.join(controlsOutputDir, filePrefix + "_points.mha"), os.path.join(trainOutputDir, "points", filePrefix + "_points.mha") )
        copy( os.path.join(controlsOutputDir, filePrefix + "_expert.mha"), os.path.join(trainOutputDir, "expert", filePrefix + "_expert.mha") )
      else :
        copy( os.path.join(controlsOutputDir, filePrefix + ".mha"), os.path.join(testOutputDir, "images", filePrefix + ".mha") )
        copy( os.path.join(controlsOutputDir, filePrefix + "_points.mha"), os.path.join(testOutputDir, "points", filePrefix + "_points.mha") )
        copy( os.path.join(controlsOutputDir, filePrefix + "_expert.mha"), os.path.join(testOutputDir, "expert", filePrefix + "_expert.mha") )
      i = i+1

  i=0
  # Process tumor files
  for file in tumorFiles:
      # Create filenames
      filePrefix = os.path.basename(os.path.splitext(file)[0])
      fileDirectory = os.path.dirname(os.path.abspath(file))
      templateFile = file
      tubeFile = os.path.join(fileDirectory, "TRE", filePrefix + ".tre")
      expertFile = os.path.join(expertDir, filePrefix + "_expert.mha")

      # Process
      convert( templateFile, tubeFile, expertFile )
      shrink( templateFile, expertFile, tumorsOutputDir + filePrefix )

      # Mixing tumors for testing and training.
      if i%2 == 0 :
        copy( os.path.join(tumorsOutputDir, filePrefix + ".mha"), os.path.join(trainOutputDir, "images", filePrefix + ".mha") )
        copy( os.path.join(tumorsOutputDir, filePrefix + "_points.mha"), os.path.join(trainOutputDir, "points", filePrefix + "_points.mha") )
        copy( os.path.join(tumorsOutputDir, filePrefix + "_expert.mha"), os.path.join(trainOutputDir, "expert", filePrefix + "_expert.mha") )
      else :
        copy( os.path.join(tumorsOutputDir, filePrefix + ".mha"), os.path.join(testOutputDir, "images", filePrefix + ".mha") )
        copy( os.path.join(tumorsOutputDir, filePrefix + "_points.mha"), os.path.join(testOutputDir, "points", filePrefix + "_points.mha") )
        copy( os.path.join(tumorsOutputDir, filePrefix + "_expert.mha"), os.path.join(testOutputDir, "expert",  filePrefix + "_expert.mha") )
      i = i+1

if __name__ == "__main__":
  main( sys.argv[0:] )
