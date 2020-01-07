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
proj_rel_path = script_params['PROJECT_REL_PATH']

caffe_proj_root = os.path.join(caffe_root, "data", proj_rel_path)
hardDrive_proj_root = os.path.join(hardDrive_root, proj_rel_path)


# Shrink images
def shrink(inputImage, expertImage, outputImagePrefix):

    subprocess.call(["ShrinkImage", "-n 512,512,10",
                     "-o 0,0,20",
                     "-p", outputImagePrefix + "_zslab_points.mha",
                     inputImage,
                     outputImagePrefix + "_zslab.mha"])
    subprocess.call(["ShrinkImage", "-n 512,512,10",
                     expertImage,
                     outputImagePrefix + "_zslab_expert.mha"])

    subprocess.call(["ImageMath",
                     "-w", outputImagePrefix + ".mha",
                     inputImage])


# Convert TRE file to Image
def convertTubesToImage(templateFile, tubeFile, outputFile):

    subprocess.call(["ConvertTubesToImage", "-r",
                     templateFile,
                     tubeFile,
                     outputFile])


# Copy infile to outFile and create dirs if not present
def copy(inFile, outFile):

    # create path if it doesnt exist
    out_path = os.path.dirname(outFile)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # copy file
    shutil.copyfile(inFile, outFile)


def main(argv):

    # Input data directories
    controls = os.path.join(caffe_proj_root, "Controls/*")
    tumors = os.path.join(caffe_proj_root, "LargeTumor/*")

    # Output data directories
    expertDir = os.path.join(caffe_proj_root, "expert")
    controlsOutputDir = os.path.join(hardDrive_proj_root, "controls")
    tumorsOutputDir = os.path.join(hardDrive_proj_root, "tumors")

    # Sanity checks
    if not os.path.exists(hardDrive_proj_root):
        os.mkdir(hardDrive_proj_root)

    if not os.path.exists(expertDir):
        os.mkdir(expertDir)

    if not os.path.exists(controlsOutputDir):
        os.mkdir(controlsOutputDir)

    if not os.path.exists(tumorsOutputDir):
        os.mkdir(tumorsOutputDir)

    # Files to process
    controlMhdFiles = glob.glob(os.path.join(controls, "*.mhd"))
    tumorMhdFiles = glob.glob(os.path.join(tumors, "*.mhd"))

    # Output directories where to copy the data for training and testing
    trainOutputDir = os.path.join(hardDrive_proj_root, "training")
    testOutputDir = os.path.join(hardDrive_proj_root, "testing")

    # Process control files
    i = 0

    for mhdFile in controlMhdFiles:

        # Create filenames
        filePrefix = os.path.basename(os.path.splitext(mhdFile)[0])
        fileDir = os.path.dirname(os.path.abspath(mhdFile))

        print("control file %d/%d : " %
              (i + 1, len(controlMhdFiles)), filePrefix)

        treFile = os.path.join(fileDir, "TRE", filePrefix + ".tre")
        expertFile = os.path.join(
            controlsOutputDir, filePrefix + "_expert.mha")

        # Process
        convertTubesToImage(mhdFile, treFile, expertFile)
        shrink(mhdFile, expertFile,
               os.path.join(controlsOutputDir, filePrefix))

        # Mixing controls for testing and training.
        if i % 2 == 0:
            curOutputDir = trainOutputDir
        else:
            curOutputDir = testOutputDir

        copy(os.path.join(controlsOutputDir, filePrefix + ".mha"),
             os.path.join(curOutputDir, "images", filePrefix + ".mha"))

        copy(os.path.join(controlsOutputDir, filePrefix + "_zslab.mha"),
             os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"))

        copy(os.path.join(controlsOutputDir, filePrefix + "_zslab_points.mha"),
             os.path.join(curOutputDir, "points", filePrefix + "_zslab_points.mha"))

        copy(os.path.join(controlsOutputDir, filePrefix + "_expert.mha"),
             os.path.join(curOutputDir, "expert", filePrefix + "_expert.mha"))

        copy(os.path.join(controlsOutputDir, filePrefix + "_zslab_expert.mha"),
             os.path.join(curOutputDir, "expert", filePrefix + "_zslab_expert.mha"))

        i += 1

    # Process tumor files
    i = 0

    for mhdFile in tumorMhdFiles:

        # Create filenames
        filePrefix = os.path.basename(os.path.splitext(mhdFile)[0])
        fileDir = os.path.dirname(os.path.abspath(mhdFile))

        print("tumor file %d / %d : " %
              (i + 1, len(tumorMhdFiles)), filePrefix)

        treFile = os.path.join(fileDir, "TRE", filePrefix + ".tre")
        expertFile = os.path.join(
            tumorsOutputDir, filePrefix + "_expert.mha")

        # Process
        convertTubesToImage(mhdFile, treFile, expertFile)
        shrink(mhdFile, expertFile,
               os.path.join(tumorsOutputDir, filePrefix))

        # Mixing controls for testing and training.
        if i % 2 == 0:
            curOutputDir = trainOutputDir
        else:
            curOutputDir = testOutputDir

        copy(os.path.join(tumorsOutputDir, filePrefix + ".mha"),
             os.path.join(curOutputDir, "images", filePrefix + ".mha"))

        copy(os.path.join(tumorsOutputDir, filePrefix + "_zslab.mha"),
             os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"))

        copy(os.path.join(tumorsOutputDir, filePrefix + "_zslab_points.mha"),
             os.path.join(curOutputDir, "points", filePrefix + "_zslab_points.mha"))

        copy(os.path.join(tumorsOutputDir, filePrefix + "_expert.mha"),
             os.path.join(curOutputDir, "expert", filePrefix + "_expert.mha"))

        copy(os.path.join(tumorsOutputDir, filePrefix + "_zslab_expert.mha"),
             os.path.join(curOutputDir, "expert", filePrefix + "_zslab_expert.mha"))

        i += 1

if __name__ == "__main__":
    main(sys.argv[0:])
