#!/usr/bin/python

###########################################################################
# saveSlabs.py :
#
# Save each slabs as a single image and compute expert training mask.
#
###########################################################################

import os
import glob
import sys
import subprocess
import json

import matplotlib.pyplot as plt

# Append ITK libs
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Modules/ThirdParty/VNL/src/vxl/lib'))
import itk


# Compute Training mask
def computeTrainingMask(expertInput, trainingMaskOutput):
    subprocess.call(["SegmentBinaryImageSkeleton",
                     expertInput,
                     expertInput + "tmp.png"])
    subprocess.call(["ImageMath", expertInput,
                     "-t", "0", "0.5", "0", "255",
                     "-W", "0", expertInput + "tmp2.png",
                     "-M", "1", "1", "255", "0",
                     "-a", "1", "-1", expertInput + "tmp2.png",
                     "-a", "0.5", "255", expertInput + "tmp.png",
                     "-W", "0", trainingMaskOutput])
    subprocess.call(["rm", "-rf", expertInput + "tmp.png"])
    subprocess.call(["rm", "-rf", expertInput + "tmp2.png"])

    # WARNING: Couldn't write to PNG using the implemented CLI
    # subprocess.call( ["ComputeTrainingMask",
    # expertInput,
    # trainingMaskOutput,
    #"--notVesselWidth","1"] )


# Save slab
def saveSlabs(mhaFileList):

    PixelType = itk.F
    Dimension = 3
    ImageType = itk.Image[PixelType, Dimension]
    ReaderType = itk.ImageFileReader[ImageType]

    for mhaFile in mhaFileList:

        print(mhaFile)

        reader = ReaderType.New()
        reader.SetFileName(str(mhaFile))
        reader.Update()
        img = reader.GetOutput()
        buf = itk.PyBuffer[ImageType].GetArrayFromImage(img)

        # filenames definition
        filePrefix = os.path.basename(os.path.splitext(mhaFile)[0])
        fileDirectory = os.path.dirname(os.path.abspath(mhaFile))

        # Iterate through each slab
        for i in range(0, buf.shape[0]):

            outputImage = os.path.join(
                fileDirectory, str(i) + "_" + filePrefix + ".png")
            plt.imsave(outputImage, buf[i, :, :], cmap='Greys_r')

            # Compute the expert training mask
            if "expert" in outputImage:
                computeTrainingMask(outputImage, outputImage)

# Path variables
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']
proj_rel_path = script_params['PROJECT_REL_PATH']

caffe_proj_root = os.path.join(caffe_root, "data", proj_rel_path)
hardDrive_proj_root = os.path.join(hardDrive_root, proj_rel_path)

trainDataDir = os.path.join(hardDrive_proj_root, "training")
valDataDir = os.path.join(hardDrive_proj_root, "testing")


saveSlabs(glob.glob(os.path.join(trainDataDir, "images", "*.mha")))
saveSlabs(glob.glob(os.path.join(trainDataDir, "expert", "*.mha")))

saveSlabs(glob.glob(os.path.join(valDataDir, "images", "*.mha")))
saveSlabs(glob.glob(os.path.join(valDataDir, "expert", "*.mha")))
