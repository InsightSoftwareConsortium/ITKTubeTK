#!/usr/bin/python

###########################################################################
# PrepareTrainingData.py :
#
# Prepares training data for CNN
#
###########################################################################

import glob
import json
import math
import os
import random
import subprocess
import shutil
import sys

import skimage.io

import utils

# Append ITK libs
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Modules/ThirdParty/VNL/src/vxl/lib'))
import itk

# Define paths
script_params = json.load(open('params.json'))

caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']
proj_rel_path = script_params['PROJECT_REL_PATH']

caffe_proj_root = os.path.join(caffe_root, "data", proj_rel_path)
hardDrive_proj_root = os.path.join(hardDrive_root, proj_rel_path)


# Create segmentation mask from tre file
def createExpertSegmentationMask(inputImageFile, treFile, outputExpertSegFile):

    # convert tre file to image
    subprocess.call(["ConvertTubesToImage", "-r",
                     inputImageFile,
                     treFile,
                     outputExpertSegFile])


# Shrink images
def shrink(inputImage, expertImage, outputImagePrefix):
    """Shrink inputImage and expertImage according to script_params.
    Output (where '*' stands for outputImagePrefix):
    - *.mha: MHA copy of inputImage
    - *_zslab_points.mha: Image of max points (vectors)
    - *_zslab.mha: Shrunk inputImage
    - *_zslab_expert.mha: Shrunk expertImage

    """

    shrinked_size = "512,512,%d" % script_params["NUM_SLABS"]
    window_overlap = "0,0,%d" % script_params["SLAB_OVERLAP"]

    subprocess.call(["ShrinkImage",
                     "-n", shrinked_size,
                     "-o", window_overlap,
                     "-p", outputImagePrefix + "_zslab_points.mha",
                     inputImage,
                     outputImagePrefix + "_zslab.mha"])

    subprocess.call(["ShrinkImage",
                     "-n", shrinked_size,
                     "-o", window_overlap,
                     "-i", outputImagePrefix + "_zslab_points.mha",
                     expertImage,
                     outputImagePrefix + "_zslab_expert.mha"])

    """
    subprocess.call(["ShrinkImage",
                     "-n", shrinked_size,
                     expertImage,
                     outputImagePrefix + "_zslab_expert.mha"])
    """

    # save a copy of inputImage in .mha format
    subprocess.call(["ImageMath",
                     "-w", outputImagePrefix + ".mha",
                     inputImage])


def createZMIPSlabsFor(name, inputDir, outputDir):
    # Sanity check
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Process files
    printSectionHeader('Creating Z-MIP slabs for %ss' % name)

    mhdFiles = glob.glob(os.path.join(inputDir, "*", "*.mhd"))

    i = 0

    for mhdFile in mhdFiles:

        print("\n%s file %d/%d : %s" %
              (name, i + 1, len(mhdFiles), mhdFile))

        fileName = os.path.basename(os.path.splitext(mhdFile)[0])
        fileDir = os.path.dirname(os.path.abspath(mhdFile))

        treFile = os.path.join(fileDir, "TRE", fileName + ".tre")
        expertSegFile = os.path.join(outputDir,
                                     fileName + "_expert.mha")

        # Process
        createExpertSegmentationMask(mhdFile, treFile, expertSegFile)

        shrink(mhdFile, expertSegFile,
               os.path.join(outputDir, fileName))

        i += 1


# create z-mip slabs
def createZMIPSlabs():

    # Input data directories where mha/mhd and associated tre files are located
    controlInputDir = os.path.join(caffe_proj_root, "Controls")
    tumorInputDir = os.path.join(caffe_proj_root, "LargeTumor")

    # Output data directories
    controlOutputDir = os.path.join(hardDrive_proj_root, "controls")
    tumorOutputDir = os.path.join(hardDrive_proj_root, "tumors")

    # Process control files
    createZMIPSlabsFor('control', controlInputDir, controlOutputDir)

    # Process tumor files
    createZMIPSlabsFor('tumor', tumorInputDir, tumorOutputDir)

# Compute Training mask
def computeTrainingMask(expertSegMask, outputTrainingMask):

    subprocess.call(["SegmentBinaryImageSkeleton",
                     expertSegMask,
                     expertSegMask + "_skel.png"])

    subprocess.call([
        "ImageMath", expertSegMask,
         # dilate vessel mask with kernel radius = 1
         "-M", "1", "1", "255", "0",
        # subtract vessel mask from dilated version to get vessel boundary
        "-a", "1", "-1", expertSegMask,
        # create training mask with vessel center-line (=255) and boundary (=128)
        "-a", "0.5", "255", expertSegMask + "_skel.png",
        # write training mask
        "-W", "0", outputTrainingMask])

    subprocess.call(["rm", "-rf", expertSegMask + "_skel.png"])

    # WARNING: Couldn't write to PNG using the implemented CLI
    # subprocess.call( ["ComputeTrainingMask",
    # expertSegMask,
    # outputTrainingMask,
    # "--notVesselWidth","1"] )


# save each of the z-mip slabs from .mha files as a .png file
def saveSlabs(mhaFileList):

    PixelType = itk.F
    Dimension = 3
    ImageType = itk.Image[PixelType, Dimension]
    ReaderType = itk.ImageFileReader[ImageType]

    for mhaFile in mhaFileList:

        print 'saving slabs of %s' % mhaFile

        reader = ReaderType.New()
        reader.SetFileName(str(mhaFile))
        reader.Update()
        img = reader.GetOutput()
        buf = itk.PyBuffer[ImageType].GetArrayFromImage(img)

        # convert to [0, 255] range
        buf = 255.0 * (buf - buf.min()) / (buf.max() - buf.min())
        buf = buf.astype('uint8')

        # file names definition
        fileName = os.path.basename(os.path.splitext(mhaFile)[0])
        fileDir = os.path.dirname(os.path.abspath(mhaFile))

        # Iterate through each slab
        for i in range(buf.shape[0]):

            outputImage = os.path.join(
                fileDir, str(i) + "_" + fileName + ".png")

            skimage.io.imsave(outputImage, buf[i, :, :])


# assign control and tumor volumes equally to training and testing
def splitControlTumorData():

    # Input data directories
    controlInputDir = os.path.join(caffe_proj_root, "Controls")
    tumorInputDir = os.path.join(caffe_proj_root, "LargeTumor")

    # Output data directories
    controlOutputDir = os.path.join(hardDrive_proj_root, "controls")
    tumorOutputDir = os.path.join(hardDrive_proj_root, "tumors")

    trainOutputDir = os.path.join(hardDrive_proj_root, "training")
    testOutputDir = os.path.join(hardDrive_proj_root, "testing")

    # Sanity checks
    if not os.path.exists(trainOutputDir):
        os.makedirs(trainOutputDir)

    if not os.path.exists(testOutputDir):
        os.makedirs(testOutputDir)

    # Process control files
    printSectionHeader('Splitting control data into training and testing')

    controlMhdFiles = glob.glob(os.path.join(controlInputDir, "*", "*.mhd"))

    i = 0

    for mhdFile in controlMhdFiles:

        print("\ncontrol file %d/%d : %s" %
              (i + 1, len(controlMhdFiles), mhdFile))

        filePrefix = os.path.basename(os.path.splitext(mhdFile)[0])

        # Split equally for training and testing
        if i % 2 == 0:
            curOutputDir = trainOutputDir
        else:
            curOutputDir = testOutputDir

        # copy input volume
        utils.copy(os.path.join(controlOutputDir, filePrefix + ".mha"),
                   os.path.join(curOutputDir, "images", filePrefix + ".mha"))

        # copy z-mip slab volume
        utils.copy(os.path.join(controlOutputDir, filePrefix + "_zslab.mha"),
                   os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"))

        # copy z-mip slab point map
        utils.copy(os.path.join(controlOutputDir, filePrefix + "_zslab_points.mha"),
                   os.path.join(curOutputDir, "points", filePrefix + "_zslab_points.mha"))

        # copy expert volume
        utils.copy(os.path.join(controlOutputDir, filePrefix + "_expert.mha"),
                   os.path.join(curOutputDir, "expert", filePrefix + "_expert.mha"))

        # copy expert z-mip slab volume
        utils.copy(os.path.join(controlOutputDir, filePrefix + "_zslab_expert.mha"),
                   os.path.join(curOutputDir, "expert", filePrefix + "_zslab_expert.mha"))

        # save slabs as pngs
        saveSlabs([
            os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"),
            os.path.join(curOutputDir, "expert",
                         filePrefix + "_zslab_expert.mha")
        ])

        i += 1

    # Process tumor files
    printSectionHeader('Splitting tumor data into training and testing')

    tumorMhdFiles = glob.glob(os.path.join(tumorInputDir, "*", "*.mhd"))

    i = 0

    for mhdFile in tumorMhdFiles:

        print("\ntumor file %d / %d : %s" %
              (i + 1, len(tumorMhdFiles), mhdFile))

        filePrefix = os.path.basename(os.path.splitext(mhdFile)[0])

        # Split equally for training and testing
        if i % 2 == 0:
            curOutputDir = trainOutputDir
        else:
            curOutputDir = testOutputDir

        # copy input volume
        utils.copy(os.path.join(tumorOutputDir, filePrefix + ".mha"),
                   os.path.join(curOutputDir, "images", filePrefix + ".mha"))

        # copy z-mip slab volume
        utils.copy(os.path.join(tumorOutputDir, filePrefix + "_zslab.mha"),
                   os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"))

        # copy z-mip slab point map
        utils.copy(os.path.join(tumorOutputDir, filePrefix + "_zslab_points.mha"),
                   os.path.join(curOutputDir, "points", filePrefix + "_zslab_points.mha"))

        # copy expert volume
        utils.copy(os.path.join(tumorOutputDir, filePrefix + "_expert.mha"),
                   os.path.join(curOutputDir, "expert", filePrefix + "_expert.mha"))

        # copy expert z-mip slab volume
        utils.copy(os.path.join(tumorOutputDir, filePrefix + "_zslab_expert.mha"),
                   os.path.join(curOutputDir, "expert", filePrefix + "_zslab_expert.mha"))

        # save slabs as pngs
        saveSlabs([
            os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"),
            os.path.join(curOutputDir, "expert",
                         filePrefix + "_zslab_expert.mha")
        ])

        i += 1


# Extracts +ve (vessel center) and -ve (background) patches from image
def extractPatchesFromImage(rootDir, imageName, outputDir, patchListFile):

    # patch/window radius
    w = script_params['PATCH_RADIUS']

    # fraction of negatives near vessel boundary
    vessel_bnd_neg_frac = script_params['NEGATIVES_NEAR_VESSEL_BOUNDARY']

    # read input image
    inputImageFile = os.path.join(rootDir, "images", imageName + '.png')
    inputImage = skimage.io.imread(inputImageFile)

    # read expert segmented label mask
    expertSegFile = os.path.join(rootDir, "expert", imageName + '_expert.png')
    trainingMaskFile = os.path.join(rootDir, "expert", imageName + '_train.png')

    computeTrainingMask(expertSegFile, trainingMaskFile)

    trainingMask = skimage.io.imread(trainingMaskFile)

    # Iterate through expert mask and find pos/neg patches
    numVesselPatches = 0

    vesselBndInd = []  # Indices of background pixels near vessel boundary
    bgndInd = []  # Indices of all other background pixels
    subsample = 1  # Increase to reduce time for debugging

    for i in range(w, trainingMask.shape[0] - w - 1, subsample):
        for j in range(w, trainingMask.shape[1] - w - 1, subsample):

            # Vessel center-line pixel (positive)
            if trainingMask[i, j] > 0.6 * 255:

                filename = os.path.join(
                    "1", imageName + "_" + str(i) + "_" + str(j) + ".png")

                patchListFile.write(filename + " " + str(1) + "\n")

                skimage.io.imsave(os.path.join(outputDir, filename),
                                  inputImage[i - w : i + w + 1, j - w : j + w + 1])

                numVesselPatches += 1

            # Vessel bound pixel (negative)
            elif trainingMask[i, j] > 0:

                vesselBndInd.append([i, j])

            # Background pixel (negative)
            else:

                bgndInd.append([i, j])

    # Pick a subset of background (negative) patches near vessel boundary
    numVesselBndPatches = int(
        math.ceil(vessel_bnd_neg_frac * numVesselPatches))

    selVesselBndInd = random.sample(vesselBndInd, numVesselBndPatches)

    for [i, j] in selVesselBndInd:

        filename = os.path.join("0", imageName +
                                "_" + str(i) + "_" + str(j) + ".png")

        patchListFile.write(filename + " " + str(0) + "\n")

        skimage.io.imsave(os.path.join(outputDir, filename),
                          inputImage[i - w : i + w + 1, j - w : j + w + 1])

    # Pick rest of background (negative) patches away from vessel boundary
    numOtherBgndPatches = \
        int(math.ceil((1.0 - vessel_bnd_neg_frac) * numVesselPatches))

    selBgndInd = random.sample(bgndInd, numOtherBgndPatches)

    for [i, j] in selBgndInd:

        filename = os.path.join("0", imageName +
                                "_" + str(i) + "_" + str(j) + ".png")

        patchListFile.write(filename + " " + str(0) + "\n")

        skimage.io.imsave(os.path.join(outputDir, filename),
                          inputImage[i - w : i + w + 1, j - w : j + w + 1])

    print 'Number of positive patches = ', numVesselPatches


# convert train/test images to patches
def createTrainTestPatches():

    # create training patches
    printSectionHeader('Creating training patches')

    trainDataDir = os.path.join(hardDrive_proj_root, "training")
    trainImageFiles = glob.glob(os.path.join(trainDataDir, "images", "*.png"))

    trainPatchesDir = os.path.join(trainDataDir, "patches")
    for i in range(2):
        if not os.path.exists(os.path.join(trainPatchesDir, str(i))):
            os.makedirs(os.path.join(trainPatchesDir, str(i)))

    trainPatchListFile = open(os.path.join(trainPatchesDir, "train.txt"), "w")
    trainPatchListFile.truncate()

    i = 0

    for imageFile in trainImageFiles:

        print('\nCreating patches for train file %d/%d - %s' %
              (i + 1, len(trainImageFiles), imageFile))

        imageName = os.path.basename(os.path.splitext(imageFile)[0])

        extractPatchesFromImage(trainDataDir, imageName, trainPatchesDir,
                                trainPatchListFile)

        i += 1

    trainPatchListFile.close()

    # create testing patches
    printSectionHeader('Creating testing patches')

    testDataDir = os.path.join(hardDrive_proj_root, "testing")
    testImageFiles = glob.glob(os.path.join(testDataDir, "images", "*.png"))

    testPatchesDir = os.path.join(testDataDir, "patches")
    for i in range(2):
        if not os.path.exists(os.path.join(testPatchesDir, str(i))):
            os.makedirs(os.path.join(testPatchesDir, str(i)))

    testPatchListFile = open(os.path.join(testPatchesDir, "val.txt"), "w+")
    testPatchListFile.truncate()  # Erase file

    i = 0
    for imageFile in testImageFiles:

        print('\nCreating patches for test file %d/%d - %s' %
              (i + 1, len(testImageFiles), imageFile))

        imageName = os.path.basename(os.path.splitext(imageFile)[0])

        extractPatchesFromImage(testDataDir, imageName, testPatchesDir,
                                testPatchListFile)

        i += 1

    testPatchListFile.close()


# create lmdb
def createTrainTestLmdb():

    printSectionHeader('Creating LMDBs for train and test data')

    caffe_tools_dir = os.path.join(caffe_root, "build", "tools")
    convert_imageset_exec = os.path.join(caffe_tools_dir, "convert_imageset")

    # create training lmdb
    print('Creating training lmdb ...\n')

    train_patches_dir = os.path.join(
        hardDrive_proj_root, "training", "patches/")

    train_lmdb_dir = os.path.join(caffe_proj_root, "Net_TrainData")
    if os.path.exists(train_lmdb_dir):
        shutil.rmtree(train_lmdb_dir)

    subprocess.call([convert_imageset_exec,
                     "--shuffle",
                     "--gray",
                     train_patches_dir,
                     os.path.join(train_patches_dir, "train.txt"),
                     train_lmdb_dir])

    # create testing lmdb
    print('Creating testing lmdb ...\n')

    test_patches_dir = os.path.join(hardDrive_proj_root, "testing", "patches/")

    test_lmdb_dir = os.path.join(caffe_proj_root, "Net_ValData")
    if os.path.exists(test_lmdb_dir):
        shutil.rmtree(test_lmdb_dir)

    subprocess.call([convert_imageset_exec,
                     "--shuffle",
                     "--gray",
                     test_patches_dir,
                     os.path.join(test_patches_dir, "val.txt"),
                     test_lmdb_dir])


def printSectionHeader(title):

    print("\n" + "#" * (len(title) + 4))
    print('# ' + title)
    print("#" * (len(title) + 4) + "\n")


def run():

    # create z-mip slabs
    createZMIPSlabs()

    # assign control and tumor volumes equally to training and testing
    # Note: this must be called after createZMIPSlabs()
    splitControlTumorData()

    # convert train/test images to patches
    createTrainTestPatches()

    # create lmdb database
    createTrainTestLmdb()


if __name__ == "__main__":
    run()
