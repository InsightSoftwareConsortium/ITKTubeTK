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

# Where the input data is to be found, to be conceptually
# distinguished from its location in the caffe root directory
input_image_root = caffe_proj_root

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


def createZMIPSlabsForFile(mhdFile, outputDir):
    """Create slabs and related files corresponding to mhdFile in outputDir.
    Input:
    - $input/*.mhd: The image file header
    - $input/TRE/*.tre: The expert TRE file
    Output:
    - $output/*_expert.mha: The expert MHA file
    - : All output from shrink with $output/* as outputImagePrefix

    """
    fileName = os.path.basename(os.path.splitext(mhdFile)[0])
    fileDir = os.path.dirname(os.path.abspath(mhdFile))

    treFile = os.path.join(fileDir, "TRE", fileName + ".tre")
    expertSegFile = os.path.join(outputDir,
                                 fileName + "_expert.mha")

    # Process
    createExpertSegmentationMask(mhdFile, treFile, expertSegFile)

    shrink(mhdFile, expertSegFile,
           os.path.join(outputDir, fileName))

def createZMIPSlabs(name, inputDir, outputDir):
    """Process all image files in immediate subdirectories of inputDir to
    correspondingly prefixed images in outputDir.  outputDir is
    created if it doesn't already exist.  The subdirectory structure
    is not replicated.

    See the documentation of createZMIPSlabsForFile for the exact
    files created.

    """
    # Sanity check
    utils.ensureDirectoryExists(outputDir)

    # Process files
    printSectionHeader('Creating Z-MIP slabs for %ss' % name)

    mhdFiles = glob.glob(os.path.join(inputDir, "*", "*.mhd"))

    for i, mhdFile in enumerate(mhdFiles):

        print("\n%s file %d/%d : %s" %
              (name, i + 1, len(mhdFiles), mhdFile))

        createZMIPSlabsForFile(mhdFile, outputDir)


# create z-mip slabs
def createControlTumorZMIPSlabs():
    """Create slabs from the directories Controls and LargeTumor in
    input_image_root via createZMIPSlabs and put the results in
    controls and tumors subdirectories, respectively, of
    hardDrive_proj_root.

    """

    # Input data directories where mha/mhd and associated tre files are located
    controlInputDir = os.path.join(input_image_root, "Controls")
    tumorInputDir = os.path.join(input_image_root, "LargeTumor")

    # Output data directories
    controlOutputDir = os.path.join(hardDrive_proj_root, "controls")
    tumorOutputDir = os.path.join(hardDrive_proj_root, "tumors")

    # Process control files
    createZMIPSlabs('control', controlInputDir, controlOutputDir)

    # Process tumor files
    createZMIPSlabs('tumor', tumorInputDir, tumorOutputDir)

# Compute Training mask
def computeTrainingMask(expertSegMask, outputTrainingMask):
    """Compute a training mask from expertSegMask, written to
    outputTrainingMask.

    Note: A temporary file -- expertSegMask + "_skel.png" -- is
    created and then removed in the process.

    """

    skeletonFile = expertSegMask + "_skel.png"

    subprocess.call(["SegmentBinaryImageSkeleton",
                     expertSegMask,
                     skeletonFile])

    subprocess.call([
        "ImageMath", expertSegMask,
         # dilate vessel mask with kernel radius = 1
         "-M", "1", "1", "255", "0",
        # subtract vessel mask from dilated version to get vessel boundary
        "-a", "1", "-1", expertSegMask,
        # create training mask with vessel center-line (=255) and boundary (=128)
        "-a", "0.5", "255", skeletonFile,
        # write training mask
        "-W", "0", outputTrainingMask])

    os.remove(skeletonFile)

    # WARNING: Couldn't write to PNG using the implemented CLI
    # subprocess.call( ["ComputeTrainingMask",
    # expertSegMask,
    # outputTrainingMask,
    # "--notVesselWidth","1"] )


# save each of the z-mip slabs from an .mha file as .png files
def saveSlabs(mhaFile):
    """Extract each slab of mhaFile to a separate file.
    Input: $input/*.mha
    Output: $input/#_*.png, where # is the slab number

    """

    print 'saving slabs of %s' % mhaFile

    reader = itk.ImageFileReader.New(FileName = str(mhaFile))
    reader.Update()
    buf = itk.GetArrayFromImage(reader.GetOutput())

    # convert to [0, 255] range
    buf = 255.0 * (buf - buf.min()) / (buf.max() - buf.min())
    buf = buf.astype('uint8')

    # file names definition
    fileName = os.path.basename(os.path.splitext(mhaFile)[0])
    fileDir = os.path.dirname(os.path.abspath(mhaFile))

    # Iterate through each slab
    for i, slab in enumerate(buf):

        outputImage = os.path.join(
            fileDir, str(i) + "_" + fileName + ".png")

        skimage.io.imsave(outputImage, slab)


def splitData(name, inputDir, outputDir, trainOutputDir, testOutputDir):
    """Split the various outputs created from the image files in inputDir,
    which reside in outputDir, between trainOutputDir and
    testOutputDir, reorganizing them in the process.

    With an input file named *.mhd, the following outputs are moved
    into the following subdirectories in the destination folder:
    - *.mha: images
    - *_zslab.mha: images (also split into PNG slices)
    - *_zslab_points.mha: points
    - *_expert.mha: expert
    - *_zslab_expert.mha: expert (also split into PNG slices)

    """
    # Process files
    printSectionHeader('Splitting %s data into training and testing' % name)

    mhdFiles = glob.glob(os.path.join(inputDir, "*", "*.mhd"))

    for i, mhdFile in enumerate(mhdFiles):

        print("\n%s file %d/%d : %s" %
              (name, i + 1, len(mhdFiles), mhdFile))

        filePrefix = os.path.basename(os.path.splitext(mhdFile)[0])

        # Split equally for training and testing
        if i % 2 == 0:
            curOutputDir = trainOutputDir
        else:
            curOutputDir = testOutputDir

        # "suffix" is surround by fileName and '.mha'
        # "dir" is a subdirectory of curOutputDir
        suffixesAndDirsForCopying = [
            ('', 'images'), # input volume
            ('_zslab', 'images'), # z-mip slab volume
            ('_zslab_points', 'points'), # z-mip slab point map
            ('_expert', 'expert'), # expert volume
            ('_zslab_expert', 'expert'), # expert z-mip slab volume
        ]

        for suffix, dir in suffixesAndDirsForCopying:
            fileName = filePrefix + suffix + '.mha'
            utils.copy(os.path.join(outputDir, fileName),
                       os.path.join(curOutputDir, dir, fileName))

        # save slabs as pngs
        saveSlabs(os.path.join(curOutputDir, "images", filePrefix + "_zslab.mha"))
        saveSlabs(os.path.join(curOutputDir, "expert", filePrefix + "_zslab_expert.mha"))

# assign control and tumor volumes equally to training and testing
def splitControlTumorData():
    """Split the data created from the images in the directories Controls
    and LargeTumor in input_image_root via splitData and put the
    results in training and testing subdirectories of
    hardDrive_proj_root.

    """

    # Input data directories
    controlInputDir = os.path.join(input_image_root, "Controls")
    tumorInputDir = os.path.join(input_image_root, "LargeTumor")

    # Output data directories
    controlOutputDir = os.path.join(hardDrive_proj_root, "controls")
    tumorOutputDir = os.path.join(hardDrive_proj_root, "tumors")

    trainOutputDir = os.path.join(hardDrive_proj_root, "training")
    testOutputDir = os.path.join(hardDrive_proj_root, "testing")

    # Sanity checks
    utils.ensureDirectoryExists(trainOutputDir)
    utils.ensureDirectoryExists(testOutputDir)

    # Process control files
    splitData('control', controlInputDir, controlOutputDir, trainOutputDir, testOutputDir)

    # Process tumor files
    splitData('tumor', tumorInputDir, tumorOutputDir, trainOutputDir, testOutputDir)


# Extracts +ve (vessel center) and -ve (background) patches from image
def extractPatchesFromImage(rootDir, imageName, outputDir, patchListFile):
    """Convert an image to a set of patches.  Patches are sorted into
    "positive" and "negative" patches, i.e. those with and without a
    vessel at the center.  Positive patches have index 1, negative
    patches index 0.  Patch relative paths and their indices are
    written to the patchListFile file object.

    Input:
    - $rootDir/images/$imageName.png: The image to extract slices from
    - $rootDir/expert/$imageName_expert.png: The corresponding expert mask

    Output:
    - $outputDir/0/$imageName_$i_$j.png: Negative patches
    - $outputDir/1/$imageName_$i_$j.png: Positive patches

    """

    def writePatch(patchSetIndex, i, j):

        """Write out a patch centered at i,j to a correspondingly named file,
        using the given patch set index.  Create a corresponding entry
        in patchListFile.

        """
        psi = str(patchSetIndex)

        filename = os.path.join(
            psi, imageName + "_" + str(i) + "_" + str(j) + ".png")

        patchListFile.write(filename + " " + psi + "\n")

        skimage.io.imsave(os.path.join(outputDir, filename),
                          inputImage[i - w : i + w + 1, j - w : j + w + 1])

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

                writePatch(1, i, j)

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
        writePatch(0, i, j)

    # Pick rest of background (negative) patches away from vessel boundary
    numOtherBgndPatches = \
        int(math.ceil((1.0 - vessel_bnd_neg_frac) * numVesselPatches))

    selBgndInd = random.sample(bgndInd, numOtherBgndPatches)

    for [i, j] in selBgndInd:
        writePatch(0, i, j)

    print 'Number of positive patches = ', numVesselPatches


def createPatches(name, dataDir, patchListFile):
    """Create patch files from images in dataDir.

    Input:
    - $dataDir/images/*.png
    - $dataDir/expert/*_expert.png

    Output:
    - $dataDir/patches/{0,1}/*_$i_$j.png
    - $dataDir/patches/$patchListFile: List of patch files and patch index

    """
    printSectionHeader('Creating %s patches' % name)

    imageFiles = glob.glob(os.path.join(dataDir, "images", "*.png"))

    patchesDir = os.path.join(dataDir, "patches")
    for i in range(2):
        utils.ensureDirectoryExists(os.path.join(patchesDir, str(i)))

    with open(os.path.join(patchesDir, patchListFile), "w") as patchListFile:

        for i, imageFile in enumerate(imageFiles):

            print('\nCreating patches for %s file %d/%d - %s' %
                  (name, i + 1, len(imageFiles), imageFile))

            imageName = os.path.basename(os.path.splitext(imageFile)[0])

            extractPatchesFromImage(dataDir, imageName, patchesDir,
                                    patchListFile)

# convert train/test images to patches
def createTrainTestPatches():
    """Create training and testing patches from the training and testing
    subdirectories of the working data directory.  See createPatches
    for the exact files read and created.

    """

    # create training patches
    createPatches('training', dataDir=os.path.join(hardDrive_proj_root, "training"),
                  patchListFile="train.txt")

    # create testing patches
    createPatches('testing', dataDir=os.path.join(hardDrive_proj_root, "testing"),
                  patchListFile="val.txt")


def createLmdb(name, patchesDir, patchListFile, lmdbDir):
    """Create an LMDB instance in lmdbDir from the patches in patchesDir,
    indexed by patchListFile

    """
    caffe_tools_dir = os.path.join(caffe_root, "build", "tools")
    convert_imageset_exec = os.path.join(caffe_tools_dir, "convert_imageset")

    print('Creating %s lmdb ...\n' % name)

    if os.path.exists(lmdbDir):
        shutil.rmtree(lmdbDir)

    subprocess.call([convert_imageset_exec,
                     "--shuffle",
                     "--gray",
                     patchesDir,
                     os.path.join(patchesDir, patchListFile),
                     lmdbDir])

# create lmdb
def createTrainTestLmdb():
    """Create LMBD instances for the training and testing data.  For
    training and testing, take the patches created by
    createTrainTestPatches and create LMDB instances in the hard drive
    project directory titled Net_TrainData and Net_ValData,
    respectively.

    """

    printSectionHeader('Creating LMDBs for train and test data')

    # create training lmdb
    createLmdb('training', os.path.join(hardDrive_proj_root, 'training', 'patches/'),
               'train.txt', os.path.join(hardDrive_proj_root, "Net_TrainData"))

    # create testing lmdb
    createLmdb('testing', os.path.join(hardDrive_proj_root, 'testing', 'patches/'),
               'val.txt', os.path.join(hardDrive_proj_root, "Net_ValData"))


def printSectionHeader(title):

    print("\n" + "#" * (len(title) + 4))
    print('# ' + title)
    print("#" * (len(title) + 4) + "\n")


def run():

    # create z-mip slabs
    createControlTumorZMIPSlabs()

    # assign control and tumor volumes equally to training and testing
    # Note: this must be called after createZMIPSlabs()
    splitControlTumorData()

    # convert train/test images to patches
    createTrainTestPatches()

    # create lmdb database
    createTrainTestLmdb()


if __name__ == "__main__":
    run()
