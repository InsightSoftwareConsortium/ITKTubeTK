#!/usr/bin/python

###########################################################################
# PrepareTrainingData.py :
#
# Prepares training data for CNN
#
###########################################################################

import glob
import math
import os
import random
import subprocess
import shutil
import sys

import skimage.io
import numpy as np

import deploy
import utils
from utils import script_params

import itk

# Define paths
output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

# Create segmentation mask from tre file
def createExpertSegmentationMask(inputImageFile, treFile, outputExpertSegFile):

    # convert tre file to image
    subprocess.call(["ConvertTubesToImage", "-r",
                     inputImageFile,
                     treFile,
                     outputExpertSegFile])


def createPreppedImagesForFile(mhdFile, outputDir):
    """Create preprocessed image and related files corresponding to mhdFile in outputDir.
    Input:
    - $input/*.mh[ad]: The image file header
    - $input/TRE/*.tre: The expert TRE file
    Output:
    - $output/*_prepped_expert.mha: The expert MHA file matching the prepped output
    - : All output from prep with $output/* as outputImagePrefix

    """
    fileName = os.path.basename(os.path.splitext(mhdFile)[0])
    fileDir = os.path.dirname(os.path.abspath(mhdFile))

    treFile = os.path.join(fileDir, "TRE", fileName + ".tre")
    expertSegFile = os.path.join(outputDir,
                                 fileName + "_prepped_expert.mha")

    # Process
    _, prepped_mhd_file = deploy.prep(mhdFile, outputDir)

    createExpertSegmentationMask(prepped_mhd_file, treFile, expertSegFile)

def createPreppedImagesForType(name, inputDir, outputDir):
    """Process all image files in immediate subdirectories of inputDir to
    correspondingly prefixed images in outputDir.  outputDir is
    created if it doesn't already exist.  The subdirectory structure
    is not replicated.

    See the documentation of createPreppedImagesForFile for the exact
    files created.

    """
    # Sanity check
    utils.ensureDirectoryExists(outputDir)

    # Process files
    printSectionHeader('Preprocessing images for %s' % name)

    mhdFiles = glob.glob(os.path.join(inputDir, script_params['TYPE_SUBDIR_STRUCTURE'], "*.mh[ad]"))

    for i, mhdFile in enumerate(mhdFiles):

        print("\n%s file %d/%d : %s" %
              (name, i + 1, len(mhdFiles), mhdFile))

        createPreppedImagesForFile(mhdFile, outputDir)


def createPreppedImages():
    """Create preprocessed images from the directories in input_data_root
    that are the keys of TYPES via createPreppedImagesForType and put
    the results in subdirectories of output_data_root named as the
    corresponding values.

    """
    for input_dir, output_dir in script_params['TYPES'].items():
        createPreppedImagesForType(
            output_dir,
            # Input data directory where mha/mhd and associated tre files are located
            os.path.join(input_data_root, input_dir),
            os.path.join(output_data_root, output_dir),
        )

def splitDataForType(name, inputDir, outputDir, trainOutputDir, testOutputDir):
    """Split the various outputs created from the image files in inputDir,
    which reside in outputDir, between trainOutputDir and
    testOutputDir.

    With an input file named *.mh[ad], the following outputs are moved
    into the destination folder:
    - *_prepped.mha
    - *_resampled.mha (If RESAMPLED_SPACING is a number)
    - *_prepped_expert.mha

    """
    # Process files
    printSectionHeader('Splitting %s data into training and testing' % name)

    mhdFiles = glob.glob(os.path.join(inputDir, script_params['TYPE_SUBDIR_STRUCTURE'], "*.mh[ad]"))

    for i, mhdFile in enumerate(mhdFiles):

        print("\n%s file %d/%d : %s" %
              (name, i + 1, len(mhdFiles), mhdFile))

        filePrefix = os.path.basename(os.path.splitext(mhdFile)[0])

        # Split equally for training and testing
        if i % 2 == 0:
            curOutputDir = trainOutputDir
        else:
            curOutputDir = testOutputDir

        for path in glob.glob(os.path.join(outputDir, filePrefix + '_*.mha')):
            fileName = os.path.basename(path)
            utils.symlink_entries_through(outputDir, curOutputDir, fileName)

# assign volumes equally to training and testing
def splitData():
    """Split the data created from the images in the directories in
    input_data_root named by the keys of TYPES via splitData and put
    the results in training and testing subdirectories of
    output_data_root.

    """
    trainOutputDir = os.path.join(output_data_root, "training")
    testOutputDir = os.path.join(output_data_root, "testing")

    # Sanity checks
    utils.ensureDirectoryExists(trainOutputDir)
    utils.ensureDirectoryExists(testOutputDir)

    for input_dir, output_dir in script_params['TYPES'].items():
        splitDataForType(
            output_dir,
            os.path.join(input_data_root, input_dir),
            os.path.join(output_data_root, output_dir),
            os.path.join(output_data_root, "training"),
            os.path.join(output_data_root, "testing"),
        )


# Extracts +ve (vessel center) and -ve (background) patches from image
def extractPatchesFromImageGenerator(rootDir, imageName):
    """Convert an image to a set of patches.  Patches are sorted into
    "positive" and "negative" patches, i.e. those with and without a
    vessel at the center.  Positive patches have index 1, negative
    patches index 0.  A path-like string containing patch information
    is yielded with the image data for each patch as part of a
    generator.

    Input:
    - $rootDir/$imageName.mha: The image to extract slices from
    - $rootDir/$imageName_expert.mha: The corresponding expert mask

    Output:
    - ("0/$imageName/$i_$j_$k.png", 0, image): Negative patches
    - ("1/$imageName/$i_$j_$k.png", 1, image): Positive patches

    """

    def patchEntry(patchSetIndex, coords):
        """Return a pseudo-path, the given index, and the patch corresponding
        to the given coordinates of the input image.

        """
        filename = os.path.join(
            str(patchSetIndex), imageName, '_'.join(map(str, coords)) + ".png")

        image = utils.extractPatch(inputImagePadded, coords)

        return filename, patchSetIndex, image

    # patch/window radius
    w = script_params['PATCH_RADIUS']

    total_pos_patches = script_params['POSITIVE_PATCHES_PER_INPUT_FILE']
    total_neg_patches = script_params['NEGATIVE_TO_POSITIVE_RATIO'] * total_pos_patches

    # read input image
    inputImageFile = os.path.join(rootDir, imageName + '.mha')
    inputImageReader = itk.ImageFileReader.New(FileName=str(inputImageFile))
    inputImageReader.Update()
    inputImage = itk.GetArrayViewFromImage(inputImageReader.GetOutput())
    inputImagePadded = utils.pad(inputImage)

    # read expert segmented label mask
    expertSegFile = os.path.join(rootDir, imageName + '_expert.mha')

    expertMaskImage = itk.imread(str(expertSegFile))
    expertMask = itk.GetArrayViewFromImage(expertMaskImage)

    lbm = deploy.locally_brightest_mask(inputImage)
    trainingMask = np.where(lbm, expertMask, -1)

    # Iterate through expert mask and find pos/neg patches
    mask = [trainingMask == 0, trainingMask == 1]

    # Linear, flat indices
    indices = [np.where(m.reshape(-1))[0] for m in mask]
    # Desired number of each type of patch (not necessarily a whole number)
    desiredFractionalPatches = np.array([total_neg_patches, total_pos_patches])
    availablePatches = np.array([len(ind) for ind in indices])
    availableFraction = availablePatches.astype(float) / desiredFractionalPatches
    # The largest fraction we can take, capped at 1
    minAvailableFraction = min(availableFraction.min(), 1.)
    if minAvailableFraction != 1.:
        print("WARNING: Too few of a desired patch type; " +
              "scaling all types down by {}% accordingly".format(minAvailableFraction * 100))
    numPatches = np.ceil(desiredFractionalPatches * minAvailableFraction - 1e-9).astype(int)

    for label, (ind, n) in enumerate(zip(indices, numPatches)):
        selInd = random.sample(ind, n)
        selInd = np.transpose(np.unravel_index(selInd, trainingMask.shape))
        for coords in selInd:
            yield patchEntry(label, coords)


def createPatchesGenerator(name, dataDir):
    """Create a generator yielding patches from the images in dataDir.

    Input:
    - $dataDir/*_prepped.mha
    - $dataDir/*_prepped_expert.mha

    Output:
    - ("{0,1}/*/$i_$j_$k.png", {0,1}, image)

    """
    printSectionHeader('Creating %s patches' % name)

    imageFiles = glob.glob(os.path.join(dataDir, "*_prepped.mha"))

    for i, imageFile in enumerate(imageFiles):

        print('\nCreating patches for %s file %d/%d - %s' %
              (name, i + 1, len(imageFiles), imageFile))

        imageName = os.path.basename(os.path.splitext(imageFile)[0])

        for result in extractPatchesFromImageGenerator(dataDir, imageName):
            yield result

def createDB(name, dataDir, dbDir):
    """Create a database in dbDir from the images in dataDir."""
    print('Creating %s DB ...\n' % name)

    if os.path.exists(dbDir):
        shutil.rmtree(dbDir)

    utils.ensureDirectoryExists(dbDir)

    db = utils.open_sqlite3_db(dbDir)
    db.execute('''create table "temp"."PatchesUnshuffled" (
        "filename" text,
        "patch_index" integer,
        -- Interpreted as a square C-major-order float16 array with each
        -- dimension (2 * PATCH_RADIUS + 1) and extra final dimension of size 3
        "image_data" blob
    )''')

    def generate_data():
        patches = createPatchesGenerator(name, dataDir)
        for i, (filename, patch_index, image_data) in enumerate(patches):
            yield filename, patch_index, buffer(image_data.astype(np.float16))
            if i % 1000 == 0:
                print("{} patches processed".format(i))

    db.executemany('''insert into "PatchesUnshuffled" values (?, ?, ?)''',
                   generate_data())

    print("Shuffling patches")

    with utils.choice(db, "PatchesUnshuffled") as select:
        db.execute('create table "Patches" as ' + select)

    db.execute('''drop table "PatchesUnshuffled"''')

    db.commit()
    db.close()

def createTrainTestDB():
    """Create databases for the training and testing data.

    For training and testing, take the patches created by
    createTrainTestPatches and create a database in the output
    directory titled Net_TrainData and Net_ValData, respectively.

    """

    printSectionHeader('Creating DBs for train and test data')

    # create training DB
    createDB('training', os.path.join(output_data_root, 'training'),
             os.path.join(output_data_root, "Net_TrainData"))

    # create testing DB
    createDB('testing', os.path.join(output_data_root, 'testing'),
             os.path.join(output_data_root, "Net_ValData"))


def printSectionHeader(title):

    print("\n" + "#" * (len(title) + 4))
    print('# ' + title)
    print("#" * (len(title) + 4) + "\n")


def run():

    # preprocess images
    createPreppedImages()

    # assign volumes equally to training and testing
    # Note: this must be called after createPreppedImages()
    splitData()

    # create database
    createTrainTestDB()


if __name__ == "__main__":
    run()
