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

import utils
from utils import script_params

# Append ITK libs
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Modules/ThirdParty/VNL/src/vxl/lib'))
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


# Smooth images
def smooth(inputImage, expertImage, outputImagePrefix):
    """Smooth inputImage and expertImage according to script_params.
    Output (where '*' stands for outputImagePrefix):
    - *_smooth.mha: Smoothed inputImage
    - *_smooth_expert.mha: Smoothed expertImage

    """
    outputImagePrefix = str(outputImagePrefix)

    smoothing_radius = script_params['SMOOTHING_RADIUS']

    reader = itk.ImageFileReader.New(FileName=str(inputImage))
    filter = itk.MedianImageFilter.New(reader.GetOutput(), Radius=smoothing_radius)
    writer = itk.ImageFileWriter.New(filter.GetOutput(),
                                     FileName=outputImagePrefix + "_smooth.mha",
                                     UseCompression=True)
    writer.Update()

    reader.SetFileName(str(expertImage))
    writer.SetFileName(outputImagePrefix + "_smooth_expert.mha")
    writer.Update()


def createSmoothedImagesForFile(mhdFile, outputDir):
    """Create smoothed image and related files corresponding to mhdFile in outputDir.
    Input:
    - $input/*.mhd: The image file header
    - $input/TRE/*.tre: The expert TRE file
    Output:
    - $output/*_expert.mha: The expert MHA file
    - : All output from smooth with $output/* as outputImagePrefix

    """
    fileName = os.path.basename(os.path.splitext(mhdFile)[0])
    fileDir = os.path.dirname(os.path.abspath(mhdFile))

    treFile = os.path.join(fileDir, "TRE", fileName + ".tre")
    expertSegFile = os.path.join(outputDir,
                                 fileName + "_expert.mha")

    # Process
    createExpertSegmentationMask(mhdFile, treFile, expertSegFile)

    smooth(mhdFile, expertSegFile,
           os.path.join(outputDir, fileName))

def createSmoothedImages(name, inputDir, outputDir):
    """Process all image files in immediate subdirectories of inputDir to
    correspondingly prefixed images in outputDir.  outputDir is
    created if it doesn't already exist.  The subdirectory structure
    is not replicated.

    See the documentation of createSmoothedImagesForFile for the exact
    files created.

    """
    # Sanity check
    utils.ensureDirectoryExists(outputDir)

    # Process files
    printSectionHeader('Creating smoothed images for %ss' % name)

    mhdFiles = glob.glob(os.path.join(inputDir, "*", "*.mhd"))

    for i, mhdFile in enumerate(mhdFiles):

        print("\n%s file %d/%d : %s" %
              (name, i + 1, len(mhdFiles), mhdFile))

        createSmoothedImagesForFile(mhdFile, outputDir)


# create z-mip slabs
def createControlTumorSmoothedImages():
    """Create smoothed images from the directories Controls and LargeTumor
    in input_data_root via createZMIPSlabs and put the results in
    controls and tumors subdirectories, respectively, of
    output_data_root.

    """

    # Input data directories where mha/mhd and associated tre files are located
    controlInputDir = os.path.join(input_data_root, "Controls")
    tumorInputDir = os.path.join(input_data_root, "LargeTumor")

    # Output data directories
    controlOutputDir = os.path.join(output_data_root, "controls")
    tumorOutputDir = os.path.join(output_data_root, "tumors")

    # Process control files
    createSmoothedImages('control', controlInputDir, controlOutputDir)

    # Process tumor files
    createSmoothedImages('tumor', tumorInputDir, tumorOutputDir)

# Compute Training mask
def computeTrainingMask(expertSegMask, outputTrainingMask):
    """Compute a training mask from expertSegMask, written to
    outputTrainingMask.  It takes on the following values:

    - 0: Background, other
    - 51: Vessel, other
    - 181: Vessel, center line
    - 255: Background, vessel edge

    Note: These values are arbitrary and were picked for the
    convenience of generating them.

    Note: Two intermediate files -- expertSegMask + "_skel.mha" and +
    "_dilated.mha" -- are created in the process.  The latter is
    removed afterwards.

    """

    skeletonFile = expertSegMask + "_skel.mha"

    subprocess.call(["SegmentBinaryImageSkeleton",
                     expertSegMask,
                     skeletonFile])

    dilated = expertSegMask + '_dilated.mha'

    subprocess.call([
        'ImageMath', expertSegMask,
        '-M', '1', '4', '1', '0',
        '-W', '0', dilated])

    subprocess.call([
        "ImageMath", dilated,
         # dilate vessel mask with kernel radius = 1
         "-M", "1", "1", "1", "0",
        # subtract vessel mask from dilated version to get vessel boundary
        "-a", "255", "-204", dilated,
        # create training mask with vessel center-line (=255) and boundary (=128)
        "-a", "1", "130", skeletonFile,
        # write training mask
        "-W", "0", outputTrainingMask])

    os.unlink(dilated)


def splitData(name, inputDir, outputDir, trainOutputDir, testOutputDir):
    """Split the various outputs created from the image files in inputDir,
    which reside in outputDir, between trainOutputDir and
    testOutputDir.

    With an input file named *.mhd, the following outputs are moved
    into the destination folder:
    - *_smooth.mha
    - *_expert.mha
    - *_smooth_expert.mha

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

        # "suffix" is surround by fileName+'_' and '.mha'
        for suffix in ['smooth', 'expert', 'smooth_expert']:
            fileName = filePrefix + '_' + suffix + '.mha'
            utils.symlink_entries_through(outputDir, curOutputDir, fileName)

# assign control and tumor volumes equally to training and testing
def splitControlTumorData():
    """Split the data created from the images in the directories Controls
    and LargeTumor in input_data_root via splitData and put the
    results in training and testing subdirectories of
    output_data_root.

    """

    # Input data directories
    controlInputDir = os.path.join(input_data_root, "Controls")
    tumorInputDir = os.path.join(input_data_root, "LargeTumor")

    # Output data directories
    controlOutputDir = os.path.join(output_data_root, "controls")
    tumorOutputDir = os.path.join(output_data_root, "tumors")

    trainOutputDir = os.path.join(output_data_root, "training")
    testOutputDir = os.path.join(output_data_root, "testing")

    # Sanity checks
    utils.ensureDirectoryExists(trainOutputDir)
    utils.ensureDirectoryExists(testOutputDir)

    # Process control files
    splitData('control', controlInputDir, controlOutputDir, trainOutputDir, testOutputDir)

    # Process tumor files
    splitData('tumor', tumorInputDir, tumorOutputDir, trainOutputDir, testOutputDir)


def dict_to_list(d):
    """Convert a dict whose keys are the ints 0..N-1 into a list of length
    N such that l[x] == d[x].  In the process, check that the keys are
    indeed such a range.

    """
    l = [None]*len(d)
    try:
        for k, v in d.items():
            l[k] = v
    except IndexError:
        raise ValueError("Argument does not have a complete set of indices!")
    return l

# Extracts +ve (vessel center) and -ve (background) patches from image
def extractPatchesFromImageGenerator(rootDir, imageName):
    # TODO update docstring
    """Convert an image to a set of patches.  Patches are sorted into
    "positive" and "negative" patches, i.e. those with and without a
    vessel at the center.  Positive patches have index 1, negative
    patches index 0.  Patch relative paths and their indices are
    yielded with the image data as part of a generator.

    Input:
    - $rootDir/$imageName.mha: The image to extract slices from
    - $rootDir/$imageName_expert.mha: The corresponding expert mask

    Output:
    - ("0/$imageName/$i_$j.png", 0, image): Negative patches
    - ("1/$imageName/$i_$j.png", 1, image): Positive patches

    """

    def patchEntry(patchSetIndex, coords):
        # TODO document
        filename = os.path.join(
            str(patchSetIndex), imageName, '_'.join(map(str, coords)) + ".png")

        image = utils.extractPatch(inputImage, coords)

        return filename, patchSetIndex, image

    # patch/window radius
    w = script_params['PATCH_RADIUS']

    total_pos_patches = script_params['POSITIVE_PATCHES_PER_INPUT_FILE']
    total_neg_patches = script_params['NEGATIVE_TO_POSITIVE_RATIO'] * total_pos_patches

    num_patch_types = 4
    vessel_ctl_pos, other_pos, vessel_bnd_neg, other_neg = range(num_patch_types)

    patch_index = np.array(dict_to_list({
        vessel_ctl_pos: 1,
        other_pos: 1,
        vessel_bnd_neg: 0,
        other_neg: 0,
    }))

    frac = np.array(dict_to_list({
        vessel_ctl_pos: script_params['POSITIVES_NEAR_VESSEL_CENTERLINE'],
        other_pos: script_params['OTHER_POSITIVES'],
        vessel_bnd_neg: script_params['NEGATIVES_NEAR_VESSEL_BOUNDARY'],
        other_neg: script_params['OTHER_NEGATIVES'],
    }))

    for i in range(2):
        if abs(frac[patch_index == i].sum() - 1.0) > 1e-9:
            raise ValueError("{} patch fraction sum must be 1.0, is {}".format(
                "Positive" if i else "Negative", patch_frac_sum))

    # read input image
    inputImageFile = os.path.join(rootDir, imageName + '.mha')
    inputImageReader = itk.ImageFileReader.New(FileName=str(inputImageFile))
    inputImageReader.Update()
    inputImage = itk.GetArrayFromImage(inputImageReader.GetOutput())

    # read expert segmented label mask
    expertSegFile = os.path.join(rootDir, imageName + '_expert.mha')
    trainingMaskFile = os.path.join(rootDir, imageName + '_train.mha')

    computeTrainingMask(expertSegFile, trainingMaskFile)

    trainingMaskReader = itk.ImageFileReader.New(FileName=str(trainingMaskFile))
    trainingMaskReader.Update()
    trainingMask = itk.GetArrayFromImage(trainingMaskReader.GetOutput())

    # Iterate through expert mask and find pos/neg patches

    # Slice that we want, which excludes edge pixels
    s = np.s_[w:-w or None]
    trainingMaskMid = trainingMask[(s,) * 3]

    mask = dict_to_list({
        # Vessel centerline pixel (positive)
        vessel_ctl_pos: trainingMaskMid == 181,
        # Other vessel pixel (positive)
        other_pos: trainingMaskMid == 51,
        # Vessel bound pixel (negative)
        vessel_bnd_neg: trainingMaskMid == 255,
        # Other background pixel (negative)
        other_neg: trainingMaskMid == 0,
    })

    # Linear, flat indices
    indices = [np.where(m.reshape(-1))[0] for m in mask]
    # Desired number of each type of patch (not necessarily a whole number)
    desiredFractionalPatches = frac * np.where(patch_index, total_pos_patches, total_neg_patches)
    availablePatches = np.array([len(ind) for ind in indices])
    availableFraction = availablePatches.astype(float) / desiredFractionalPatches
    # The largest fraction we can take, capped at 1
    minAvailableFraction = min(availableFraction.min(), 1.)
    if minAvailableFraction != 1.:
        print("WARNING: Too few of a desired patch type; " +
              "scaling all types down by {}% accordingly".format(minAvailableFraction * 100))
    numPatches = np.ceil(desiredFractionalPatches * minAvailableFraction - 1e-9).astype(int)

    for key, label in enumerate(patch_index):
        selInd = random.sample(indices[key], numPatches[key])
        selInd = np.transpose(np.unravel_index(selInd, trainingMaskMid.shape)) + w
        for coords in selInd:
            yield patchEntry(label, coords)


def createPatchesGenerator(name, dataDir):
    # TODO update documentation
    """Create patch files from images in dataDir.

    Input:
    - $dataDir/images/*.png
    - $dataDir/expert/*_expert.png

    Output:
    - $dataDir/patches/{0,1}/*_$i_$j.png
    - $dataDir/patches/$patchListFile: List of patch files and patch index

    """
    printSectionHeader('Creating %s patches' % name)

    imageFiles = glob.glob(os.path.join(dataDir, "*_smooth.mha"))

    for i, imageFile in enumerate(imageFiles):

        print('\nCreating patches for %s file %d/%d - %s' %
              (name, i + 1, len(imageFiles), imageFile))

        imageName = os.path.basename(os.path.splitext(imageFile)[0])

        for result in extractPatchesFromImageGenerator(dataDir, imageName):
            yield result

def createDB(name, dataDir, dbDir):
    # TODO update documentation
    """Create a database in dbDir from the patches in patchesDir, indexed
    by patchListFile

    """
    print('Creating %s DB ...\n' % name)

    if os.path.exists(dbDir):
        shutil.rmtree(dbDir)

    utils.ensureDirectoryExists(dbDir)

    db = utils.open_sqlite3_db(dbDir)
    db.execute('''create table "temp"."PatchesUnshuffled" (
        "filename" text,
        "patch_index" integer,
        -- Interpreted as a square C-major-order uint8 array with each
        -- dimension (2 * PATCH_RADIUS + 1) and extra final dimension of size 3
        "image_data" blob
    )''')

    def generate_data():
        patches = createPatchesGenerator(name, dataDir)
        for i, (filename, patch_index, image_data) in enumerate(patches):
            yield filename, patch_index, buffer(image_data.copy())
            if i % 1000 == 0:
                print("{} patches processed".format(i))

    db.executemany('''insert into "PatchesUnshuffled" values (?, ?, ?)''',
                   generate_data())

    print("Shuffling patches")

    db.execute('''create table "main"."Patches" as
                  select * from "PatchesUnshuffled"
                  order by random()''')

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

    # create z-mip slabs
    createControlTumorSmoothedImages()

    # assign control and tumor volumes equally to training and testing
    # Note: this must be called after createControlTumorSmoothedImages()
    splitControlTumorData()

    # create database
    createTrainTestDB()


if __name__ == "__main__":
    run()
