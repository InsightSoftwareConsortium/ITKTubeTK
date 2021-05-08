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
import numpy as np

import utils

# Append ITK libs
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build',
                             'Modules/ThirdParty/VNL/src/vxl/lib'))
import itk

# Define paths
script_params = json.load(open('params.json'))

output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

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
    input_data_root via createZMIPSlabs and put the results in
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
    - $outputDir/0/$imageName/$i_$j.png: Negative patches
    - $outputDir/1/$imageName/$i_$j.png: Positive patches

    """

    def writePatch(patchSetIndex, i, j):

        """Write out a patch centered at i,j to a correspondingly named file,
        using the given patch set index.  Create a corresponding entry
        in patchListFile.

        """
        psi = str(patchSetIndex)

        filename = os.path.join(
            psi, imageName, str(i) + "_" + str(j) + ".png")

        patchListFile.write(filename + " " + psi + "\n")

        skimage.io.imsave(os.path.join(outputDir, filename),
                          inputImage[i - w : i + w + 1, j - w : j + w + 1])

    for i in range(2):
        utils.ensureDirectoryExists(os.path.join(outputDir, str(i), imageName))

    # patch/window radius
    w = script_params['PATCH_RADIUS']

    total_pos_patches = script_params['POSITIVE_PATCHES_PER_INPUT_FILE'] / script_params['NUM_SLABS']
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
    inputImageFile = os.path.join(rootDir, "images", imageName + '.png')
    inputImage = skimage.io.imread(inputImageFile)

    # read expert segmented label mask
    expertSegFile = os.path.join(rootDir, "expert", imageName + '_expert.png')
    trainingMaskFile = os.path.join(rootDir, "expert", imageName + '_train.png')

    computeTrainingMask(expertSegFile, trainingMaskFile)

    expertSeg = skimage.io.imread(expertSegFile)
    trainingMask = skimage.io.imread(trainingMaskFile)

    # Iterate through expert mask and find pos/neg patches

    # Slice that we want, which excludes edge pixels
    s = np.s_[w:-w-1]
    s = (s, s)
    trainingMaskMid = trainingMask[s]
    expertSegMid = expertSeg[s]

    vessel_ctl_pos_mask = trainingMaskMid > 0.6 * 255
    mask = dict_to_list({
        # Vessel centerline pixel (positive)
        vessel_ctl_pos: vessel_ctl_pos_mask,
        # Other vessel pixel (positive)
        other_pos: (expertSegMid > 0) & ~vessel_ctl_pos_mask,
        # Vessel bound pixel (negative)
        vessel_bnd_neg: (trainingMaskMid > 0) & ~vessel_ctl_pos_mask,
        # Other background pixel (negative)
        other_neg: (expertSegMid == 0) & (trainingMaskMid == 0),
    })

    indices = [np.array(np.where(m)).T + w for m in mask]
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
        for i, j in selInd:
            writePatch(label, i, j)


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
    createPatches('training', dataDir=os.path.join(output_data_root, "training"),
                  patchListFile="train.txt")

    # create testing patches
    createPatches('testing', dataDir=os.path.join(output_data_root, "testing"),
                  patchListFile="val.txt")


def createDB(name, patchesDir, patchListFile, dbDir):
    """Create a database in dbDir from the patches in patchesDir, indexed
    by patchListFile

    """
    print('Creating %s DB ...\n' % name)

    if os.path.exists(dbDir):
        shutil.rmtree(dbDir)

    utils.ensureDirectoryExists(dbDir)

    patchListFile = os.path.join(patchesDir, patchListFile)

    db = utils.open_sqlite3_db(dbDir)
    db.execute('''create table "Patches" (
        "filename" text,
        "patch_index" integer,
        -- Interpreted as a square C-major-order uint8 array with each
        -- dimension (2 * PATCH_RADIUS + 1)
        "image_data" blob
    )''')

    lines = list(open(patchListFile))
    random.shuffle(lines)

    def generate_data():
        for i, l in enumerate(lines):
            filename, patch_index = l.rstrip().rsplit(' ', 1)
            patch_index = int(patch_index)
            image_data = buffer(skimage.io.imread(os.path.join(patchesDir, filename)).copy())

            yield filename, patch_index, image_data

            if i % 1000 == 0:
                print("{} patches processed".format(i))

    db.executemany('''insert into "Patches" values (?, ?, ?)''',
                   generate_data())

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
    createDB('training', os.path.join(output_data_root, 'training', 'patches/'),
             'train.txt', os.path.join(output_data_root, "Net_TrainData"))

    # create testing DB
    createDB('testing', os.path.join(output_data_root, 'testing', 'patches/'),
             'val.txt', os.path.join(output_data_root, "Net_ValData"))


def printSectionHeader(title):

    print("\n" + "#" * (len(title) + 4))
    print('# ' + title)
    print("#" * (len(title) + 4) + "\n")


def run():

    # create z-mip slabs
    #createControlTumorZMIPSlabs()

    # assign control and tumor volumes equally to training and testing
    # Note: this must be called after createZMIPSlabs()
    #splitControlTumorData()

    # convert train/test images to patches
    createTrainTestPatches()

    # create database
    createTrainTestDB()


if __name__ == "__main__":
    run()
