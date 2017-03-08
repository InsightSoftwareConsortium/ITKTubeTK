#!/usr/bin/python

###########################################################################
# PreProcessing.py :
#
# Iterate through the expert labelmap and create 65x65 patches around the
# central pixel. All positive pixels are used as positives input cases.
# The same amount of negatives is randomly picked. For each input patch,
# the corresponding filename and expected output are written to a text file
# and will be used later to create the database.
#
###########################################################################

import os
import glob
import json

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

import random
import math


# Create image set and save expected output
def createImgSet(expertImg, inputImg, filenamePrefix, fileOutputDir, textFile):

    w = 32  # Patch size
    count = 0  # Input count
    negativeImageIndex = []  # Whole image index giving negative output
    negativeIndex = []  # Vessel bound index giving negative output

    # Write filename and expected output
    textFile = open(textFile, "a")
    textFile.truncate()  # Erase file

    resample = 1  # WARNING: resample pixels to reduce training size for Debug

    # Iterate through the expert label map
    for i in range(0, expertImg.shape[0], resample):
        for j in range(0, expertImg.shape[1], resample):

            if j > w and j + w + 1 < inputImg.shape[1]:
                if i > w and i + w + 1 < inputImg.shape[0]:

                    # Centerline pixel (positive)
                    if expertImg[i, j] > 0.5:
                        count += 1
                        filename = os.path.join(
                            "1", filenamePrefix + "_" + str(i) + "_" + str(j) + ".png")
                        textFile.write(filename + " " + str(1) + "\n")
                        plt.imsave(os.path.join(fileOutputDir, filename),
                                   inputImg[i - w:i + w + 1, j - w:j + w + 1], cmap='Greys_r')

                    # Vessel bound pixel (negative)
                    elif expertImg[i, j] > 0:
                        negativeIndex.append([i, j])

                    # Background pixel (negative)
                    else:
                        negativeImageIndex.append([i, j])

    # Pick random negatives from vessel bound
    rndmNegativeInd = random.sample(negativeIndex, int(math.ceil(0.8 * count)))
    for [i, j] in rndmNegativeInd:
        filename = os.path.join("0", filenamePrefix +
                                "_" + str(i) + "_" + str(j) + ".png")
        textFile.write(filename + " " + str(0) + "\n")
        plt.imsave(os.path.join(fileOutputDir, filename),
                   inputImg[i - w:i + w + 1, j - w:j + w + 1], cmap='Greys_r')

    # Pick random negatives from the entire image
    rndmNegativeImageInd = random.sample(
        negativeImageIndex, int(math.ceil(0.2 * count)))
    for [i, j] in rndmNegativeImageInd:
        filename = os.path.join("0", filenamePrefix +
                                "_" + str(i) + "_" + str(j) + ".png")
        textFile.write(filename + " " + str(0) + "\n")
        plt.imsave(os.path.join(fileOutputDir, filename),
                   inputImg[i - w:i + w + 1, j - w:j + w + 1], cmap='Greys_r')

    textFile.close()
    print(count)

########
# Main #
########

# Path variable
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']
proj_rel_path = script_params['PROJECT_REL_PATH']

caffe_proj_root = os.path.join(caffe_root, "data", proj_rel_path)
hardDrive_proj_root = os.path.join(hardDrive_root, proj_rel_path)

trainDataDir = os.path.join(hardDrive_proj_root, "training")
valDataDir = os.path.join(hardDrive_proj_root, "testing")

# Text file
trainFilename = os.path.join(caffe_proj_root, "train.txt")
trainFile = open(trainFilename, "w+")
trainFile.truncate()  # Erase file
trainFile.close()

valFilename = os.path.join(caffe_proj_root, "val.txt")
valFile = open(valFilename, "w+")
valFile.truncate()  # Erase file
valFile.close()

# Output patches directories
trainFileOutputDir = os.path.join(trainDataDir, "out")
if not os.path.exists(trainFileOutputDir):
    os.mkdir(trainFileOutputDir)

for label in range(2):
    curLabelOutputDir = os.path.join(trainFileOutputDir, str(label))
    if not os.path.exists(curLabelOutputDir):
        os.mkdir(curLabelOutputDir)

valFileOutputDir = os.path.join(valDataDir, "out")
if not os.path.exists(valFileOutputDir):
    os.mkdir(valFileOutputDir)

for label in range(2):
    curLabelOutputDir = os.path.join(valFileOutputDir, str(label))
    if not os.path.exists(curLabelOutputDir):
        os.mkdir(curLabelOutputDir)

# Images directories
trainExpertDir = os.path.join(trainDataDir, "expert")
trainImgDir = os.path.join(trainDataDir, "images")

valExpertDir = os.path.join(valDataDir, "expert")
valImgDir = os.path.join(valDataDir, "images")

# Create train set
trainImages = glob.glob(os.path.join(trainImgDir, "*.png"))

for trainImage in trainImages:

    print(trainImage)

    # Get image ID
    trainImagePrefix = os.path.basename(os.path.splitext(trainImage)[0])

    # Set filename
    trainExpertFile = os.path.join(
        trainExpertDir, trainImagePrefix + "_expert.png")
    trainImageFile = os.path.join(trainImgDir, trainImagePrefix + ".png")

    # Load images
    trainExpert = mpimg.imread(trainExpertFile)

    # print trainExpert.shape
    # trainExpert=trainExpert[:,:,0]

    trainImg = mpimg.imread(trainImageFile)
    # trainImg=trainImg[:,:,0]

    # Write images and text files
    createImgSet(trainExpert, trainImg, trainImagePrefix,
                 trainFileOutputDir, trainFilename)

# Create validation set
valImages = glob.glob(os.path.join(valImgDir, "*.png"))
for valImage in valImages:

    print(valImage)

    # Get image ID
    valImagePrefix = os.path.basename(os.path.splitext(valImage)[0])

    # Set filename
    valExpertFilename = os.path.join(
        valExpertDir, valImagePrefix + "_expert.png")
    valImgFilename = os.path.join(valImgDir, valImagePrefix + ".png")

    # Load images
    valExpert = mpimg.imread(valExpertFilename)

    # valExpert=valExpert[:,:,0]
    valImg = mpimg.imread(valImgFilename)
    # valImg=valImg[:,:,0]

    # Write images and text files
    createImgSet(valExpert, valImg, valImagePrefix,
                 valFileOutputDir, valFilename)
