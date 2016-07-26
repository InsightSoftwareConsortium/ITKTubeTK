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

import os, glob, sys

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

import numpy as np

import random, math

# Create image set and save expected output
def createImgSet( expertImg, inputImg, filenamePrefix, fileOutputDir,
  textFilename ):
  w = 32; #Patch size
  count = 0; #Input count
  negativeImageIndex = [] #Whole image index giving negative output
  negativeIndex = [] #Vessel bound index giving negative output

  # Write filename and expected output
  textFile = open(textFilename, "a")
  textFile.truncate() # Erase file

  resample = 1;#WARNING: resample pixels to reduce training size for Debug

  # Iterate through the expert label map
  for i in range(0,expertImg.shape[0],resample):
    for j in range(0,expertImg.shape[1],resample):
      if j>w and j+w+1<inputImg.shape[1] :
        if i>w and i+w+1<inputImg.shape[0]:
              # Centerline pixel (positive)
              if expertImg[i,j] > 0.5:
                count = count + 1
                filename = filenamePrefix + "_" + str(i) + "_" + str(j) +".png"
                textFile.write(filename + " " + str(1) + "\n")
                plt.imsave(fileOutputDir + filename, inputImg[i-w:i+w+1,j-w:j+w+1],cmap='Greys_r')
              # Vessel bound pixel (negative)
              elif expertImg[i,j] > 0:
                negativeIndex.append([i,j])
              # Background pixel (negative)
              else :
                negativeImageIndex.append([i,j])

  # Pick random negatives from vessel bound
  rndmNegativeInd = random.sample(negativeIndex, int(math.ceil(0.8*count)))
  for [i,j] in rndmNegativeInd :
    filename = filenamePrefix + "_" + str(i) + "_" + str(j) + ".png"
    textFile.write(filename + " " + str(0) + "\n")
    plt.imsave(fileOutputDir + filename, inputImg[i-w:i+w+1,j-w:j+w+1],cmap='Greys_r')
  # Pick random negatives from the entire image
  rndmNegativeImageInd = random.sample(negativeImageIndex, int(math.ceil(0.2*count)))
  for [i,j] in rndmNegativeImageInd :
    filename = filenamePrefix + "_" + str(i) + "_" + str(j) + ".png"
    textFile.write(filename + " " + str(0) + "\n")
    plt.imsave(fileOutputDir + filename, inputImg[i-w:i+w+1,j-w:j+w+1],cmap='Greys_r')

  textFile.close()
  print count

########
# Main #
########
# Path variable
hardDrive_root = "/media/lucas/krs0014/"
caffe_root = "./"

# Text file
trainFilename = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/train.txt"
trainFile = open(trainFilename, "w+")
trainFile.truncate() # Erase file
trainFile.close()

valFilename = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/val.txt"
valFile = open(valFilename, "w+")
valFile.truncate() # Erase file
valFile.close()

# Output patches directories
trainFileOutputDir= hardDrive_root + "SegmentVesselsUsingNeuralNetworks/training/out/"
valFileOutputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/out/"

# Images directories
trainExpertDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/training/expert/"
trainImgDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/training/images/"

valExpertDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/expert/"
valImgDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/images/"

# Create train set
trainImages = glob.glob( os.path.join( trainImgDir, "*.png" ) )
for trainImage in trainImages:
  print trainImage
  # Get image ID
  trainImagePrefix = os.path.basename(os.path.splitext(trainImage)[0])
  # Set filename
  trainExpertFilename = trainExpertDir + trainImagePrefix + "_expert.png"
  trainImgFilename = trainImgDir + trainImagePrefix + ".png"
  # Load images
  trainExpert=mpimg.imread(trainExpertFilename)
  #print trainExpert.shape
  #trainExpert=trainExpert[:,:,0]
  trainImg=mpimg.imread(trainImgFilename)
  #trainImg=trainImg[:,:,0]

  # Write images and text files
  createImgSet( trainExpert, trainImg, trainImagePrefix, trainFileOutputDir, trainFilename )

#Create validation set
valImages = glob.glob( os.path.join( valImgDir, "*.png" ) )
for valImage in valImages:
  print valImage
  # Get image ID
  valImagePrefix = os.path.basename(os.path.splitext(valImage)[0])
  # Set filename
  valExpertFilename = valExpertDir + valImagePrefix + "_expert.png"
  valImgFilename = valImgDir + valImagePrefix + ".png"
  # Load images
  valExpert=mpimg.imread(valExpertFilename)
  #valExpert=valExpert[:,:,0]
  valImg=mpimg.imread(valImgFilename)
  #valImg=valImg[:,:,0]

  # Write images and text files
  createImgSet( valExpert, valImg, valImagePrefix, valFileOutputDir, valFilename )
