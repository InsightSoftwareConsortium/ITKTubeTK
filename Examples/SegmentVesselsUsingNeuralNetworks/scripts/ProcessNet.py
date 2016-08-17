#!/usr/bin/python

###########################################################################
# ProcessNet.py
###########################################################################

from __future__ import division
import os, sys

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pylab import *

sys.path.insert(0, caffe_root + 'python') #Append pycaffe to sys.path

import caffe
from caffe import layers as L, params as P

# Path variable
caffe_root = "./"  # this file should be run from {caffe_root} (otherwise change this line)

## Process one input image
def segmentImage( net, inputFilename, outputFilename ):
  img=mpimg.imread( inputFilename )
  # WARNING : Useful Debug plot and print to check data pixel value and shape
  #img=img.astype(float)
  #plt.imshow(img[:,:,0], cmap='gray')
  #plt.show()
  img=img[:,:,0]
  #print img[250:260,250]
  #plt.imshow(img, cmap='gray')
  #plt.show()

  w = 32; # Patch size
  index = [] # Save pixel index to reconstruct output
  inputImg = [] # Input patches to classify
  outputImg = np.zeros((img.shape[0],img.shape[1])) # Output segmented slab

  for i in range(img.shape[0]):
      for j in range(img.shape[1]):
	#img[i,j]/=255.0
	  if j>w and j+w<img.shape[1] :
	    if i>w and i+w<img.shape[0]:
	      index.append([i,j]) #Save index
	      inputImg.append(img[i-w:i+w+1,j-w:j+w+1]) #Save input patch

  print "Number of pixels to classify : ", len(inputImg)

  status = 0# Status report variable
  # Classify pixels
  for i in range(100*512,len(inputImg)-100*512,1) : #WARNING: Reduce images size for speed
    net.blobs['data'].data[0,0] = inputImg[i]
    ## perform classification
    output = net.forward()
    #pixelScore = net.blobs['score'].data[0][1] # the output score vector for the first image in the batch
    pixelClass = output['loss'][0][1] # the output probability vector for the first image in the batch
    ## set output
    outputImg[index[i][0],index[i][1]]=pixelClass
    ## Print status
    #print str(i)+" : ",net.blobs['score'].data[0].argmax(), pixelClass, net.blobs['score'].data[0][0],net.blobs['score'].data[0][1]
    if i % (len(inputImg)//20) == 0:
      print ".." + str(status) + "%"
      status = status + 5

  # Save output
  plt.imsave( outputFilename, outputImg, cmap='Greys_r' )


########
# Main #
########
# Model weights
caffeModelWeights = caffe_root + 'data/SegmentVesselsUsingNeuralNetworks/NetProto/net_iter_10000.caffemodel'

# Sanity check
if os.path.isfile( caffeModelWeights ):
  print("Caffe model weights " + caffeModelWeights + " found." )
else :
  print('Error : Caffe model weights "' + caffeModelWeights + '" not found.\n'+
    '  Train neural network before processing it.')
  sys.exit(0)

# Model definition
caffeModelDef = caffe_root + 'data/SegmentVesselsUsingNeuralNetworks/NetProto/deploy.prototxt'

# Sanity check
if os.path.isfile( caffeModelDef ):
  print("Caffe model definition " + caffeModelDef + " found." )
else :
  print('Error : Caffe model definition "' + caffeModelDef + '" not found.\n'+
    '  Train neural network before processing it.')
  sys.exit(0)

# Create network
caffe.set_mode_gpu()

net = caffe.Net(caffeModelDef,      # defines the structure of the model
                caffeModelWeights,  # contains the trained weights
                caffe.TEST)     # use test mode (e.g., don't perform dropout)

# set the size of the input (we can skip this if we're happy
#  with the default; we can also change it later, e.g., for different batch sizes)
w = 32;
net.blobs['data'].reshape(1,         # batch size
                          1,         # 1-channel
                          2*w+1, 2*w+1)  # image size

# Create and segment inputs
hardDrive_root = "/media/lucas/krs0014/"

inputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/images/"
outputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/output/"
inputAnimal = "pp07_A34_left.png" #WARNING Hardcoded

## Segment all animal slabs
for i in range(10):
  inputFilename = inputDir + str(i) + "_" + inputAnimal
  outputFilename = outputDir + str(i) + "_" + inputAnimal
  print(inputFilename)
  segmentImage( net, inputFilename, outputFilename )
