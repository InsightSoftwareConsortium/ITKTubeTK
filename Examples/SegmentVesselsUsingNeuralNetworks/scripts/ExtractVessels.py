#!/usr/bin/python

###########################################################################
# ExtractVessels.py:
#
# Extract vessels using SegmentTubes CLI. The network reconstructed output
# is used as extraction seeds. The .mtp vascular model has been computed
# from a previous extraction using EnhancetubeUsingDiscriminantAnalysis CLI.
# The scale image can also be obtained that way, or can be replaced by a
# constant.
#
###########################################################################

import sys
import os
import subprocess

def computeSeeds( inputImage, outputImage ):
  subprocess.call( ["ImageMath", inputImage,
                    "-t", "0", "8", "0", "1",
                    "-W","0", outputImage] )

def segmentTubes( inputImage, outputTRE, seedImage, scaleImage, vascularModel ):
  subprocess.call( ["SegmentTubes",
                    #"-o", tmpDir + "/out_VesselMask.mha",
                    "-P", vascularModel,
                    "-M", seedImage,
                    "-s", "0.1",       #Scale constant
                    #"-S", scaleImage, # or scale image
                    inputImage,
                    outputTRE ] )
  # Fill gaps
  subprocess.call( ["ConvertTubesToTubeTree",
                     "--maxTubeDistanceToRadiusRatio", "3",
                     "--removeOrphanTubes",
                    outputTRE,
                    outputTRE + "Tree.tre"] )
  subprocess.call( ["TreeMath",
                    "-f", "S",
                    outputTRE + "Tree.tre",
                    "-w", outputTRE + "Tree.tre"] )

########
# Main #
########
# Path variable
caffe_root = "./"
hardDrive_root = "/media/lucas/krs0014/"

extractionDir = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/Extraction/"

# Images path
seedImage = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/output/3d.mha"
inputImage = caffe_root + "data/SegmentVesselsUsingNeuralNetworks/Controls/A34/pp07_A34_left.mhd"
outputTRE = extractionDir + "out.tre"
scaleImage = extractionDir + "input_ExpertEnhancedScales.mha"
vascularModel = extractionDir + "vascularModel.mtp"

# Process extraction
computeSeeds( seedImage, extractionDir + "seed.mha" )
segmentTubes( inputImage, outputTRE, extractionDir+"seed.mha", scaleImage, vascularModel )