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

import json
import os
import subprocess


def computeSeeds( inputImage, outputImage ):
    subprocess.call(["ImageMath", inputImage,
                     "-t", "0", "8", "0", "1",
                     "-W", "0", outputImage])


def segmentTubes(inputImage, outputTRE, seedImage, scaleImage, vascularModel):

    subprocess.call(["SegmentTubes",
                    #"-o", tmpDir + "/out_VesselMask.mha",
                    "-P", vascularModel,
                    "-M", seedImage,
                    "-s", "0.1",       #Scale constant
                    #"-S", scaleImage, # or scale image
                    inputImage,
                    outputTRE])
    # Fill gaps
    subprocess.call(["ConvertTubesToTubeTree",
                     "--maxTubeDistanceToRadiusRatio", "3",
                     "--removeOrphanTubes",
                    outputTRE,
                    outputTRE + "Tree.tre"])
    subprocess.call(["TreeMath",
                    "-f", "S",
                    outputTRE + "Tree.tre",
                    "-w", outputTRE + "Tree.tre"])

########
# Main #
########

# Set path variables
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']

proj_name = "SegmentVesselsUsingNeuralNetworks"
caffe_vseg_root = os.path.join(caffe_root, "data", proj_name)
hardDrive_vseg_root = os.path.join(hardDrive_root, proj_name)

extractionDir = os.path.join(caffe_vseg_root, "Extraction")

# Images path
seedImage = os.path.join(hardDrive_vseg_root, "output/3d.mha")
inputImage = os.path.join(caffe_vseg_root, "Controls/A36/pp07_A36_left.mhd")
outputTRE = os.path.join(extractionDir, "out.tre")
scaleImage = os.path.join(extractionDir, "input_ExpertEnhancedScales.mha")
vascularModel = os.path.join(extractionDir, "vascularModel.mtp")

# Process extraction
computeSeeds(seedImage, os.path.join(extractionDir, "seed.mha"))
segmentTubes(inputImage, outputTRE, os.path.join(extractionDir, "seed.mha"), scaleImage, vascularModel)
