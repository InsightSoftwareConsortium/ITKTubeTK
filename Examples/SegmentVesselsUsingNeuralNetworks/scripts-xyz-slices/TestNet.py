#!/usr/bin/python
"""
Tests the trained CNN vessel segmentation mode on all test images
"""
import itertools
import os
import sys
import glob
import subprocess
import time

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import skimage.io
import skimage.filters

import utils
from utils import script_params

# Append ITK libs
sys.path.append(os.path.join(
    os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Wrapping/Generators/Python'))
sys.path.append(os.path.join(
    os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Modules/ThirdParty/VNL/src/vxl/lib'))
import itk

# TODO factor out the itk import procedure and the dependency this creates
import deploy

# Define paths
output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

testDataDir = os.path.join(output_data_root, "testing")

import keras.models as M


def run():

    # set output directory
    outputDir = os.path.join(output_data_root, "testing", "cnn")

    utils.ensureDirectoryExists(outputDir)

    sys.stdout = utils.Logger(os.path.join(outputDir, 'net_test.log'))

    # Model definition
    modelDef = os.path.join(
        output_data_root, "NetProto", "net_best.hdf5")

    if os.path.isfile(modelDef):
        print("Model definition " + modelDef + " found.")
    else:
        print('Error : model definition "' + modelDef + '" not found.\n' +
              '  Train neural network before processing it.')
        sys.exit(1)

    # Create network
    model = M.load_model(modelDef)

    # get list of test mha files
    testMhaFiles = glob.glob(os.path.join(testDataDir, "*_prepped.mha"))

    # segment all test .mha files
    for mhaFile in testMhaFiles:

        testAnimal = os.path.basename(os.path.splitext(mhaFile)[0])[:-8]

        # segment image
        deploy.segmentPreppedImage(
            model, mhaFile,
            os.path.join(outputDir, testAnimal + '_vess_prob.mha')
        )

        # segment tubes using ridge traversal
        vascularModelFile = os.path.join(input_data_root, 'vascularModel.mtp')

        originalImage, = glob.glob(os.path.join(
            input_data_root, '*', '*', testAnimal + ".mhd"))

        deploy.segmentTubes(originalImage, vascularModelFile, outputDir,
                            script_params['VESSEL_SEED_PROBABILITY'],
                            script_params['VESSEL_SCALE'])

        break


if __name__ == "__main__":
    run()
