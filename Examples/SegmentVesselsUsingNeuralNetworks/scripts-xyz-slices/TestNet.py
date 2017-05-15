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

import deploy
import utils
from utils import script_params

import itk

# Define paths
output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

testDataDir = os.path.join(output_data_root, "testing")


def run():

    # set output directory
    outputDir = os.path.join(output_data_root, "testing", "cnn")

    utils.ensureDirectoryExists(outputDir)

    sys.stdout = utils.Logger(os.path.join(outputDir, 'net_test.log'))

    try:
        # Model definition
        model = utils.load_best_model()
    except IOError:
        print("Could not load model from expected path!")
        raise

    print("Model definition found.")

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

        originalImage = utils.original_image(testAnimal)

        deploy.segmentTubes(originalImage, vascularModelFile, outputDir,
                            script_params['VESSEL_SEED_PROBABILITY'],
                            script_params['VESSEL_SCALE'])

        break


if __name__ == "__main__":
    run()
