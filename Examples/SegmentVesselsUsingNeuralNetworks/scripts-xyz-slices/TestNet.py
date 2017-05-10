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

# Define paths
output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

testDataDir = os.path.join(output_data_root, "testing")

import keras.models as M

def duplicate(im):
    """Duplicate an itk image"""
    f = itk.ImageDuplicator.New(im)
    f.Update()
    return f.GetOutput()

# Segment an image
def segmentImage(net, input_file, output_file):

    print "Segmenting image", input_file

    data_shape = net.input_shape

    print data_shape

    # read input slab image
    input_image_itk = itk.imread(str(input_file))
    input_image = itk.GetArrayViewFromImage(input_image_itk)

    # get foreground mask
    input_revcum = np.cumsum(np.bincount(input_image.reshape(-1))[::-1])[::-1]
    th = np.count_nonzero(input_revcum >= input_revcum[0] * script_params['DEPLOY_TOP_FRAC']) - 2
    fgnd_mask = input_image > th

    # get test_batch_size and patch_size used for cnn net
    test_batch_size = script_params['DEPLOY_BATCH_SIZE']
    patch_size = data_shape[0][1]

    print 'Test batch shape = ', data_shape

    # collect all patches
    print "Extracting patches ... "

    start_time = time.time()

    w = np.int(patch_size / 2)

    patch_indices = np.stack(np.where(fgnd_mask[(np.s_[w:-w],) * input_image.ndim]), axis=-1) + w

    end_time = time.time()

    print '\tTook %s seconds' % (end_time - start_time)

    num_patches = patch_indices.shape[0]

    print "\tNo of patches = %s" % num_patches

    # Classify patches using cnn and write result in output image
    print "Classifying patches ... "

    start_time = time.time()

    output_image_itk = duplicate(input_image_itk)
    output_image = itk.GetArrayViewFromImage(output_image_itk)
    output_image.fill(0)

    for i in range(0, num_patches, test_batch_size):

        # get current batch of patches
        cur_pind = patch_indices[i:i + test_batch_size]
        cur_patches = utils.prepareInputArray(np.stack(
            utils.extractPatch(input_image, indices) for indices in cur_pind
        ))

        # perform classification using cnn
        prob_vessel = net.predict_on_batch(cur_patches)[:, 1]

        output_image[tuple(cur_pind.T)] = (prob_vessel * 255).round()

        # Print progress
        print '\t %.2f%%' % (100.0 * (i + len(prob_vessel)) / num_patches),
        print "%.4f, %.4f, %.4f" % (prob_vessel.min(),
                                    prob_vessel.max(),
                                    prob_vessel.mean())

    end_time = time.time()
    print '\tTook %s seconds' % (end_time - start_time)

    # Save output
    itk.imwrite(output_image_itk, str(output_file), compression=True)


def segmentTubes(inputImageName, vascularModelFile, outputDir,
                 vess_seed_prob=0.95, vess_scale=0.1):

    inputImageFile, = glob.glob(os.path.join(
        input_data_root, '*', '*', inputImageName + ".mhd"))

    # compute seed image
    vessProbImageFile = os.path.join(
        outputDir, inputImageName + "_vess_prob.mha")
    outSeedImageFile = os.path.join(
        outputDir, inputImageName + "_vess_seeds.mha")

    subprocess.call(["ImageMath", vessProbImageFile,
                     "-t", str(255 * vess_seed_prob), "255", "1", "0",
                     "-W", "0", outSeedImageFile])

    # segment tubes using ridge traversal
    outVsegMaskFile = os.path.join(outputDir, inputImageName + "_vseg.mha")
    outVsegTreFile = os.path.join(outputDir, inputImageName + "_vseg.tre")

    subprocess.call(["SegmentTubes",
                     "-o", outVsegMaskFile,
                     "-P", vascularModelFile,
                     "-M", outSeedImageFile,
                     "-s", str(vess_scale),
                     inputImageFile, outVsegTreFile])

    # Fill gaps and convert to a tree
    subprocess.call(["ConvertTubesToTubeTree",
                     "--maxTubeDistanceToRadiusRatio", "3",
                     "--removeOrphanTubes",
                     outVsegTreFile,
                     outVsegTreFile])

    subprocess.call(["TreeMath",
                     "-f", "S",
                     "-w", outVsegTreFile,
                     outVsegTreFile])


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
        sys.exit(0)

    # Create network
    net = M.load_model(modelDef)

    # get list of test mha files
    testMhaFiles = glob.glob(os.path.join(testDataDir, "*_prepped.mha"))

    # segment all test .mha files
    for mhaFile in testMhaFiles:

        testAnimal = os.path.basename(os.path.splitext(mhaFile)[0])[:-8]

        # segment image
        segmentImage(net, mhaFile, os.path.join(outputDir, testAnimal + '_vess_prob.mha'))

        # segment tubes using ridge traversal
        vascularModelFile = os.path.join(input_data_root, 'vascularModel.mtp')

        segmentTubes(testAnimal, vascularModelFile, outputDir,
                     script_params['VESSEL_SEED_PROBABILITY'],
                     script_params['VESSEL_SCALE'])

        break

if __name__ == "__main__":
    run()
