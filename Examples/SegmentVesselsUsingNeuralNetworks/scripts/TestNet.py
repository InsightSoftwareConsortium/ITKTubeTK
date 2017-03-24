#!/usr/bin/python
"""
Tests the trained CNN vessel segmentation mode on all test images
"""
import os
import sys
import json
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

# Append ITK libs
sys.path.append(os.path.join(
    os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Wrapping/Generators/Python'))
sys.path.append(os.path.join(
    os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Modules/ThirdParty/VNL/src/vxl/lib'))
import itk

# Define paths
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

testDataDir = os.path.join(output_data_root, "testing")

# import caffe
sys.path.insert(0, os.path.join(caffe_root, 'python'))
import caffe


# Segment a slab
def segmentSlab(net, input_file, output_file):

    print "Segmenting slab ", input_file

    print net.blobs['data'].data.shape

    # read input slab image
    input_image = skimage.io.imread(input_file)

    # get foreground mask
    th = skimage.filters.threshold_otsu(input_image)
    fgnd_mask = input_image > th

    # get test_batch_size and patch_size used for cnn net
    test_batch_size = net.blobs['data'].data.shape[0]
    patch_size = net.blobs['data'].data.shape[2]

    print 'Test batch shape = ', net.blobs['data'].data.shape

    # collect all patches
    print "Extracting patches ... "

    start_time = time.time()

    w = np.int(patch_size / 2)

    patches = []
    patch_indices = []

    for i in range(w, input_image.shape[0] - w):
        for j in range(w, input_image.shape[1] - w):

            # check if pixel is in foreground
            # if ~fgnd_mask[i, j]:
            #    continue

            # store patch center index
            patch_indices.append([i, j])

            # store patch
            cur_patch = input_image[(i - w):(i + w + 1), (j - w):(j + w + 1)]
            patches.append(cur_patch)

    end_time = time.time()

    print '\tTook %s seconds' % (end_time - start_time)

    num_patches = len(patch_indices)
    patches = np.array(patches)
    patch_indices = np.array(patch_indices)

    print "\tNo of patches = %s" % num_patches

    # Classify patches using cnn and write result in output image
    print "Classifying patches ... "

    start_time = time.time()

    output_image = np.zeros_like(input_image)  # Output segmented slab
    output_image_flat = np.ravel(output_image)

    pw = -num_patches % test_batch_size
    patches = np.pad(patches, ((0, pw), (0, 0), (0, 0)), 'constant')
    patch_indices = np.pad(patch_indices, ((0, pw), (0, 0)), 'constant')

    patches = patches[:, np.newaxis]
    print patches.shape

    for i in range(0, num_patches, test_batch_size):

        # get current batch of patches
        cur_pind = patch_indices[i:i + test_batch_size]
        cur_patches = patches[i:i + test_batch_size]

        # perform classification using cnn
        net.blobs['data'].data[...] = cur_patches

        output = net.forward()

        prob_vessel = output['loss'][:, 1]

        if i + test_batch_size > num_patches:

            # remove padded part for last batch
            num_valid = num_patches - i
            prob_vessel = prob_vessel[0:num_valid]
            cur_pind = cur_pind[0:num_valid, :]

        # set output
        cur_pind_lin = np.ravel_multi_index((cur_pind[:, 0], cur_pind[:, 1]),
                                            output_image.shape)

        output_image_flat[cur_pind_lin] = prob_vessel * 255

        # Print progress
        print '\t %.2f%%' % (100.0 * (i + len(prob_vessel)) / num_patches),
        if len(prob_vessel) > 0:
            print "%.4f, %.4f, %.4f" % (prob_vessel.min(),
                                        prob_vessel.max(),
                                        prob_vessel.mean())
        else:
            print ""

    end_time = time.time()
    print '\tTook %s seconds' % (end_time - start_time)

    # Save output
    skimage.io.imsave(output_file, output_image)


def combineSlabs(inputImageName, outputDir):

    PixelType = itk.UC
    Dimension = 3
    ImageType = itk.Image[PixelType, Dimension]

    NameGeneratorType = itk.NumericSeriesFileNames
    nameGenerator = NameGeneratorType.New()
    nameGenerator.SetStartIndex(0)
    nameGenerator.SetEndIndex(9)
    nameGenerator.SetIncrementIndex(1)
    nameGenerator.SetSeriesFormat(
        str(os.path.join(outputDir, "%d_" + inputImageName + '_zslab.png')))

    SeriesReaderType = itk.ImageSeriesReader[ImageType]
    seriesReader = SeriesReaderType.New()
    seriesReader.SetFileNames(nameGenerator.GetFileNames())
    #  seriesReader.SetImageIO( itk.PNGImageIO.New() )

    WriterType = itk.ImageFileWriter[ImageType]
    writer = WriterType.New()
    writer.SetFileName(
        str(os.path.join(outputDir, inputImageName + "_zslab.mha")))
    writer.SetInput(seriesReader.GetOutput())
    writer.Update()


def reconstructSegVolume(inputImageName, outputDir):

    PixelType = itk.F
    Dimension = 3
    ImageType = itk.Image[PixelType, Dimension]

    # combine slab images into a .mha volume file
    combineSlabs(inputImageName, outputDir)

    # Read segmented slab volume
    slabSegFile = str(os.path.join(outputDir, inputImageName + "_zslab.mha"))

    SlabVolumeReaderType = itk.ImageFileReader[ImageType]
    slabVolumeReader = SlabVolumeReaderType.New()
    slabVolumeReader.SetFileName(slabSegFile)
    slabVolumeReader.Update()
    slabSeg = slabVolumeReader.GetOutput()
    slabSegBuf = itk.PyBuffer[ImageType].GetArrayFromImage(slabSeg)

    # Read point map image
    pointMapFile = str(os.path.join(
        testDataDir, "points", inputImageName + "_zslab_points.mha"))

    VectorImageType = itk.VectorImage[PixelType, Dimension]
    PointMapReaderType = itk.ImageFileReader[VectorImageType]
    pointMapReader = PointMapReaderType.New()
    pointMapReader.SetFileName(pointMapFile)
    pointMapReader.Update()
    pointMap = pointMapReader.GetOutput()
    pointMapBuf = itk.PyBuffer[VectorImageType].GetArrayFromImage(pointMap)

    # Read input volume
    inputVolumeFile = str(os.path.join(os.path.join(
        testDataDir, "images", inputImageName + ".mha")))

    InputVolumeReaderType = itk.ImageFileReader[ImageType]
    inputVolumeReader = InputVolumeReaderType.New()
    inputVolumeReader.SetFileName(inputVolumeFile)
    inputVolumeReader.Update()

    # reconstruct segmentation volume in original space by inverse mapping
    img3d = inputVolumeReader.GetOutput()

    img3d.FillBuffer(0.0)

    numSlabs = slabSegBuf.shape[0]

    for slabIdx in range(numSlabs):

        print "Reconstructing slab ", slabIdx + 1, "/", numSlabs, "..."

        for i in range(slabSegBuf.shape[1]):
            for j in range(slabSegBuf.shape[2]):

                predProb = slabSegBuf[slabIdx, i, j] / 255.0

                src_x = float(pointMapBuf[slabIdx, i, j, 0])
                src_y = float(pointMapBuf[slabIdx, i, j, 1])
                src_z = float(pointMapBuf[slabIdx, i, j, 2])
                src_index = img3d.TransformPhysicalPointToIndex(
                    [src_x, src_y, src_z])

                img3d.SetPixel(src_index, predProb)

    # Blur image to fill gaps between painted pixels.
    BlurFilterType = itk.DiscreteGaussianImageFilter[ImageType, ImageType]
    blurFilter = BlurFilterType.New()
    blurFilter.SetInput(img3d)
    blurFilter.SetVariance(1)
    blurFilter.SetMaximumKernelWidth(1)

    print("Writing output ..")
    WriterType = itk.ImageFileWriter[ImageType]
    writer = WriterType.New()
    writer.UseCompressionOn()
    writer.SetFileName(
        str(os.path.join(outputDir, inputImageName + "_vess_prob.mha")))
    writer.SetInput(blurFilter.GetOutput())
    writer.Update()


def segmentTubes(inputImageName, vascularModelFile, outputDir,
                 vess_seed_prob=0.95, vess_scale=0.1):

    inputImageFile = os.path.join(
        testDataDir, "images", inputImageName + ".mha")

    # compute seed image
    vessProbImageFile = os.path.join(
        outputDir, inputImageName + "_vess_prob.mha")
    outSeedImageFile = os.path.join(
        outputDir, inputImageName + "_vess_seeds.mha")

    subprocess.call(["ImageMath", vessProbImageFile,
                     "-t", str(vess_seed_prob), "1.0", "1", "0",
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

    # Model weights
    caffeModelWeights = os.path.join(
        output_data_root, "NetProto", "net_best.caffemodel")

    if os.path.isfile(caffeModelWeights):
        print("Caffe model weights " + caffeModelWeights + " found.")
    else:
        print('Error : Caffe model weights "' + caffeModelWeights + '" not found.\n' +
              '  Train neural network before processing it.')
        sys.exit(0)

    # Model definition
    caffeModelDef = os.path.join(
        output_data_root, "NetProto", "net_deploy.prototxt")

    if os.path.isfile(caffeModelDef):
        print("Caffe model definition " + caffeModelDef + " found.")
    else:
        print('Error : Caffe model definition "' + caffeModelDef + '" not found.\n' +
              '  Train neural network before processing it.')
        sys.exit(0)

    # Create network
    caffe.set_mode_gpu()

    net = caffe.Net(str(caffeModelDef),  # defines the structure of the model
                    str(caffeModelWeights),  # contains the trained weights
                    caffe.TEST)  # use test mode (e.g., don't perform dropout)

    # set the size of the input (we can skip this if we're happy
    # with the default; we can also change it later, e.g., for different batch
    # sizes)
    w = script_params['PATCH_RADIUS']

    net.blobs['data'].reshape(script_params['DEPLOY_BATCH_SIZE'],  # batch size
                              1,  # 1-channel
                              2 * w + 1, 2 * w + 1)  # image size

    # get list of test mha files
    testMhaFiles = glob.glob(os.path.join(
        testDataDir, "images", "*_zslab.mha"))

    # segment all test .mha files
    for mhaFile in testMhaFiles:

        testAnimal = os.path.basename(os.path.splitext(mhaFile)[0])[:-6]
        print "Segmenting ", testAnimal + ".mha", " ..."

        # segment all slabs
        for i in range(script_params['NUM_SLABS']):

            inputFilename = os.path.join(
                testDataDir, "images", str(i) + "_" + testAnimal + '_zslab.png')

            outputFilename = os.path.join(
                outputDir, str(i) + "_" + testAnimal + '_zslab.png')

            segmentSlab(net, inputFilename, outputFilename)

        # reconstruct volume from segmented slabs by mapping back to their MIP
        # location
        reconstructSegVolume(testAnimal, outputDir)

        # segment tubes using ridge traversal
        vascularModelFile = os.path.join(input_data_root, 'vascularModel.mtp')

        segmentTubes(testAnimal, vascularModelFile, outputDir,
                     script_params['VESSEL_SEED_PROBABILITY'],
                     script_params['VESSEL_SCALE'])

        break

if __name__ == "__main__":
    run()
