"""Routines used for applying an already trained network to images"""

import os.path
import subprocess
import time

import itk
import numpy as np

import utils
from utils import script_params


def duplicate(im):
    """Duplicate an itk image"""
    f = itk.ImageDuplicator.New(im)
    f.Update()
    return f.GetOutput()


# Preprocess ("prep") images
def prep(inputImage, outputDir, expertImage=None):
    """Preprocess inputImage and expertImage (if not None) according to
    script_params.  Output (where '*' stands for outputDir +
    basename(inputImage) (without extension)):
    - *_prepped.mha: Preprocessed inputImage
    - *_prepped_expert.mha: Preprocessed expertImage

    """
    outputImagePrefix = os.path.join(outputDir, os.path.splitext(os.path.basename(inputImage))[0])
    outputImagePrefix = str(outputImagePrefix)

    smoothing_radius = script_params['SMOOTHING_RADIUS']

    reader = itk.ImageFileReader.New(FileName=str(inputImage))
    smoothing_filter = itk.MedianImageFilter.New(reader.GetOutput(),
                                                 Radius=smoothing_radius)
    equalization_filter = itk.AdaptiveHistogramEqualizationImageFilter.New(
        smoothing_filter.GetOutput(),
        Radius=script_params['PATCH_RADIUS'],
        Alpha=0, Beta=0,
    )
    writer = itk.ImageFileWriter.New(equalization_filter.GetOutput(),
                                     FileName=outputImagePrefix + "_prepped.mha",
                                     UseCompression=True)
    writer.Update()

    if expertImage is None:
        return writer.GetFileName()
    else:
        utils.symlink_through(expertImage, outputImagePrefix + '_prepped_expert.mha')


def segmentPreppedImage(model, input_file, output_file):
    """Segment (really, generate seed points from) a preprocessed image"""

    print "Segmenting image", input_file

    data_shape = model.input_shape

    print data_shape

    # read input slab image
    input_image_itk = itk.imread(str(input_file))
    input_image = itk.GetArrayViewFromImage(input_image_itk)

    # get foreground mask
    input_revcum = np.cumsum(np.bincount(input_image.reshape(-1))[::-1])[::-1]
    th = np.count_nonzero(input_revcum >= input_revcum[0] * script_params['DEPLOY_TOP_FRAC']) - 2
    fgnd_mask = input_image > th

    # get test_batch_size and patch_size used for cnn model
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

    prob_vessel = utils.predict_on_indices(model, input_image, patch_indices, test_batch_size)
    output_image[tuple(patch_indices.T)] = (prob_vessel * 255).round()

    end_time = time.time()
    print '\tTook %s seconds' % (end_time - start_time)

    # Save output
    itk.imwrite(output_image_itk, str(output_file), compression=True)


def segmentTubes(originalImage, vascularModelFile, outputDir,
                 vess_seed_prob=0.95, vess_scale=0.1):
    inputImageName = os.path.splitext(os.path.basename(originalImage))[0]

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
                     originalImage, outVsegTreFile])

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
