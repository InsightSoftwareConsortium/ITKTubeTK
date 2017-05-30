"""Routines used for applying an already trained network to images"""

import os.path
import subprocess
import time

import itk
import numpy as np

import utils
from utils import script_params


# Preprocess ("prep") images
def prep(inputImage, outputDir):
    """Preprocess inputImage and expertImage (if not None) according to
    script_params.  Output (where '*' stands for outputDir +
    basename(inputImage) (without extension)):
    - *_prepped.mha: Preprocessed inputImage

    """
    outputImagePrefix = os.path.join(outputDir, os.path.splitext(os.path.basename(inputImage))[0])
    outputImagePrefix = str(outputImagePrefix)

    smoothing_radius = script_params['SMOOTHING_RADIUS']

    reader = itk.ImageFileReader.New(FileName=str(inputImage))
    if script_params['RESAMPLE_SPACING'] is not None:
        reader.UpdateOutputInformation()
        size, spacing = itk.size(reader), itk.spacing(reader)
        new_spacing = np.full(len(spacing), script_params['RESAMPLE_SPACING'])
        spacing_ratio = new_spacing / spacing
        # Preserve the physical location of index -0.5 (corner of pixel 0)
        new_origin = reader.GetOutput().TransformContinuousIndexToPhysicalPoint(
            itk.ContinuousIndex[itk.D, reader.GetOutput().GetImageDimension()](0.5 * (spacing_ratio - 1))
        )
        new_size = (size / spacing_ratio).round().astype(int)
        Interpolator = itk.BSplineInterpolateImageFunction[
            reader.GetOutput(),
            itk.D,
            itk.template(reader.GetOutput())[1][0],
        ]
        resample = itk.ResampleImageFilter.New(
            reader,
            Size=list(new_size),
            OutputSpacing=new_spacing,
            OutputOrigin=new_origin,
            OutputDirection=reader.GetOutput().GetDirection(),
            Interpolator=Interpolator.New(UseImageDirection=True),
        )
        to_smoothing = resample
    else:
        to_smoothing = reader

    smoothing_filter = itk.MedianImageFilter.New(to_smoothing,
                                                 Radius=smoothing_radius)

    CIF = itk.CastImageFilter[type(smoothing_filter.GetOutput()),
                              itk.Image[itk.F, reader.GetOutput().GetImageDimension()]]
    smoothed_single = CIF.New(smoothing_filter)

    writer = itk.ImageFileWriter.New(smoothed_single,
                                     FileName=outputImagePrefix + "_prepped.mha",
                                     UseCompression=True)
    writer.Update()

    return writer.GetFileName()


def chunked_argmax(arr, window):
    """Compute a tuple (length arr.ndim) of arrays representing the
    indices of the max values in window-shaped chunks of arr, which
    must evenly divide into windows.

    """
    split_dim = arr.reshape(tuple(x for s, w in zip(arr.shape, window) for x in (s / w, w)))
    transpose_dim = split_dim.transpose(tuple(range(0, split_dim.ndim, 2)) +
                                        tuple(range(1, split_dim.ndim, 2)))
    flat_dim = transpose_dim.reshape(transpose_dim.shape[:arr.ndim] + (-1,))
    argmaxes = np.argmax(flat_dim, axis=-1)
    flat_argmaxes = argmaxes.reshape(-1)
    in_chunk_indices = np.unravel_index(flat_argmaxes, window)
    chunk_corner_indices = (np.indices(argmaxes.shape).reshape((arr.ndim, -1)).T * window).T
    return tuple(chunk_corner_indices + in_chunk_indices)


def segmentPreppedImage(model, input_file, output_file):
    """Segment (really, generate seed points from) a preprocessed image"""

    print "Segmenting image", input_file

    data_shape = model.input_shape

    print data_shape

    # read input slab image
    input_image_itk = itk.imread(str(input_file))
    input_image = itk.GetArrayViewFromImage(input_image_itk)

    # get foreground mask
    tw = script_params['DEPLOY_TOP_WINDOW']
    input_image_padded = np.pad(input_image, tuple(
        (0, -s % tw) for s in input_image.shape
    ), 'constant', constant_values=-np.inf)
    fgnd_mask = np.zeros_like(input_image, dtype=bool)
    fgnd_mask[chunked_argmax(input_image_padded, (tw,) * input_image.ndim)] = True

    # get test_batch_size and patch_size used for cnn model
    test_batch_size = script_params['DEPLOY_BATCH_SIZE']

    design = script_params['NETWORK_DESIGN']
    if design == 'xyz':
        patch_size = data_shape[0][1]
    elif design == 'full3d':
        patch_size = data_shape[1]
    else:
        raise ValueError("Unknown NETWORK_DESIGN")

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

    output_image = np.zeros_like(input_image, dtype=np.uint8)

    prob_vessel = utils.predict_on_indices(model, input_image, patch_indices, test_batch_size)
    output_image[tuple(patch_indices.T)] = (prob_vessel * 255).round()

    end_time = time.time()
    print '\tTook %s seconds' % (end_time - start_time)

    # Save output
    output_image_itk = itk.GetImageViewFromArray(output_image)
    output_image_itk.CopyInformation(input_image_itk)
    itk.imwrite(output_image_itk, str(output_file), compression=True)


def segmentTubes(prepped_image, vascularModelFile, output_prefix,
                 vess_seed_prob=0.95, vess_scale=0.1):
    # compute seed image
    vessProbImageFile = output_prefix + "_vess_prob.mha"
    outSeedImageFile = output_prefix + "_vess_seeds.mha"

    subprocess.call(["ImageMath", vessProbImageFile,
                     "-t", str(255 * vess_seed_prob), "255", "1", "0",
                     "-W", "0", outSeedImageFile])

    # segment tubes using ridge traversal
    outVsegMaskFile = output_prefix + "_vseg.mha"
    outVsegTreFile = output_prefix + "_vseg.tre"

    subprocess.call(["SegmentTubes",
                     "-o", outVsegMaskFile,
                     "-P", vascularModelFile,
                     "-M", outSeedImageFile,
                     "-s", str(vess_scale),
                     prepped_image, outVsegTreFile])

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
