#!/usr/bin/env python

from glob import glob
import os
import random

import itk
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np

import utils
from utils import script_params

# Define paths
output_data_root = script_params['OUTPUT_DATA_ROOT']
input_data_root = script_params['INPUT_DATA_ROOT']

test_data_dir = os.path.join(output_data_root, "testing")
test_output_dir = os.path.join(test_data_dir, 'cnn')
stats_base = os.path.join(test_data_dir, 'stats')

def whole_image_confusion():
    """Compute true/false positive/negative rates, normalizing so that the
    underlying negatives and positives are in equal proportion.

    """
    print("Generating whole-image confusion matrices")
    base = os.path.join(stats_base, 'whole_image_confusion')
    utils.ensureDirectoryExists(base)
    # TODO figure out the output shape of SegmentTubes (the shape of
    # the _vseg image).  It will at any rate be a different size than
    # the expert image; an interesting question is whether it will be
    # the same as the original image.
    name_keys = [os.path.basename(x)[:-9] for x in glob(os.path.join(test_output_dir, '*_vseg.mha'))]
    for name in name_keys:
        print(name)
        # TODO this line will be a problem; cf. previous TODO
        expert_im = itk.imread(str(os.path.join(test_data_dir, name + '_prepped_expert.mha')))
        network_im = itk.imread(str(os.path.join(test_output_dir, name + '_vseg.mha')))
        expert_arr, network_arr = map(itk.GetArrayViewFromImage, (expert_im, network_im))
        network_arr = np.where(network_arr, 1, 0)
        bins = np.bincount((2 * expert_arr + network_arr).reshape(-1), minlength=4).astype(float)
        bins[:2] /= bins[:2].sum()
        bins[2:] /= bins[2:].sum()
        bins /= 2
        with open(os.path.join(base, name + '.txt'), 'w') as f:
            for b, n in zip(bins, [
                    'True negatives',
                    'False positives',
                    'False negatives',
                    'True positives',
            ]):
                f.write('{}: {}\n'.format(n, b))

def whole_image_roc():
    """Plot an ROC curve for the whole-image probability image."""
    print("Generating whole-image ROC curves")
    base = os.path.join(stats_base, 'whole_image_roc')
    utils.ensureDirectoryExists(base)
    name_keys = [os.path.basename(x)[:-14] for x in glob(os.path.join(test_output_dir, '*_vess_prob.mha'))]
    for name in name_keys:
        print(name)
        expert_im = itk.imread(str(os.path.join(test_data_dir, name + '_prepped_expert.mha')))
        network_im = itk.imread(str(os.path.join(test_output_dir, name + '_vess_prob.mha')))
        expert_arr, network_arr = map(itk.GetArrayViewFromImage, (expert_im, network_im))
        all_bins = np.bincount((256 * expert_arr + network_arr).reshape(-1), minlength=512).astype(float)
        write_roc_plot_from_all_bins(all_bins, os.path.join(base, name + '.png'))

def sampling_roc(samples_per_class=5000):
    """Sample so many positive and negative patches per image, and use the
    model prediction on these to estimate the ROC curve.  Also
    generate the confusion matrix using VESSEL_SEED_PROBABILITY as the
    threshold for reference.

    """
    print("Generating sampled ROC curves")
    base = os.path.join(stats_base, 'sampling_roc')
    utils.ensureDirectoryExists(base)
    w = script_params['PATCH_RADIUS']
    model = utils.load_best_model()
    name_keys = [os.path.basename(x)[:-12] for x in glob(os.path.join(test_data_dir, '*_prepped.mha'))]
    for name in name_keys:
        print(name)
        expert_im = itk.imread(str(os.path.join(test_data_dir, name + '_prepped_expert.mha')))
        prepped_im = itk.imread(str(os.path.join(test_data_dir, name + '_prepped.mha')))
        expert_arr, prepped_arr = map(itk.GetArrayViewFromImage, (expert_im, prepped_im))
        sample_ind = [random.sample(np.transpose(np.where(
            (expert_arr if i else 1 - expert_arr)[(np.s_[w:-w],) * expert_arr.ndim]
        )) + w, samples_per_class) for i in 0, 1]
        neg_pred, pos_pred = (utils.predict_on_indices(
            model, prepped_arr, sample_ind[0] + sample_ind[1],
            script_params['DEPLOY_BATCH_SIZE'],
        ).reshape((2, -1)) * 255).round().astype(int)
        all_bins = np.concatenate((np.bincount(neg_pred, minlength=256),
                                   np.bincount(pos_pred, minlength=256))).astype(float)
        write_roc_plot_from_all_bins(all_bins, os.path.join(base, name + '.png'))
        # TODO ROC curve

def write_roc_plot_from_all_bins(all_bins, path):
    """From a flat array of all expert-negative, the all expert-positive
    bincounts, create an ROC plot

    """
    pairs = []
    for i in range(257):
        bins = np.array([all_bins[:i].sum(), all_bins[i:256].sum(),
                         all_bins[256:i + 256].sum(), all_bins[i + 256:].sum()])
        pairs.append((bins[1]/bins[:2].sum(), bins[3]/bins[2:].sum()))
    pairs = np.stack(pairs, -1)
    plt.figure()
    plt.plot(pairs[0], pairs[1], 'k.-')
    plt.savefig(path)

def main():
    whole_image_confusion()
    whole_image_roc()
    sampling_roc()

if __name__ == '__main__':
    main()
