#!/usr/bin/env python

from glob import glob
import os

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
    name_keys = [os.path.basename(x)[:-9] for x in glob(os.path.join(test_output_dir, '*_vseg.mha'))]
    for name in name_keys:
        print(name)
        expert_im = itk.imread(str(os.path.join(test_data_dir, name + '_expert.mha')))
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
        expert_im = itk.imread(str(os.path.join(test_data_dir, name + '_expert.mha')))
        network_im = itk.imread(str(os.path.join(test_output_dir, name + '_vess_prob.mha')))
        expert_arr, network_arr = map(itk.GetArrayViewFromImage, (expert_im, network_im))
        pairs = []
        for i in range(257):
            thresh = network_arr >= i
            bins = np.bincount((2 * expert_arr + thresh).reshape(-1), minlength=4).astype(float)
            pairs.append((bins[1]/bins[:2].sum(), bins[3]/bins[2:].sum()))
        pairs = np.stack(pairs, -1)
        plt.figure()
        plt.plot(pairs[0], pairs[1], 'k.-')
        plt.savefig(os.path.join(base, name + '.png'))

def main():
    whole_image_confusion()
    whole_image_roc()

if __name__ == '__main__':
    main()
