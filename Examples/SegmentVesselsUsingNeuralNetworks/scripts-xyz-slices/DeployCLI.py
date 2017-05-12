#!/usr/bin/env python

# TODO SEM CLI-ize

import argparse
import os.path

import keras.models as M

import deploy
import utils
from utils import script_params

parser = argparse.ArgumentParser(description='Segment vessels using a neural network')
parser.add_argument('params', help='Path to the params.json file used to train the model')
parser.add_argument('model', help='Path to a saved Keras model file')
parser.add_argument('vascularModelFile', help='Path to file with vascular model parameters')
parser.add_argument('inputImage')
parser.add_argument('outputDir', nargs='?', default='.', help='Directory to write all output to.  Output file names all start with the basename of inputImage without its extension')
parser.add_argument('--preprocessed', help='Path to an appropriately preprocessed version of inputImage.  Automatically generated if not given.')


def main(args):
    utils.set_params_path(args.params)
    if args.preprocessed is None:
        args.preprocessed = deploy.prep(args.inputImage, args.outputDir)
    model = M.load_model(args.model)
    prefix = os.path.join(args.outputDir, os.path.splitext(os.path.basename(args.inputImage))[0])
    deploy.segmentPreppedImage(model, args.preprocessed, prefix + '_vess_prob.mha')
    deploy.segmentTubes(args.inputImage, args.vascularModelFile, args.outputDir,
                        script_params['VESSEL_SEED_PROBABILITY'],
                        script_params['VESSEL_SCALE'])


if __name__ == '__main__':
    main(parser.parse_args())
