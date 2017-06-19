#!/usr/bin/env python

import os.path

import ctk_cli
import keras.models as M

import deploy
import utils
from utils import script_params


def main(args):
    utils.set_params_path(args.params)
    if (args.resampled is None) ^ (script_params['RESAMPLE_SPACING'] is None or args.preprocessed is None):
        raise ValueError("A resampled image should be supplied iff resampling is"
                         " enabled in the parameters file and a preprocessed"
                         " image is given.")
    if args.preprocessed is None:
        args.resampled, args.preprocessed = deploy.prep(args.inputImage, args.outputDir)
    model = M.load_model(args.model)
    prefix = os.path.join(args.outputDir, os.path.splitext(os.path.basename(args.inputImage))[0])
    deploy.generate_seed_points(model, args.preprocessed, prefix + '_vess_prob.mha')
    deploy.segmentTubes(args.resampled, args.vascularModelFile, prefix,
                        script_params['VESSEL_SEED_PROBABILITY'],
                        script_params['VESSEL_SCALE'])


if __name__ == '__main__':
    main(ctk_cli.CLIArgumentParser().parse_args())
