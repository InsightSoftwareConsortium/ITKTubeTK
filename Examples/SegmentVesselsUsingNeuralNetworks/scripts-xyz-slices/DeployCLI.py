#!/usr/bin/env python

import os.path

import ctk_cli
import keras.models as M

import deploy
import utils
from utils import script_params


def main(args):
    utils.set_params_path(args.params)
    if args.preprocessed is None:
        args.preprocessed = deploy.prep(args.inputImage, args.outputDir)
    model = M.load_model(args.model)
    prefix = os.path.join(args.outputDir, os.path.splitext(os.path.basename(args.inputImage))[0])
    deploy.segmentPreppedImage(model, args.preprocessed, prefix + '_vess_prob.mha')
    deploy.segmentTubes(args.preprocessed, args.vascularModelFile, prefix,
                        script_params['VESSEL_SEED_PROBABILITY'],
                        script_params['VESSEL_SCALE'])


if __name__ == '__main__':
    main(ctk_cli.CLIArgumentParser().parse_args())
