#!/usr/bin/env python

import errno
import os

import ctk_cli
import keras.models as M

from tubetk.vseg.cnn import deploy, utils

script_params = utils.script_params

def main(args):
    utils.set_params_path(args.params)
    if (args.resampled is None) ^ (script_params['RESAMPLE_SPACING'] is None or args.preprocessed is None):
        raise ValueError("A resampled image should be supplied iff resampling is"
                         " enabled in the parameters file and a preprocessed"
                         " image is given.")
    try:
        os.mkdir(args.outputDir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    if args.preprocessed is None:
        args.resampled, args.preprocessed = deploy.prep(args.inputImage, args.outputDir)
    elif args.resampled is None:
        args.resampled = args.inputImage
    model = M.load_model(args.model)
    prefix = os.path.join(args.outputDir, os.path.splitext(os.path.basename(args.inputImage))[0])
    deploy.generate_seed_points(model, args.preprocessed, prefix)
    deploy.segmentTubes(args.resampled, args.vascularModelFile, prefix,
                        script_params['VESSEL_SEED_PROBABILITY'],
                        script_params['VESSEL_SCALE'])


if __name__ == '__main__':
    main(ctk_cli.CLIArgumentParser().parse_args())
