#!/usr/bin/env python

"""Tune the free parameters, debug, and characterize analysis performed by
RegisterImageToTubesUsingRigidTransform."""

import argparse
import json
import subprocess
import sys


def run_analysis(config, config_file):
    io_params = config['ParameterGroups'][0]['Parameters']
    print(io_params[2]['Value'])
    subprocess.check_call([config['Executables']['Analysis'],
                          '--parameterstorestore', config_file.name,
                          io_params[0]['Value'],
                          io_params[1]['Value'],
                          io_params[2]['Value']])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=argparse.FileType('r'),
                        help='Configuration file for tuning analysis')
    args = parser.parse_args()
    config_file = args.configuration

    config = json.load(config_file)
    run_analysis(config, config_file)
