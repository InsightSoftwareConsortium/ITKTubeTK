#!/usr/bin/env python

from __future__ import print_function

import argparse
import json
import os
import shutil
from subprocess import check_call
import sys

from .utils import script_params, symlink_entries_through

stages = [
    'PrepareTrainingData',
    'TrainNet',
    'TestNet',
    'ComputeStatistics',
]
stages_dict = {k: i for i, k in enumerate(stages)}
commands = [
    'PrepareTrainingData.py',
    'TrainNet.py',
    'TestNet.py',
    'ComputeStatistics.py',
]

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('from_', metavar='from', nargs='?')
    p.add_argument('to', nargs='?')
    p.add_argument('--base')
    p.add_argument('--outputDir', '-o')
    return p.parse_args()

def dispatch():
    a = parse_args()
    if a.to is None:
        a.to = a.from_

    if a.from_ is None:
        raise NotImplementedError

    if any(tf not in stages_dict for tf in (a.to, a.from_)):
        raise ValueError('to and from must be one of the following: ' + ' '.join(stages))

    odr = a.outputDir or script_params['OUTPUT_DATA_ROOT']
    os.mkdir(odr)

    script_dir = os.path.dirname(__file__)
    source = os.path.join(odr, 'source')
    shutil.copytree(script_dir, source)

    sp = script_params.copy()
    sp['OUTPUT_DATA_ROOT'] = '..'
    with open(os.path.join(source, 'params.json'), 'w') as f:
        json.dump(sp, f, indent=4, separators=(',', ': '))

    def setup_links():
        if a.from_ == 'PrepareTrainingData':
            if a.base is not None:
                raise ValueError("--base must not be passed if starting from beginning")
            return

        if a.base is None:
            raise ValueError("--base required if not starting from beginning")

        symlink_entries_through(a.base, odr, 'training', 'Net_TrainData', 'Net_ValData', *sp['TYPES'].values())
        os.mkdir(os.path.join(odr, 'testing'))
        symlink_entries_through(a.base, odr,
                                *(os.path.join('testing', x)
                                  for x in os.listdir(os.path.join(a.base, 'testing'))
                                  if x.endswith('.mha')))

        if a.from_ == 'TrainNet':
            return

        symlink_entries_through(a.base, odr, 'NetProto')

        if a.from_ == 'TestNet':
            return

        symlink_entries_through(a.base, odr, 'testing/cnn')

    setup_links()

    with open(os.path.join(source, 'script'), 'a') as f:
        f.write(repr(sys.argv)+'\n')

        for i in range(stages_dict[a.from_], stages_dict[a.to] + 1):
            check_call([os.path.abspath(os.path.join(source, commands[i]))], cwd=source)
            f.write(repr([os.path.join('.', commands[i])])+'\n')

if __name__=='__main__':
    dispatch()
