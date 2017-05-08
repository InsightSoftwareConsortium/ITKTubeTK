#!/usr/bin/env python

from __future__ import print_function

import argparse
from collections import OrderedDict
import json
import os
import shutil
from subprocess import call
import sys

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
    if not (a.from_ == a.to == 'TrainNet'):
        raise NotImplementedError

    with open('params.json') as f:
        # So that the output resembles the input
        script_params = json.load(f, object_pairs_hook=OrderedDict)
    odr = a.outputDir or script_params['OUTPUT_DATA_ROOT']
    os.mkdir(odr)

    source = os.path.join(odr, 'source')
    shutil.copytree('.', source)

    sp = script_params.copy()
    sp['OUTPUT_DATA_ROOT'] = '..'
    with open(os.path.join(source, 'params.json'), 'w') as f:
        json.dump(sp, f, indent=4, separators=(',', ': '))

    if a.base is None:
        raise NotImplementedError

    for x in 'Net_TrainData', 'Net_ValData':
        os.symlink(os.path.relpath(os.path.join(a.base, x), odr), os.path.join(odr, x))

    with open(os.path.join(source, 'script'), 'a') as f:
        f.write(repr(sys.argv)+'\n')

        call([os.path.join(source, 'TrainNet.py')], cwd=source)
        f.write(repr(['TrainNet.py'])+'\n')

if __name__=='__main__':
    dispatch()
