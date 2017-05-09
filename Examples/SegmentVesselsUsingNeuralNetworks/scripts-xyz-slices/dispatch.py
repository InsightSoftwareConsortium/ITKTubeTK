#!/usr/bin/env python

from __future__ import print_function

import argparse
from collections import OrderedDict
import json
import os
import shutil
from subprocess import call
import sys

stages = [
    'TODO starting value',
    'TrainNet',
    'TestNet',
]
stages_dict = {k: i for i, k in enumerate(stages)}
commands = [
    None,
    'TrainNet.py',
    'TestNet.py',
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

    if a.from_ != 'TODO starting value' and a.base is None:
        raise ValueError("--base required if not starting from beginning")

    if a.from_ == 'TrainNet':
        for x in 'Net_TrainData', 'Net_ValData':
            symlink_entries_through(a.base, odr, x)
    elif a.from_ == 'TestNet':
        symlink_entries_through(a.base, odr, 'NetProto')
        os.mkdir(os.path.join(odr, 'testing'))
        symlink_entries_through(a.base, odr,
                                *(os.path.join('testing', x)
                                  for x in os.listdir(os.path.join(a.base, 'testing'))
                                  if x.endswith('.mha')))
    else:
        raise NotImplementedError

    with open(os.path.join(source, 'script'), 'a') as f:
        f.write(repr(sys.argv)+'\n')

        for i in range(stages_dict[a.from_], stages_dict[a.to] + 1):
            call([os.path.join(source, commands[i])], cwd=source)
            f.write(repr([commands[i]])+'\n')

def symlink_entries_through(source, dest, *args):
    for arg in args:
        symlink_through(os.path.join(source, arg), os.path.join(dest, arg))

def symlink_through(source, dest):
    """Create a symlink at dest to source, or what source links to if
    source is itself a symlink.

    """
    if os.path.islink(source):
        source = readlink_absolute(source)
    os.symlink(os.path.relpath(source, os.path.dirname(dest)), dest)

def readlink_absolute(path):
    rl = os.readlink(path)
    if not os.path.isabs(rl):
        rl = os.path.normpath(os.path.join(os.path.dirname(path), rl))
    return rl

if __name__=='__main__':
    dispatch()
