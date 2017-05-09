#!/usr/bin/env python

from __future__ import print_function

import argparse
from collections import OrderedDict
import json
import os
import shutil
from subprocess import check_call
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

    def setup_links():
        if a.from_ == 'TODO starting value':
            if a.base is not None:
                raise ValueError("--base must not be passed if starting from beginning")
            return

        if a.base is None:
            raise ValueError("--base required if not starting from beginning")

        symlink_entries_through(a.base, odr, 'Net_TrainData', 'Net_ValData')
        os.mkdir(os.path.join(odr, 'testing'))
        symlink_entries_through(a.base, odr,
                                *(os.path.join('testing', x)
                                  for x in os.listdir(os.path.join(a.base, 'testing'))
                                  if x.endswith('.mha')))

        if a.from_ == 'TrainNet':
            return

        symlink_entries_through(a.base, odr, 'NetProto')

    setup_links()

    with open(os.path.join(source, 'script'), 'a') as f:
        f.write(repr(sys.argv)+'\n')

        for i in range(stages_dict[a.from_], stages_dict[a.to] + 1):
            check_call([os.path.join(source, commands[i])], cwd=source)
            f.write(repr([os.path.join('.', commands[i])])+'\n')

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
