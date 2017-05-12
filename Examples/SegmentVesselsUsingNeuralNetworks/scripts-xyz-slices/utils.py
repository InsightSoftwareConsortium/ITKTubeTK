import json
import numpy as np
import os
import sys
import shutil

import sqlite3

class Logger(object):

    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass


with open(os.path.join(os.path.dirname(__file__), 'params.json')) as f:
    script_params = json.load(f)

def ensureDirectoryExists(path):
    """Create the directory named by path and any necessary parents if it
    doesn't exist.

    """
    if not os.path.exists(path):
        os.makedirs(path)


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

def open_sqlite3_db(dir):
    """Open the sqlite3 database contained in dir.  We use "data.sqlite3"."""
    return sqlite3.connect(os.path.join(dir, 'data.sqlite3'))

def scale_net_input_data(data):
    """Rescale data, presumably obtained from an 8-bit grayscale image, to
    the range [0, 1] for feeding into the network.

    """
    return data / 255.

def extractPatch(im, indices):
    """Return a patch extracted from im at indices.

    Currently, this means returning an array of size (2W+1 x 2W+1 x
    3), where W is the patch radius.

    """
    w = script_params['PATCH_RADIUS']
    # Return N (N-1)-dimensional slices
    return np.stack((im[tuple(np.s_[x - w : x + w + 1] if i != j else x
                              for j, x in enumerate(indices))]
                     for i in range(len(indices))),
                    axis=-1)

def separateChannels(im):
    """Transpose an ...xN image into an Nx...x1 image"""
    return np.moveaxis(im, -1, 0)[..., np.newaxis]

def prepareInputArray(im):
    """Convert a Bx...xC array to a list of C Bx...x1 arrays, converting
    data types appropriately in the process.

    """
    return list(scale_net_input_data(separateChannels(im)))
