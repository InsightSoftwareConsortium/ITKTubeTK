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


def ensureDirectoryExists(path):
    """Create the directory named by path and any necessary parents if it
    doesn't exist.

    """
    if not os.path.exists(path):
        os.makedirs(path)

# Copy infile to outFile and create dirs if not present
def copy(inFile, outFile):

    # create path if it doesnt exist
    out_path = os.path.dirname(outFile)
    ensureDirectoryExists(out_path)

    # copy file
    shutil.copyfile(inFile, outFile)

def open_sqlite3_db(dir):
    """Open the sqlite3 database contained in dir.  We use "data.sqlite3"."""
    return sqlite3.connect(os.path.join(dir, 'data.sqlite3'))

def scale_net_input_data(data):
    """Rescale data, presumably obtained from an 8-bit grayscale image, to
    the range [0, 1] for feeding into the network.

    """
    return data / 255.
