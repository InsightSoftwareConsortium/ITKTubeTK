import os
import sys
import shutil

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


# Copy infile to outFile and create dirs if not present
def copy(inFile, outFile):

    # create path if it doesnt exist
    out_path = os.path.dirname(outFile)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # copy file
    shutil.copyfile(inFile, outFile)