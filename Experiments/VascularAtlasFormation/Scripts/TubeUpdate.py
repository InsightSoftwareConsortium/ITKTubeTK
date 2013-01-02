import os, sys, glob, fileinput
from optparse import OptionParser
import shutil as sh
import tempfile


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = OptionParser()
    parser.add_option( "-l", "--list",
        help="List file of .tre files to be mofified\n" )

    (options, args) = parser.parse_args()

    if ( options.list is None ):
        print "Error: Need to specify file with .tre filenames!"
        return -1

    with open( options.list ) as f:
        content = f.readlines()

    print "Processing %d .tre files" % len( content )

    for file in content:
        file = file.rstrip()

        fd, temp_file_path = tempfile.mkstemp()
        with open( file ) as a:
            data = a.readlines()

        found = False
        for line in data:
            if not "ElementSpacing" in line:
                os.write(fd, line )
                continue
            if found:
                os.write(fd, line )
            else:
                os.write(fd, "ElementSpacing = 1 1 1\n" )
                found = True

        fileinput.close()
        os.close(fd)

        sh.copyfile(temp_file_path, file )
        os.remove(temp_file_path)


if __name__ == "__main__":
    sys.exit( main() )
