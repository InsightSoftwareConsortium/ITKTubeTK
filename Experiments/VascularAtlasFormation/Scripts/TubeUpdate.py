import os, sys, glob, fileinput
from optparse import OptionParser
import shutil as sh


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

    tmp_fid = open( "/tmp/tmpfile", "w" )
    with open( file ) as a:
      data = a.readlines()

    found = False
    for line in data:
      if not "ElementSpacing" in line:
        tmp_fid.write( line )
        continue

      if found:
        tmp_fid.write( line )
      else:
        # Assuming that the first occ. is the WRONG group spacing
        tmp_fid.write( "ElementSpacing = 1 1 1\n" )
        found = True

    fileinput.close()
    tmp_fid.close()

    sh.copyfile( "/tmp/tmpfile", file )


# Entry
if __name__ == "__main__":
  sys.exit( main() )
