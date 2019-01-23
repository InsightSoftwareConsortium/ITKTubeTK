#!/usr/bin/env python
from __future__ import with_statement

# Verify that each class's PrintSelf prints all ivars that have
# a Set/Get.
#
#
# Verbose Levels:
#
#  0:  No extra printing
#  1:  Print basic extra information
#  2:  Print lots of details
#
VERBOSE = 2

import glob
import os
import re
import sys

class ClassErrors(object):
    """Keep a running count of PrintSelf related errors.  To analyze an
    additional file, call check_header_file."""
    def __init__(self):
        self.missing_printselfs = 0
        self.missing_superclass_calls = 0
        self.missing_ivars = 0


    def get_total_defects(self):
        return self.missing_printselfs + self.missing_superclass_calls \
                + self.missing_ivars


    def check_header_file(self, fname):
        """Check the header file for the appropriate PrintSelf declaration and
        the associated implementation file for a complete implementation."""
        if VERBOSE >= 2:
            print('\nProcessing file: ' + fname)

        class_name, superclass_name, printself_found, ivars = self._analyze_header(fname)

        if class_name and len(ivars) > 0 and not printself_found:
            self.missing_printselfs += 1
            if VERBOSE >= 1:
                print('PrintSelf declaration missing for ' + class_name)

        if class_name and printself_found:
            class_path = fname[:-2]
            if os.path.exists(class_path + '.hxx'):
                self._check_printself(class_path + '.hxx', ivars,
                        superclass_name)
            elif os.path.exists(class_path + '.cxx'):
                self._check_printself(class_path + '.cxx', ivars,
                        superclass_name)
            else:
                self._check_printself(class_path + '.h', ivars,
                        superclass_name)

    def _analyze_header(self, fname):
        """Determine the class name, superclass name, if PrintSelf is declared
        and the ivars."""
        class_name = None
        superclass_name = None
        printself_found = False
        ivars = set()

        with open(fname, 'r') as fileid:
            printself_re = re.compile(r'PrintSelf.*[(]')
            class_re = re.compile(r'class .*_EXPORT\s+([^ \n:]+)')
            superclass_re = re.compile(r'public\s+([^:<\n {]+)')
            set_get_macro_re = re.compile(r'itk[SG]et(?!StaticConst)\w*Macro\s*[(]\s*([^,)\s]+)')

            for line in fileid.readlines():
                if printself_re.search(line):
                    printself_found = True
                    continue

                class_match = class_re.search(line)
                if class_match:
                    class_name = class_match.group(1)
                    if VERBOSE >= 2:
                        print('    Class name: ' + class_name)

                superclass_match = superclass_re.search(line)
                if superclass_match:
                    superclass_name = superclass_match.group(1)
                    if VERBOSE >= 2:
                        print('    Superclass name: ' + superclass_name)
                    continue

                if class_name:
                    set_get_match = set_get_macro_re.search(line)
                    if set_get_match:
                        ivars.add(set_get_match.group(1))

        return class_name, superclass_name, printself_found, ivars


    def _check_printself(self, fname, ivars, superclass_name):
        """Check that the PrintSelf implementation calls the superclass
        PrintSelf and that all ivars are printed."""
        if VERBOSE >= 2:
            print('    Checking PrintSelf implementation in file: ' + fname)

        # State variables
        inside_printself = False
        first_curly_open_found = False
        curly_open = 0
        curly_close = 0
        # Result variables
        superclass_called = False
        ivars_printed = {}
        for ivar in ivars:
            ivars_printed[ivar] = False
        with open(fname, 'r') as fileid:
            printself_re = re.compile(r'PrintSelf.*[(]')
            if superclass_name:
                superclass_re = re.compile('(' + superclass_name + \
                    r')|(Superclass)::PrintSelf')
            else:
                superclass_re = re.compile('Superclass::PrintSelf')
            for line in fileid.readlines():
                if not inside_printself and printself_re.search(line):
                    inside_printself = True

                if inside_printself:
                    if not first_curly_open_found:
                        curly_open = line.count('{')
                        curly_close = line.count('}')
                    if curly_open > 0:
                        first_curly_open_found = True

                    if first_curly_open_found:
                        if superclass_re.search(line):
                            superclass_called = True

                        for ivar in ivars:
                            if line.find('m_' + ivar) > -1 or \
                                    line.find('Get' + ivar) > -1:
                                ivars_printed[ivar] = True

                        curly_open = curly_open + line.count('{')
                        curly_close = curly_close + line.count('}')
                        if curly_open == curly_close:
                            break

        if superclass_name and not superclass_called:
            if VERBOSE >= 1:
                print('     Superclass PrintSelf not called!')
            self.missing_superclass_calls += 1

        for ivar, found in ivars_printed.items():
            if VERBOSE >= 1:
                prefix = '        ' + ivar + ' printed:'
                print(prefix.ljust(60) + str(found))
            if not found:
                self.missing_ivars += 1


def check_directories(header_directories):
    """Check the header files in the given directories for defects and return
    the number of defects."""

    class_errors = ClassErrors()
    for directory in header_directories:
        itk_header_files = glob.glob(os.path.join(directory, 'itk*.h'))
        for header in itk_header_files:
            class_errors.check_header_file(header)

    total_defects = class_errors.get_total_defects()
    if VERBOSE >= 1:
        print('\n\nTotal defects: ' + str(total_defects))
    return total_defects


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: PrintSelfCheck.py <source directory 1> <source directory 2> ...')
        sys.exit(1)

    total_defects = check_directories(sys.argv[1:])
    sys.exit(total_defects)
