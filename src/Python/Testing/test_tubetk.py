##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

"""Tests for the TubeTK tubetk Python package."""

import os
import sys

TubeTK_BUILD_DIR=None
if 'TubeTK_BUILD_DIR' in os.environ:
    TubeTK_BUILD_DIR = os.environ['TubeTK_BUILD_DIR']
else:
    print( 'TubeTK_BUILD_DIR not found!' )
    print( '  Set environment variable' )
    sys.exit(1)
if not os.path.exists(TubeTK_BUILD_DIR):
    print( 'TubeTK_BUILD_DIR set but directory not found!' )
    print( '  TubeTK_BUILD_DIR = ' + TubeTK_BUILD_DIR )
    sys.exit(1)

sys.path.append( os.path.join(TubeTK_BUILD_DIR, 'lib') )
sys.path.append( os.path.join(TubeTK_BUILD_DIR, 'Wrapping/Generators/Python') )
sys.path.append( os.path.join(TubeTK_BUILD_DIR, 'python') )

def VesselTubeToNumPyTest(tubes, baseline_array):
    import numpy as np
    from tubetk.numpy import tubes_from_file

    print( 'Reading file: '+tubes )
    array = tubes_from_file(tubes)
    print(array.dtype)
    print(array)

    baseline = np.load(baseline_array)

    all_fields_close = True
    for field in baseline.dtype.fields.keys():
        print( "Testing field = " + field )
        if not np.allclose(array[field], baseline[field]):
            all_fields_close = False
            print('The array field: ' + field + ' does not match!')
            print('array = '+str(array[field]))
            print('baseline = '+str(baseline[field]))

    return all_fields_close

if __name__ == '__main__':
    usage = 'Usage: ' + sys.argv[0] + \
            ' <TestName> [TestArg1 TestArg2 ...  TestArgN]'
    if len(sys.argv) < 2:
        print( usage )
        sys.exit(1)

    test_to_call = locals()[sys.argv[1]]
    # Pass arguments without an '=' as args.

    args = [arg for arg in sys.argv[2:] if arg.find('=') == -1]
    args = tuple(args)
    # Pass arguments with an '=' as kwargs.
    kwargs = [arg.split('=') for arg in sys.argv[2:] if arg.find('=') != -1]
    kwargs = dict(kwargs)

    if not test_to_call(*args, **kwargs):
        sys.exit(1)
    else:
        sys.exit(0)
