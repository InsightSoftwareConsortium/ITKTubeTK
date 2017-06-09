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
sys.path.append( os.path.join(TubeTK_BUILD_DIR, 'PythonModules') )

def VesselTubeToNumPyTest(tubes, baseline_array):
    import numpy as np
    from tubetk.numpy import tubes_from_file

    array = tubes_from_file(tubes)
    print(array.dtype)
    print(array)

    baseline = np.load(baseline_array)

    all_fields_close = True
    for field in baseline.dtype.fields.iterkeys():
        if not np.allclose(array[field], baseline[field]):
            all_fields_close = False
            print('The array field: ' + field + ' does not match!')

    return all_fields_close

def PyQtGraphTubesAsCirclesTest(tube_file, screenshot):
    import pyqtgraph as pg
    import pyqtgraph.opengl as gl
    from tubetk.numpy import tubes_from_file
    from tubetk.pyqtgraph import tubes_as_circles

    tubes = tubes_from_file(tube_file)
    subsample = 30
    tubes = tubes[::subsample]

    # Setup QApplication
    qapp = pg.mkQApp()

    view = gl.GLViewWidget()
    view.setCameraPosition(distance=200)
    view.orbit(0., 20.)
    view.setWindowTitle('PyQtGraphTubePointsTest')
    view.show()

    x_grid = gl.GLGridItem()
    grid_scale = 20
    x_grid.rotate(90, 0, 1, 0)
    x_grid.translate(-10 * grid_scale, 0, 0)
    y_grid = gl.GLGridItem()
    y_grid.rotate(90, 1, 0, 0)
    y_grid.translate(0, -10 * grid_scale, 0)
    z_grid = gl.GLGridItem()
    z_grid.translate(0, 0, -10 * grid_scale)
    for grid in x_grid, y_grid, z_grid:
        grid.scale(grid_scale, grid_scale, grid_scale)
        view.addItem(grid)

    circles = tubes_as_circles(tubes)
    circles_mesh = gl.GLMeshItem(meshdata=circles, smooth=False)
    view.addItem(circles_mesh)

    #return qapp.exec_()
    view.paintGL()
    framebuffer = view.grabFrameBuffer()
    return framebuffer.save(screenshot)

if __name__ == '__main__':
    usage = 'Usage: ' + sys.argv[0] + \
            ' <TestName> [TestArg1 TestArg2 ...  TestArgN]'
    if len(sys.argv) < 2:
        print usage
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
