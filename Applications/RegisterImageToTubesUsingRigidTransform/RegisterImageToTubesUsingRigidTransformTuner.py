#!/usr/bin/env python

"""Tune the free parameters, debug, and characterize analysis performed by
RegisterImageToTubesUsingRigidTransform."""

import argparse
import json
import subprocess
import sys

import numpy as np
import pyqtgraph as pg
import pyqtgraph.console
import pyqtgraph.dockarea
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore, QtGui
import SimpleITK as sitk

from tubetk.pyqtgraph import tubes_as_circles
from tubetk.numpy import tubes_from_file


class RegistrationTuner(QtGui.QMainWindow):

    def __init__(self, config, config_filename):
        super(RegistrationTuner, self).__init__()

        self.config = config
        self.config_filename = config_filename
        self.initializeUI()

        self.iteration = 0
        self.number_of_iterations = 0

    def initializeUI(self):
        self.resize(1024, 768)
        self.setWindowTitle('Registration Tuner')

        exit_action = QtGui.QAction('&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(QtGui.QApplication.instance().quit)
        self.addAction(exit_action)

        self.dock_area = pg.dockarea.DockArea()
        Dock = pyqtgraph.dockarea.Dock

        run_button = QtGui.QPushButton('Run Analysis')
        run_button.setToolTip('Run the analysis with ' +
                              'current parameter settings.')
        run_button.resize(run_button.sizeHint())
        QtCore.QObject.connect(run_button, QtCore.SIGNAL('clicked()'),
                               self.run_analysis)
        run_action = QtGui.QAction('&Run Analysis', self)
        run_action.setShortcut('Ctrl+R')
        run_action.setStatusTip('Run Analysis')
        QtCore.QObject.connect(run_action, QtCore.SIGNAL('triggered()'),
                               self.run_analysis)
        self.addAction(run_action)

        console_namespace = {'pg': pg, 'np': np, 'config': self.config}
        console_text = """
This is an interactive Python console.  The numpy and pyqtgraph modules have
already been imported as 'np' and 'pg'.  The parameter configuration is
available as 'config'.
"""
        self.console = pyqtgraph.console.ConsoleWidget(
            namespace=console_namespace,
            text=console_text)

        run_console_layout = QtGui.QVBoxLayout()
        run_console_layout.addWidget(run_button)
        run_console_layout.addWidget(self.console)
        run_console = QtGui.QWidget()
        run_console.setLayout(run_console_layout)
        run_console_dock = Dock("Console", size=(500, 300))
        run_console_dock.addWidget(run_console)
        self.dock_area.addDock(run_console_dock, 'right')

        self.image_tubes = gl.GLViewWidget()
        self.image_tubes.setCameraPosition(distance=200)
        x_grid = gl.GLGridItem()
        grid_scale = 20
        x_grid.rotate(90, 0, 1, 0)
        x_grid.translate(-10 * grid_scale, 0, 0)
        y_grid = gl.GLGridItem()
        y_grid.rotate(90, 1, 0, 0)
        y_grid.translate(0, -10 * grid_scale, 0)
        z_grid = gl.GLGridItem()
        z_grid.translate( 0, 0, -10 * grid_scale)
        for grid in x_grid, y_grid, z_grid:
            grid.scale(grid_scale, grid_scale, grid_scale)
            self.image_tubes.addItem(grid)
        axis = gl.GLAxisItem()
        axis_scale = 20
        axis.scale(axis_scale, axis_scale, axis_scale)
        self.image_tubes.addItem(axis)
        image_tubes_dock = Dock("Image and Tubes", size=(640, 480))
        image_tubes_dock.addWidget(self.image_tubes)
        self.dock_area.addDock(image_tubes_dock, 'left')

        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_volume = io_params[0]['Value']
        image = sitk.ReadImage(str(input_volume))
        #image_plane0 = self._image_plane(image, 0)
        #self.image_tubes.addItem(image_plane0)
        image_plane1 = self._image_plane(image, 1)
        self.image_tubes.addItem(image_plane1)
        image_plane2 = self._image_plane(image, 2)
        self.image_tubes.addItem(image_plane2)
        input_vessel = io_params[1]['Value']
        tubes = tubes_from_file(input_vessel)
        if self.config.has_key('SubSampleTubeTree'):
            sampling = self.config['SubSampleTubeTree']['Sampling']
            tubes = tubes[::sampling]
        circles = tubes_as_circles(tubes)
        circles_mesh = gl.GLMeshItem(meshdata=circles, smooth=False)
        self.image_tubes.addItem(circles_mesh)

        self.setCentralWidget(self.dock_area)
        self.show()

    def run_analysis(self):
        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_volume = io_params[0]['Value']
        input_vessel = io_params[1]['Value']
        output_volume = io_params[2]['Value']
        subprocess.check_call([config['Executables']['Analysis'],
                              '--parameterstorestore', self.config_filename,
                              input_volume,
                              input_vessel,
                              output_volume])

        self.set_iteration(0)

    def set_iteration(self, iteration):
        """Set the iteration to visualize."""
        if iteration > self.number_of_iterations:
            raise ValueError("Invalid iteration")
        self.iteration = iteration

    @staticmethod
    def _image_plane(image, direction):
        """Create an image plane Item from the center plane in the given
        direction for the given SimpleITK Image."""
        size = image.GetSize()
        index = size[direction] / 2
        image_content = sitk.GetArrayFromImage(image)
        if direction == 0:
            plane = image_content[index, :, :]
        elif direction == 1:
            plane = image_content[:, index, :]
        else:
            plane = image_content[:, :, index]
        texture = pg.makeRGBA(plane)[0]
        image_item = gl.GLImageItem(texture)
        spacing = image.GetSpacing()
        if direction == 0:
            image_item.translate(0, 0, spacing[2] * index)
        elif direction == 1:
            image_item.rotate(-90, 0, 1, 0)
            image_item.rotate(-90, 0, 0, 1)
            image_item.translate(0, spacing[1] * index, 0)
        else:
            image_item.rotate(-90, 0, 1, 0)
            image_item.translate(spacing[0] * index, 0, 0)
        return image_item

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=argparse.FileType('r'),
                        help='Configuration file for tuning analysis')
    parser.add_argument('--non-interactive', '-n', action='store_true',
                        help='Run the analysis non-interactively.')
    args = parser.parse_args()
    config_file = args.configuration

    config = json.load(config_file)

    app = pg.mkQApp()
    tuner = RegistrationTuner(config, config_file.name)
    if(args.non_interactive):
        tuner.run_analysis()
    else:
        sys.exit(app.exec_())
