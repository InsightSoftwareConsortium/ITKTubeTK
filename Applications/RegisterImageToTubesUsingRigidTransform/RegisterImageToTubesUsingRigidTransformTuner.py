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
import tables

from tubetk.pyqtgraph import tubes_as_circles
from tubetk.numpy import tubes_from_file


class RegistrationTuner(QtGui.QMainWindow):

    def __init__(self, config, config_filename):
        super(RegistrationTuner, self).__init__()

        self.config = config
        self.config_filename = config_filename

        self.iteration = 0
        self.number_of_iterations = 0

        self.image_planes = []
        self.subsampled_tubes = None
        self.input_image = None
        self.tubes_circles = {}
        self.dock_area = None
        self.iteration_slider = None
        self.iteration_spinbox = None
        self.progression = None

        self.initializeUI()

    def initializeUI(self):
        self.resize(1024, 768)
        self.setWindowTitle('Registration Tuner')

        exit_action = QtGui.QAction('&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(QtGui.QApplication.instance().quit)
        self.addAction(exit_action)

        iteration_dock_area_layout = QtGui.QVBoxLayout()
        iteration_dock_area = QtGui.QWidget()
        iteration_dock_area.setLayout(iteration_dock_area_layout)
        self.dock_area = pg.dockarea.DockArea()
        iteration_dock_area_layout.addWidget(self.dock_area, stretch=1)
        iteration_layout = QtGui.QHBoxLayout()
        iteration_widget = QtGui.QWidget()
        iteration_widget.setLayout(iteration_layout)
        iteration_dock_area_layout.addWidget(iteration_widget)
        iteration_layout.addWidget(QtGui.QLabel('Iteration:'))
        iteration_layout.addSpacing(1)
        self.iteration_spinbox = QtGui.QSpinBox()
        self.iteration_spinbox.setMinimum(0)
        self.iteration_spinbox.setMaximum(self.number_of_iterations)
        QtCore.QObject.connect(self.iteration_spinbox,
                               QtCore.SIGNAL('valueChanged(int)'),
                               self._iteration_spinbox_changed)
        iteration_layout.addWidget(self.iteration_spinbox)
        self.iteration_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.iteration_slider.setMinimum(0)
        self.iteration_slider.setMaximum(self.number_of_iterations)
        self.iteration_slider.setTickPosition(QtGui.QSlider.TicksBelow)
        QtCore.QObject.connect(self.iteration_slider,
                               QtCore.SIGNAL('valueChanged()'),
                               self._iteration_slider_changed)
        QtCore.QObject.connect(self.iteration_slider,
                               QtCore.SIGNAL('sliderReleased()'),
                               self._iteration_slider_changed)
        iteration_layout.addWidget(self.iteration_slider, stretch=1)

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

        import pprint
        console_namespace = {'pg': pg,
                             'np': np,
                             'config': self.config,
                             'tuner': self,
                             'pp': pprint.pprint}
        console_text = """
This is an interactive Python console.  The numpy and pyqtgraph modules have
already been imported as 'np' and 'pg'.  The parameter configuration is
available as 'config'.  The RegistrationTuner instance is available as 'tuner'.
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
        if 'Visualization' in self.config:
            visualization = self.config['Visualization']
            if 'ImageTubes' in visualization:
                imageTubes = visualization['ImageTubes']
                if 'CameraPosition' in imageTubes:
                    # The current position can be obtained with
                    # tuner.image_tubes.opts
                    cameraPosition = imageTubes['CameraPosition']
                    it = self.image_tubes
                    it.setCameraPosition(distance=cameraPosition['distance'],
                                         elevation=cameraPosition['elevation'],
                                         azimuth=cameraPosition['azimuth'])
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
        self.input_image = sitk.ReadImage(str(input_volume))
        self.add_image_planes()
        input_tubes = io_params[1]['Value']
        self.subsampled_tubes = input_tubes
        if 'SubSampleTubeTree' in self.config:
            sampling = self.config['SubSampleTubeTree']['Sampling']
            # TODO: do this correctly (use tempfile, remove on exit, etc)
            self.subsampled_tubes = '/tmp/tuner_subsampled_tubes.tre'
            subprocess.check_call([config['Executables']['SubSampleTubes'],
                                  '--samplingFactor', str(sampling),
                                  input_tubes, self.subsampled_tubes])
        self.add_tubes()

        self.setCentralWidget(iteration_dock_area)
        self.show()

    def add_image_planes(self):
        """Add the image planes.  The input image is smoothed."""
        image = self.input_image
        sigma = self.config['ParameterGroups'][1]['Parameters'][0]['Value']
        if sigma > 0.0:
            for ii in range(3):
                image = sitk.RecursiveGaussian(image, sigma, direction=ii)

        for plane in self.image_planes:
            self.image_tubes.removeItem(plane)

        self.image_planes = []
        for direction in 1, 2:
            self.image_planes.append(self._image_plane(image, direction))

        for plane in self.image_planes:
            self.image_tubes.addItem(plane)

    def add_tubes(self):
        """Add the transformed tubes for the visualization."""
        memory = 3
        target_iterations = tuple(range(self.iteration,
                                  self.iteration - memory,
                                  -1))

        for it in target_iterations:
            if it >= 0 and it <= self.number_of_iterations:
                if not it in self.tubes_circles:
                    tubes = tubes_from_file(self.subsampled_tubes)
                    circles = tubes_as_circles(tubes)
                    circles_mesh = gl.GLMeshItem(meshdata=circles,
                                                 smooth=False)
                    self.tubes_circles[it] = [circles, circles_mesh]
                else:
                    circles = self.tubes_circles[it][0]
                    self.image_tubes.removeItem(self.tubes_circles[it][1])
                    circles_mesh = gl.GLMeshItem(meshdata=circles,
                                                 smooth=False)
                    self.tubes_circles[it] = [circles, circles_mesh]

                self.image_tubes.addItem(circles_mesh)

        for it, circs in self.tubes_circles.items():
            if not it in target_iterations:
                self.image_tubes.removeItem(circs[1])
                self.tubes_circles.pop(it)

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
        progression_file = io_params[3]['Value']
        progression = tables.open_file(progression_file)
        self.progression = progression.root.OptimizationParameterProgression
        self._set_number_of_iterations(self.progression[-1]['Iteration'])

        self.add_image_planes()
        self.set_iteration(self.number_of_iterations)

    def _iteration_slider_changed(self):
        iteration = self.iteration_slider.value()
        self.set_iteration(iteration)

    def _iteration_spinbox_changed(self):
        iteration = self.iteration_spinbox.value()
        self.set_iteration(iteration)

    def set_iteration(self, iteration):
        """Set the iteration to visualize."""
        if iteration > self.number_of_iterations:
            raise ValueError("Invalid iteration")
        if iteration != self.iteration:
            self.iteration = iteration
            self.add_tubes()
        if iteration != self.iteration_slider.value():
            self.iteration_slider.setValue(iteration)
        if iteration != self.iteration_spinbox.value():
            self.iteration_spinbox.setValue(iteration)

    def _set_number_of_iterations(self, iterations):
        self.number_of_iterations = iterations
        self.iteration_slider.setMaximum(iterations)
        self.iteration_spinbox.setMaximum(iterations)

    def get_number_of_iterations(self):
        return self.number_of_iterations

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
