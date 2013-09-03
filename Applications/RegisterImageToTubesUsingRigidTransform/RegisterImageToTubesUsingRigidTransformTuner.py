#!/usr/bin/env python

"""Tune the free parameters, debug, and characterize analysis performed by
RegisterImageToTubesUsingRigidTransform."""

import argparse
import json
import os
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
import matplotlib.cm

from tubetk.pyqtgraph import tubes_as_circles
from tubetk.numpy import tubes_from_file


class RegistrationTuner(QtGui.QMainWindow):

    def __init__(self, config):
        super(RegistrationTuner, self).__init__()

        self.config = config

        self.iteration = 0
        self.number_of_iterations = 0

        self.image_planes = {}
        self.subsampled_tubes = None
        self.input_image = None
        self.image_content = None
        self.tubes_circles = {}
        self.dock_area = None
        self.iteration_slider = None
        self.iteration_spinbox = None
        self.progression = None
        self.tubes_center = [0.0, 0.0, 0.0]
        self.translations = None
        self.progression_colors = None
        self.metric_values_plot = None
        self.image_controls = {}
        self.ultrasound_probe_origin = None

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
                               QtCore.SIGNAL('valueChanged(int)'),
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

        video_button = QtGui.QPushButton('Make Video Frames')
        video_button.resize(video_button.sizeHint())
        QtCore.QObject.connect(video_button, QtCore.SIGNAL('clicked()'),
                               self.make_video)
        # todo do properly
        self.video_frame_dir = '/tmp/registration_tuner'
        if not os.path.exists(self.video_frame_dir):
            os.makedirs(self.video_frame_dir)

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
        run_console_layout.addWidget(video_button)
        run_console_layout.addWidget(self.console)
        run_console = QtGui.QWidget()
        run_console.setLayout(run_console_layout)
        run_console_dock = Dock("Console", size=(400, 300))
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
                    # HACK: pyqtgraph needs to fix its api so this can be set
                    # with setCameraPosition
                    if 'center' in cameraPosition:
                        center = cameraPosition['center']
                        it.opts['center'] = QtGui.QVector3D(center[0],
                                                            center[1],
                                                            center[2])
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

        image_controls_dock = Dock("Image Display Controls", size=(640, 60))
        x_check = QtGui.QCheckBox("X Plane")
        x_check.setChecked(False)
        QtCore.QObject.connect(x_check,
                               QtCore.SIGNAL('stateChanged(int)'),
                               self._x_check_changed)
        self.image_controls['x_check'] = x_check
        image_controls_dock.addWidget(x_check, row=0, col=0)
        x_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.image_controls['x_slider'] = x_slider
        QtCore.QObject.connect(x_slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               self._x_slider_changed)
        image_controls_dock.addWidget(x_slider, row=0, col=1)
        y_check = QtGui.QCheckBox("Y Plane")
        y_check.setChecked(True)
        self.image_controls['y_check'] = y_check
        QtCore.QObject.connect(y_check,
                               QtCore.SIGNAL('stateChanged(int)'),
                               self._y_check_changed)
        image_controls_dock.addWidget(y_check, row=1, col=0)
        y_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.image_controls['y_slider'] = y_slider
        QtCore.QObject.connect(y_slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               self._y_slider_changed)
        image_controls_dock.addWidget(y_slider, row=1, col=1)
        z_check = QtGui.QCheckBox("Z Plane")
        z_check.setChecked(True)
        self.image_controls['z_check'] = z_check
        QtCore.QObject.connect(z_check,
                               QtCore.SIGNAL('stateChanged(int)'),
                               self._z_check_changed)
        image_controls_dock.addWidget(z_check, row=2, col=0)
        z_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.image_controls['z_slider'] = z_slider
        QtCore.QObject.connect(z_slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               self._z_slider_changed)
        image_controls_dock.addWidget(z_slider, row=2, col=1)
        self.dock_area.addDock(image_controls_dock, 'bottom', image_tubes_dock)

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

        metric_value_dock = Dock("Metric Value Inverse", size=(500, 300))
        self.metric_values_plot = pg.PlotWidget()
        self.metric_values_plot.setLabel('bottom', text='Iteration')
        self.metric_values_plot.setLabel('left', text='Metric Value Inverse')
        self.metric_values_plot.setLogMode(False, True)
        metric_value_dock.addWidget(self.metric_values_plot)
        self.dock_area.addDock(metric_value_dock, 'left', run_console_dock)

        self.setCentralWidget(iteration_dock_area)
        self.show()

    def add_image_planes(self):
        """Add the image planes.  The input image is smoothed."""
        image = self.input_image
        sigma = self.config['ParameterGroups'][1]['Parameters'][0]['Value']
        if sigma > 0.0:
            for ii in range(3):
                image = sitk.RecursiveGaussian(image, sigma, direction=ii)
        self.image_content = sitk.GetArrayFromImage(image)
        self.image_controls['x_slider'].setMaximum(
            self.image_content.shape[2] - 1)
        self.image_controls['y_slider'].setMaximum(
            self.image_content.shape[1] - 1)
        self.image_controls['z_slider'].setMaximum(
            self.image_content.shape[0] - 1)

        for plane in self.image_planes.keys():
            self.image_tubes.removeItem(self.image_planes[plane])
            self.image_planes.pop(plane)

        for direction in 'x', 'y', 'z':
            self.image_planes[direction] = self._image_plane(direction)
            checkbox = direction + '_check'
            if self.image_controls[checkbox].isChecked():
                self.image_planes[direction].setVisible(True)
            else:
                self.image_planes[direction].setVisible(False)

        for plane in self.image_planes:
            self.image_tubes.addItem(self.image_planes[plane])

    def add_ultrasound_probe_origin(self):
        """Add a sphere indicating the ultrasound probe origin."""
        with open(self.config['UltrasoundProbeGeometryFile'], 'r') as fp:
            line = fp.readline()
            line = line.strip()
            entries = line.split()
            origin = [float(xx) for xx in entries[1:]]
        if self.ultrasound_probe_origin:
            self.image_tubes.removeItem(self.ultrasound_probe_origin)
        sphere = gl.MeshData.sphere(rows=10,
                                    cols=20,
                                    radius=3.0)
        center_mesh = gl.GLMeshItem(meshdata=sphere,
                                    smooth=False,
                                    color=(1.0, 1.0, 0.3, 0.7),
                                    glOptions='translucent')
        center_mesh.translate(origin[0], origin[1], origin[2])
        self.ultrasound_probe_origin = center_mesh
        self.image_tubes.addItem(self.ultrasound_probe_origin)

    def add_tubes(self):
        """Add the transformed tubes for the visualization."""
        memory = 4
        target_iterations = tuple(range(self.iteration,
                                  self.iteration - memory,
                                  -1))

        for it in target_iterations:
            if it >= 0 and it <= self.number_of_iterations:
                alpha = np.exp(0.8*(it - self.iteration))
                center_color = (0.8, 0.1, 0.25, alpha)
                center = self.tubes_center
                if self.progression != None:
                    parameters = self.progression[it]['Parameters']
                else:
                    parameters = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                tubes = tubes_from_file(self.subsampled_tubes)
                #tubes_color = (0.2, 0.25, 0.75, alpha)
                tube_weights = 2./(1. + np.exp(-2 * tubes['Radius']))
                tube_weights = tube_weights - tube_weights.min()
                tube_weights = tube_weights / tube_weights.max()
                tubes_colors = matplotlib.cm.PuBuGn(tube_weights)
                tubes_colors[:, 3] = tube_weights**0.5 * alpha
                circles = tubes_as_circles(tubes, point_colors=tubes_colors)
                circles_mesh = gl.GLMeshItem(meshdata=circles,
                                             glOptions='translucent',
                                             smooth=False)
                if self.progression != None:
                    # TODO: need to verify that this is correct
                    circles_mesh.translate(-center[0],
                                           -center[1],
                                           -center[2])
                    r2d = 180. / np.pi
                    circles_mesh.rotate(parameters[0] * r2d, 1, 0, 0,
                                        local=False)
                    circles_mesh.rotate(parameters[1] * r2d, 0, 1, 0,
                                        local=False)
                    circles_mesh.rotate(parameters[2] * r2d, 0, 0, 1,
                                        local=False)
                    circles_mesh.translate(center[0] + parameters[3],
                                           center[1] + parameters[4],
                                           center[2] + parameters[5])
                if it in self.tubes_circles:
                    self.image_tubes.removeItem(self.tubes_circles[it][1])
                self.image_tubes.addItem(circles_mesh)

                sphere = gl.MeshData.sphere(rows=10,
                                            cols=20,
                                            radius=1.0)
                center_mesh = gl.GLMeshItem(meshdata=sphere,
                                            smooth=False,
                                            color=center_color,
                                            glOptions='translucent')
                center_mesh.translate(center[0] + parameters[3],
                                      center[1] + parameters[4],
                                      center[2] + parameters[5])
                if it in self.tubes_circles:
                    self.image_tubes.removeItem(self.tubes_circles[it][2])
                self.image_tubes.addItem(center_mesh)

                self.tubes_circles[it] = [circles,
                                          circles_mesh,
                                          center_mesh]

        for it, circs in self.tubes_circles.items():
            if not it in target_iterations:
                self.image_tubes.removeItem(circs[1])
                self.image_tubes.removeItem(circs[2])
                self.tubes_circles.pop(it)

    def add_translations(self):
        """Add a plot of the translation of the tube centers throughout the
        registration."""
        if self.translations:
            self.image_tubes.removeItem(self.translations)
        translations = np.array(self.progression[:]['Parameters'][:, 3:]) + \
            self.tubes_center
        colors = self.progression_colors
        colors[:, 3] = 0.9
        translation_points = gl.GLScatterPlotItem(pos=translations,
                                                  color=colors)
        translation_points.setGLOptions('translucent')
        self.translations = translation_points
        self.image_tubes.addItem(self.translations)

    def plot_metric_values(self):
        iterations = self.progression[:]['Iteration']
        metric_value_inverse = 1. / self.progression[:]['CostFunctionValue']
        min_metric = np.min(metric_value_inverse)
        self.metric_values_plot.plot({'x': iterations,
                                     'y': metric_value_inverse},
                                     pen={'color': (100, 100, 200, 200),
                                          'width': 2.5},
                                     fillLevel=min_metric,
                                     fillBrush=(255, 255, 255, 30),
                                     antialias=True,
                                     clear=True)
        self.metric_values_plot.plot({'x': [self.iteration],
                                     'y': [metric_value_inverse[self.iteration]]},
                                     symbolBrush='r',
                                     symbol='o',
                                     symbolSize=8.0)

    def run_analysis(self):
        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_volume = io_params[0]['Value']
        input_vessel = io_params[1]['Value']
        output_volume = io_params[2]['Value']
        # TODO: Do this properly with tempfile, cleanup, etc
        config_file = '/tmp/registration_tuner.json'
        with open(config_file, 'w') as fp:
            json.dump(self.config, fp)
        subprocess.check_call([config['Executables']['Analysis'],
                              '--parameterstorestore', config_file,
                              input_volume,
                              input_vessel,
                              output_volume])
        progression_file = io_params[3]['Value']
        progression = tables.open_file(progression_file)
        self.progression = np.array(progression.root.OptimizationParameterProgression)
        self.tubes_center = np.array(progression.root.FixedParameters)
        progression.close()
        self._set_number_of_iterations(self.progression[-1]['Iteration'])

        self.add_image_planes()
        if 'UltrasoundProbeGeometryFile' in self.config:
            self.add_ultrasound_probe_origin()
        self.add_translations()
        self.set_iteration(self.number_of_iterations)

    def make_video(self):
        for it in range(self.number_of_iterations + 1):
            filename = os.path.join(self.video_frame_dir,
                                    'frame_{0:04d}.png'.format(it))
            self.set_iteration(it)
            self.update()
            QtCore.QCoreApplication.processEvents()
            pixmap = QtGui.QPixmap.grabWindow(self.winId())
            print('Saving ' + filename)
            pixmap.save(filename, 'png')

    def _x_check_changed(self):
        if self.image_planes['x']:
            if self.image_controls['x_check'].isChecked():
                self.image_planes['x'].setVisible(True)
            else:
                self.image_planes['x'].setVisible(False)

    def _x_slider_changed(self):
        if self.image_controls['x_check'].isChecked():
            if self.image_planes['x']:
                self.image_tubes.removeItem(self.image_planes['x'])
            index = self.image_controls['x_slider'].value()
            self.image_planes['x'] = self._image_plane('x', index)
            self.image_tubes.addItem(self.image_planes['x'])

    def _y_check_changed(self):
        if self.image_planes['y']:
            if self.image_controls['y_check'].isChecked():
                self.image_planes['y'].setVisible(True)
            else:
                self.image_planes['y'].setVisible(False)

    def _y_slider_changed(self):
        if self.image_controls['y_check'].isChecked():
            if self.image_planes['y']:
                self.image_tubes.removeItem(self.image_planes['y'])
            index = self.image_controls['y_slider'].value()
            self.image_planes['y'] = self._image_plane('y', index)
            self.image_tubes.addItem(self.image_planes['y'])

    def _z_check_changed(self):
        if self.image_planes['z']:
            if self.image_controls['z_check'].isChecked():
                self.image_planes['z'].setVisible(True)
            else:
                self.image_planes['z'].setVisible(False)

    def _z_slider_changed(self):
        if self.image_controls['z_check'].isChecked():
            if self.image_planes['z']:
                self.image_tubes.removeItem(self.image_planes['z'])
            index = self.image_controls['z_slider'].value()
            self.image_planes['z'] = self._image_plane('z', index)
            self.image_tubes.addItem(self.image_planes['z'])

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
            self.plot_metric_values()
        if iteration != self.iteration_slider.value():
            self.iteration_slider.setValue(iteration)
        if iteration != self.iteration_spinbox.value():
            self.iteration_spinbox.setValue(iteration)

    def _set_number_of_iterations(self, iterations):
        self.number_of_iterations = iterations
        self.iteration_slider.setMaximum(iterations)
        self.iteration_spinbox.setMaximum(iterations)
        iterations_normalized = np.arange(0, iterations, dtype=np.float) / \
            iterations
        self.progression_colors = matplotlib.cm.summer(iterations_normalized)

    def get_number_of_iterations(self):
        return self.number_of_iterations

    def _image_plane(self, direction, index=None):
        """Create an image plane Item from the center plane in the given
        direction for the given SimpleITK Image."""
        if index is None:
            shape = self.image_content.shape
            if direction == 'x':
                index = shape[2] / 2
            elif direction == 'y':
                index = shape[1] / 2
            elif direction == 'z':
                index = shape[0] / 2
        if direction == 'x':
            plane = self.image_content[:, :, index]
        elif direction == 'y':
            plane = self.image_content[:, index, :]
        else:
            plane = self.image_content[index, :, :]
            plane = plane.transpose()
        texture = pg.makeRGBA(plane)[0]
        image_item = gl.GLImageItem(texture)
        spacing = self.input_image.GetSpacing()
        origin = self.input_image.GetOrigin()
        if direction == 'x':
            image_item.scale(spacing[2], spacing[1], 1)
            image_item.rotate(-90, 0, 1, 0)
            image_item.translate(origin[0] + spacing[2] * index,
                                 origin[1],
                                 origin[2])
        elif direction == 'y':
            image_item.scale(spacing[2], spacing[0], 1)
            image_item.rotate(-90, 0, 1, 0)
            image_item.rotate(-90, 0, 0, 1)
            image_item.translate(origin[0],
                                 origin[1] + spacing[1] * index,
                                 origin[2])
        else:
            image_item.scale(spacing[0], spacing[1], 1)
            image_item.translate(origin[0],
                                 origin[1],
                                 origin[2] + spacing[0] * index)
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
    tuner = RegistrationTuner(config)
    if(args.non_interactive):
        tuner.run_analysis()
    else:
        sys.exit(app.exec_())
