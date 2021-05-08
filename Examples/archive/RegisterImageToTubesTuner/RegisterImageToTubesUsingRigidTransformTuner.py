#!/usr/bin/env python

"""Tune the free parameters, debug, and characterize analysis performed by
RegisterImageToTubesUsingRigidTransform."""

import argparse
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pyqtgraph as pg
import pyqtgraph.console
import pyqtgraph.dockarea
import pyqtgraph.opengl as gl
from pyqtgraph.widgets.HistogramLUTWidget import HistogramLUTWidget
from pyqtgraph.graphicsItems.ImageItem import ImageItem
from pyqtgraph.Qt import QtCore, QtGui
import itk
import tables
import matplotlib.cm

from tubetk.pyqtgraph import tubes_as_circles
from tubetk.numpy import tubes_from_file


class HistogramWidget(HistogramLUTWidget):
    """Like HistogramLUTWidget, but hide the LUT."""

    def __init__(self, *args, **kwargs):
        super(HistogramWidget, self).__init__(*args, **kwargs)
        self.gradient.setVisible(False)


class SharedIterationLogic(QtCore.QObject):
    """Shares the current iteration between multiple registration instances."""

    _iteration = 0

    def __init__(self, parent=None):
        super(SharedIterationLogic, self).__init__(parent)

    def get_iteration(self):
        return self._iteration

    def set_iteration(self, value):
        """Set the iteration to visualize."""
        if value != self._iteration:
            self._iteration = value
            self.iteration_changed.emit(value)

    iteration = property(get_iteration,
                         set_iteration,
                         doc='Currently examined iteration.')

    iteration_changed = QtCore.pyqtSignal(int, name='iterationChanged')


class RegistrationTunerLogic(QtCore.QObject):
    """Controls the business logic for registration tuning."""

    _iteration = 0
    _number_of_iterations = 0
    _iterations_colors = None
    _progression = None
    _tubes_center = [0.0, 0.0, 0.0]
    _subsampled_tubes = None
    _metric_image = None
    _tube_point_weights = None
    _tuner_clones = []

    def __init__(self, config, shared_iteration_logic=None,
                 image_display_logic=None, parent=None):
        super(RegistrationTunerLogic, self).__init__(parent)
        self.config = config

        if shared_iteration_logic is None:
            self.shared_iteration_logic = SharedIterationLogic(self)
        else:
            self.shared_iteration_logic = shared_iteration_logic
        self.shared_iteration_logic.iteration_changed.connect(self.set_iteration)
        self.iteration_changed.connect(self.shared_iteration_logic.set_iteration)

        if image_display_logic is None:
            image_display_logic = ImageDisplayLogic()
        self.image_display_logic = image_display_logic

        TemporaryFile = tempfile.NamedTemporaryFile
        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_tubes = io_params[1]['Value']
        if 'SubSampleTubeTree' in self.config:
            sampling = self.config['SubSampleTubeTree']['Sampling']
            subsampled_tubes_fp = TemporaryFile(suffix='SubsampledTubes.tre',
                                                delete=False)
            self._subsampled_tubes = subsampled_tubes_fp.name
            subsampled_tubes_fp.close()
            subprocess.check_call([config['Executables']['SubSampleTubes'],
                                  '--samplingFactor', str(sampling),
                                  input_tubes, self._subsampled_tubes])
        else:
            self._subsampled_tubes = input_tubes

        metric_image_fp = TemporaryFile(suffix='MetricImage.mha',
                                        delete=False)
        self._metric_image = metric_image_fp.name
        metric_image_fp.close()
        tube_point_weights_fp = TemporaryFile(suffix='TubePointWeights.json',
                                              delete=False)
        self._tube_point_weights = tube_point_weights_fp.name
        tube_point_weights_fp.close()
        self.config['TubePointWeightsFile'] = self._tube_point_weights

    def initialize(self):
        self.video_frame_dir = tempfile.mkdtemp()

    def __enter__(self):
        self.initialize()
        return self

    def close(self):
        shutil.rmtree(self.video_frame_dir)
        if 'SubSampleTubeTree' in self.config:
            os.remove(self.subsampled_tubes)
        os.remove(self._metric_image)
        os.remove(self._tube_point_weights)

    def __exit__(self, type, value, traceback):
        self.close()

    def create_new_tuner(self):
        """Create a new tuner with a copy of the config, and the same shared
        set of target iterations."""
        new_config = copy.deepcopy(self.config)
        tuner_clone = RegistrationTuner(new_config,
                                        self.shared_iteration_logic,
                                        self.image_display_logic)
        self._tuner_clones.append(tuner_clone)

    def get_subsampled_tubes(self):
        return self._subsampled_tubes

    subsampled_tubes = property(get_subsampled_tubes,
                                doc='Optionally subsampled input tubes ' +
                                'filename')

    def get_metric_image(self):
        return self._metric_image

    metric_image = property(get_metric_image,
                            doc='Metric image filename')

    iteration_changed = QtCore.pyqtSignal(int, name='iterationChanged')

    def get_iteration(self):
        return self._iteration

    def set_iteration(self, value):
        """Set the iteration to visualize."""
        if value > self.number_of_iterations:
            return
        if value != self._iteration:
            self._iteration = value
            self.iteration_changed.emit(value)

    iteration = property(get_iteration,
                         set_iteration,
                         doc='Currently examined iteration.')

    number_of_iterations_changed = \
        QtCore.pyqtSignal(int,
                          name='numberOfIterationsChanged')

    def get_number_of_iterations(self):
        return self._number_of_iterations

    def set_number_of_iterations(self, value):
        if value != self._number_of_iterations:
            self._number_of_iterations = value
            iterations_normalized = np.arange(0,
                                              self._number_of_iterations,
                                              dtype=np.float) / \
                self._number_of_iterations
            colors = matplotlib.cm.summer(iterations_normalized)
            colors[:, 3] = 0.95
            self._iterations_colors = colors
            self.number_of_iterations_changed.emit(value)

    number_of_iterations = property(get_number_of_iterations,
                                    set_number_of_iterations,
                                    doc='Iterations in the optimization')

    def get_iterations_colors(self):
        return self._iterations_colors

    def get_progression(self):
        return self._progression

    progression = property(get_progression,
                           doc='Parameter optimization progression')

    def get_tubes_center(self):
        return self._tubes_center

    tubes_center = property(get_tubes_center,
                            doc='Center of the tubes')

    analysis_run = QtCore.pyqtSignal(name='analysisRun')

    def run_analysis(self):
        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_volume = io_params[0]['Value']
        input_vessel = io_params[1]['Value']
        output_transform = io_params[2]['Value']
        NamedTemporaryFile = tempfile.NamedTemporaryFile
        config_file = NamedTemporaryFile(suffix='TunerConfig.json',
                                         delete=False)
        json.dump(self.config, config_file)
        config_file.close()
        command = [config['Executables']['Analysis'],
                   '--deserialize', config_file.name,
                   input_volume,
                   input_vessel,
                   output_transform]
        # The next statement is to keep the unicorns happy. Without it,
        #   OSError: [Errno 9] Bad file descriptor
        with open(os.path.devnull, 'wb') as unicorn:
            subprocess.call(command)
        os.remove(config_file.name)
        progression_file = io_params[3]['Value']
        # TODO: is this being created correctly on a fresh run?
        print(progression_file)
        progression = tables.open_file(progression_file)
        self._progression = np.array(progression.root.
                                     OptimizationParameterProgression)
        self._tubes_center = np.array(progression.root.FixedParameters)
        progression.close()

        self.set_number_of_iterations(self.progression[-1]['Iteration'])
        self.analysis_run.emit()
        self.iteration = self.number_of_iterations

    metric_sampled = QtCore.pyqtSignal(name='metricSampled')

    def sample_metric(self):
        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_volume = io_params[0]['Value']
        input_vessel = io_params[1]['Value']
        NamedTemporaryFile = tempfile.NamedTemporaryFile
        config_file = NamedTemporaryFile(suffix='TunerConfig.json',
                                         delete=False)
        json.dump(self.config, config_file)
        config_file.close()
        command = [config['Executables']['MetricSampler'],
                   '--deserialize', config_file.name,
                   input_volume,
                   input_vessel,
                   self.metric_image]
        # The next statement is to keep the unicorns happy. Without it,
        #   OSError: [Errno 9] Bad file descriptor
        with open(os.path.devnull, 'wb') as unicorn:
            subprocess.call(command)
        os.remove(config_file.name)
        self.metric_sampled.emit()

    def center_metric_sampling_on_converged(self):
        """Center the metric sampling based on the most recent optimized
        parameter progression."""
        progression = self.progression
        if progression == None:
            return
        last_params = progression[-1]['Parameters']
        print('last_params', last_params)
        number_of_params = len(last_params)
        metric_sampler_settings = self.config['MetricSampler']
        upper_bound = metric_sampler_settings['UpperBound']
        lower_bound = metric_sampler_settings['LowerBound']
        for ii in range(number_of_params):
            half_bound = (upper_bound[ii] - lower_bound[ii]) / 2
            center = last_params[ii]
            upper_bound[ii] = center + half_bound
            lower_bound[ii] = center - half_bound
        metric_sampler_settings['UpperBound'] = upper_bound
        metric_sampler_settings['LowerBound'] = lower_bound
        self.sample_metric()


class ImageDisplayLogic(QtCore.QObject):
    """Controls the business logic for viewing the image."""

    _planes_visible = [0, 1, 1]
    _plane_indices = [0, 0, 0]
    _number_of_planes = [1, 1, 1]

    def __init__(self, parent=None):
        super(ImageDisplayLogic, self).__init__(parent)

    plane_visibility_changed = QtCore.pyqtSignal(int, int,
                                                 name='planeVisibilityChanged')

    def get_plane_visibility(self, direction):
        return self._planes_visible[direction]

    def set_plane_visibility(self, direction, visible):
        if self._planes_visible[direction] != visible:
            self._planes_visible[direction] = visible
            self.plane_visibility_changed.emit(direction, visible)

    def set_plane_visibility_x(self, visible):
        self.set_plane_visibility(0, visible)

    def set_plane_visibility_y(self, visible):
        self.set_plane_visibility(1, visible)

    def set_plane_visibility_z(self, visible):
        self.set_plane_visibility(2, visible)

    plane_visibility = property(get_plane_visibility,
                                set_plane_visibility,
                                doc='Whether the image plane is visible')

    plane_indices_changed = QtCore.pyqtSignal(int, int,
                                              name='planeIndicesChanged')

    def get_plane_indices(self, direction):
        return self._plane_indices[direction]

    def set_plane_indices(self, direction, index):
        if self._plane_indices[direction] != index:
            self._plane_indices[direction] = index
            self.plane_indices_changed.emit(direction, index)

    def set_plane_index_x(self, index):
        self.set_plane_indices(0, index)

    def set_plane_index_y(self, index):
        self.set_plane_indices(1, index)

    def set_plane_index_z(self, index):
        self.set_plane_indices(2, index)

    plane_indices = property(get_plane_indices,
                             set_plane_indices,
                             doc='The index of the plane visualized')

    number_of_planes_changed = QtCore.pyqtSignal(int, int,
                                                 name='numberOfPlanesChanged')

    def get_number_of_planes(self, direction):
        return self._number_of_planes[direction]

    def set_number_of_planes(self, direction, size):
        if self._number_of_planes[direction] != size:
            self._number_of_planes[direction] = size
            self.number_of_planes_changed.emit(direction, size)

    number_of_planes = property(get_number_of_planes,
                                set_number_of_planes,
                                doc='Number of planes in a given direction')


class IterationWidget(QtGui.QWidget):
    """Control the currently viewed iteration."""

    def __init__(self):
        super(IterationWidget, self).__init__()
        self.initializeUI()

    def initializeUI(self):
        layout = QtGui.QHBoxLayout()
        self.setLayout(layout)

        layout.addWidget(QtGui.QLabel('Iteration:'))

        layout.addSpacing(1)

        self.spinbox = QtGui.QSpinBox()
        self.spinbox.setMinimum(0)
        self.spinbox.setMaximum(0)
        layout.addWidget(self.spinbox)

        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)
        self.slider.setTickPosition(QtGui.QSlider.TicksBelow)
        layout.addWidget(self.slider, stretch=1)

    def set_iteration(self, value):
        if value != self.slider.value():
            self.slider.setValue(value)
        if value != self.spinbox.value():
            self.spinbox.setValue(value)

    def set_number_of_iterations(self, value):
        if self.spinbox.maximum() != value:
            self.spinbox.setMaximum(value)
        if self.slider.maximum() != value:
            self.slider.setMaximum(value)


class RunConsoleWidget(QtGui.QWidget):
    """The interactive Python console and buttons to run the analysis."""

    def __init__(self, config, tuner):
        super(RunConsoleWidget, self).__init__()
        self.config = config
        self.tuner = tuner
        self.initializeUI()

    def initializeUI(self):
        layout = QtGui.QVBoxLayout()
        self.setLayout(layout)

        run_action = QtGui.QAction('&Run Analysis', self)
        run_action.setShortcut('Ctrl+R')
        run_action.setStatusTip('Run Analysis')
        QtCore.QObject.connect(run_action, QtCore.SIGNAL('triggered()'),
                               self.tuner.logic.run_analysis)
        self.tuner.addAction(run_action)

        run_button = QtGui.QPushButton('Run Analysis')
        run_button.setToolTip('Run the analysis with ' +
                              'current parameter settings.')
        run_button.resize(run_button.sizeHint())
        layout.addWidget(run_button)
        QtCore.QObject.connect(run_button,
                               QtCore.SIGNAL('clicked()'),
                               self.tuner.logic.run_analysis)

        sample_metric_action = QtGui.QAction('&Sample Metric', self)
        sample_metric_action.setShortcut('Ctrl+M')
        sample_metric_action.setStatusTip('Sample Metric')
        QtCore.QObject.connect(sample_metric_action, QtCore.SIGNAL('triggered()'),
                               self.tuner.logic.sample_metric)
        self.tuner.addAction(sample_metric_action)
        sample_metric_button = QtGui.QPushButton('Sample Metric')
        sample_metric_button.setToolTip('Sample the metric space with ' +
                              'current parameter settings.')
        sample_metric_button.resize(sample_metric_button.sizeHint())
        layout.addWidget(sample_metric_button)
        QtCore.QObject.connect(sample_metric_button,
                               QtCore.SIGNAL('clicked()'),
                               self.tuner.logic.sample_metric)

        center_metric_action = QtGui.QAction('Center Metric &Sampling', self)
        center_metric_action.setShortcut('Ctrl+S')
        center_metric_action.setStatusTip('Center Metric Sampling')
        QtCore.QObject.connect(center_metric_action, QtCore.SIGNAL('triggered()'),
                               self.tuner.logic.center_metric_sampling_on_converged)
        self.tuner.addAction(center_metric_action)
        center_metric_button = QtGui.QPushButton('Center Metric Sampling')
        center_metric_button.setToolTip('Center the metric space sampling domain ' +
                              'on the last optimized values.')
        center_metric_button.resize(center_metric_button.sizeHint())
        layout.addWidget(center_metric_button)
        QtCore.QObject.connect(center_metric_button,
                               QtCore.SIGNAL('clicked()'),
                               self.tuner.logic.center_metric_sampling_on_converged)

        new_tuner_action = QtGui.QAction('&New Tuner', self)
        new_tuner_action.setShortcut('Ctrl+N')
        new_tuner_action.setStatusTip(
            'Create a new tuner based on the current configuration.')
        QtCore.QObject.connect(new_tuner_action, QtCore.SIGNAL('triggered()'),
                               self.tuner.logic.create_new_tuner)
        self.tuner.addAction(new_tuner_action)

        new_tuner_button = QtGui.QPushButton('New Tuner')
        new_tuner_button.setToolTip(
            'Create a new tuner based on the current configuration.')
        new_tuner_button.resize(new_tuner_button.sizeHint())
        layout.addWidget(new_tuner_button)
        QtCore.QObject.connect(new_tuner_button,
                               QtCore.SIGNAL('clicked()'),
                               self.tuner.logic.create_new_tuner)

        video_button = QtGui.QPushButton('Make Video Frames')
        video_button.resize(video_button.sizeHint())
        layout.addWidget(video_button)
        QtCore.QObject.connect(video_button,
                               QtCore.SIGNAL('clicked()'),
                               self.tuner.make_video)

        import pprint
        console_namespace = {'pg': pg,
                             'np': np,
                             'config': self.config,
                             'tuner': self.tuner,
                             'pp': pprint.pprint}
        console_text = """
This is an interactive Python console.  The numpy and pyqtgraph modules have
already been imported as 'np' and 'pg'.  The parameter configuration is
available as 'config'.  The RegistrationTuner instance is available as 'tuner'.
"""
        self.console = pyqtgraph.console.ConsoleWidget(
            namespace=console_namespace,
            text=console_text)
        layout.addWidget(self.console)


class ImageTubesWidget(gl.GLViewWidget):
    """The visualization of the image and tubes in 3D space."""

    def __init__(self, logic, image_display_logic, config):
        super(ImageTubesWidget, self).__init__()
        self.logic = logic
        self.image_display_logic = image_display_logic
        self.config = config

        io_params = self.config['ParameterGroups'][0]['Parameters']
        input_volume = io_params[0]['Value']
        self.input_image = itk.ImageFileReader(FileName=str(input_volume))

        if 'Visualization' in config:
            visualization = config['Visualization']
            if 'ImageTubes' in visualization:
                imageTubes = visualization['ImageTubes']
                if 'CameraPosition' in imageTubes:
                    # The current position can be obtained with
                    # tuner.image_tubes.opts
                    position = imageTubes['CameraPosition']
                    self.setCameraPosition(distance=position['distance'],
                                           elevation=position['elevation'],
                                           azimuth=position['azimuth'])
                    # HACK: pyqtgraph needs to fix its api so this can be set
                    # with setCameraPosition
                    if 'center' in position:
                        center = position['center']
                        self.opts['center'] = QtGui.QVector3D(center[0],
                                                              center[1],
                                                              center[2])
        self.initializeUI()
        visibility_signal = 'planeVisibilityChanged(int, int)'
        QtCore.QObject.connect(image_display_logic,
                               QtCore.SIGNAL(visibility_signal),
                               self.change_plane_visibility)
        indices_signal = 'planeIndicesChanged(int, int)'
        QtCore.QObject.connect(image_display_logic,
                               QtCore.SIGNAL(indices_signal),
                               self.change_plane_indices)

    def initializeUI(self):
        self.add_grids()
        self.image_planes = [None, None, None]
        self.add_image_planes()
        self.change_plane_visibility(0, False)
        self.tubes_circles = {}
        self.add_tubes()
        self.ultrasound_probe_origin = None
        self.translations = None

    def add_grids(self):
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
            self.addItem(grid)
        axis = gl.GLAxisItem()
        axis_scale = 20
        axis.scale(axis_scale, axis_scale, axis_scale)
        self.addItem(axis)

    def add_image_planes(self):
        """Add the image planes.  The input image is smoothed."""
        image = self.input_image
        sigma = self.config['ParameterGroups'][1]['Parameters'][0]['Value']
        if sigma > 0.0:
            for ii in range(3):
                image = itk.RecursiveGaussianImageFilter(image, sigma, direction=ii)
        image_content = itk.PyBuffer[image.__class__].GetArrayFromImage(image)
        self.image_display_logic.set_number_of_planes(0,
                                                      image_content.shape[2])
        self.image_display_logic.set_number_of_planes(1,
                                                      image_content.shape[1])
        self.image_display_logic.set_number_of_planes(2,
                                                      image_content.shape[0])
        self.image_content = image_content

        for direction in range(3):
            if self.image_planes[direction]:
                self.removeItem(self.image_planes[direction])
            self.image_planes[direction] = self._image_plane(direction)
            if self.image_display_logic.get_plane_visibility(direction):
                self.image_planes[direction].setVisible(True)
            else:
                self.image_planes[direction].setVisible(False)
            self.addItem(self.image_planes[direction])

    def _image_plane(self, direction, index=None):
        """Create an image plane Item from the center plane in the given
        direction for the given ITK Image."""
        if index is None:
            shape = self.image_content.shape
            if direction == 0:
                index = shape[2] / 2
            elif direction == 1:
                index = shape[1] / 2
            elif direction == 2:
                index = shape[0] / 2
        if direction == 0:
            plane = self.image_content[:, :, index]
        elif direction == 1:
            plane = self.image_content[:, index, :]
        else:
            plane = self.image_content[index, :, :]
            plane = plane.transpose()
        texture = pg.makeRGBA(plane)[0]
        image_item = gl.GLImageItem(texture)
        spacing = self.input_image.GetSpacing()
        origin = self.input_image.GetOrigin()
        if direction == 0:
            image_item.scale(spacing[2], spacing[1], 1)
            image_item.rotate(-90, 0, 1, 0)
            image_item.translate(origin[0] + spacing[2] * index,
                                 origin[1],
                                 origin[2])
        elif direction == 1:
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

    def change_plane_indices(self, direction, index):
        if self.image_display_logic.get_plane_visibility(direction):
            if self.image_planes[direction]:
                self.removeItem(self.image_planes[direction])
            self.image_planes[direction] = self._image_plane(direction, index)
            self.addItem(self.image_planes[direction])

    def change_plane_visibility(self, plane, visible):
        if self.image_planes[plane]:
            if visible:
                self.image_planes[plane].setVisible(True)
            else:
                self.image_planes[plane].setVisible(False)

    def add_tubes(self):
        """Add the transformed tubes for the visualization."""
        memory = 4
        target_iterations = tuple(range(self.logic.iteration,
                                  self.logic.iteration - memory,
                                  -1))

        for it in target_iterations:
            if it >= 0 and it <= self.logic.number_of_iterations:
                alpha = np.exp(0.8*(it - self.logic.iteration))
                center_color = (0.8, 0.1, 0.25, alpha)
                center = self.logic.tubes_center
                if self.logic.progression != None:
                    parameters = self.logic.progression[it]['Parameters']
                else:
                    parameters = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                tubes = tubes_from_file(self.logic.subsampled_tubes)
                #tubes_color = (0.2, 0.25, 0.75, alpha)
                have_point_weights_file = False
                if 'TubePointWeightsFile' in self.config and \
                        os.path.exists(self.config['TubePointWeightsFile']):
                    tube_point_weights = self.config['TubePointWeightsFile']
                    statinfo = os.stat(tube_point_weights)
                    if statinfo.st_size != 0:
                        with open(tube_point_weights, 'r') as fp:
                            tube_weights = json.load(fp)['TubePointWeights']
                            tube_weights = np.array(tube_weights)
                            have_point_weights_file = True
                if not have_point_weights_file:
                    tube_weights = 2./(1. + np.exp(-2 * tubes['Radius']))
                tube_weights = tube_weights - tube_weights.min()
                tube_weights = tube_weights / tube_weights.max()
                tubes_colors = matplotlib.cm.PuBuGn(tube_weights)
                tubes_colors[:, 3] = tube_weights**0.5 * alpha
                circles = tubes_as_circles(tubes, point_colors=tubes_colors)
                circles_mesh = gl.GLMeshItem(meshdata=circles,
                                             glOptions='translucent',
                                             smooth=False)
                if self.logic.progression != None:
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
                    self.removeItem(self.tubes_circles[it][1])
                self.addItem(circles_mesh)

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
                    self.removeItem(self.tubes_circles[it][2])
                self.addItem(center_mesh)

                self.tubes_circles[it] = [circles,
                                          circles_mesh,
                                          center_mesh]

        for it, circs in self.tubes_circles.items():
            if not it in target_iterations:
                self.removeItem(circs[1])
                self.removeItem(circs[2])
                self.tubes_circles.pop(it)

    def add_ultrasound_probe_origin(self):
        """Add a sphere indicating the ultrasound probe origin."""
        with open(self.config['UltrasoundProbeGeometryFile'], 'r') as fp:
            line = fp.readline()
            line = line.strip()
            entries = line.split()
            origin = [float(xx) for xx in entries[1:]]
        if self.ultrasound_probe_origin:
            self.removeItem(self.ultrasound_probe_origin)
        sphere = gl.MeshData.sphere(rows=10,
                                    cols=20,
                                    radius=3.0)
        center_mesh = gl.GLMeshItem(meshdata=sphere,
                                    smooth=False,
                                    color=(1.0, 1.0, 0.3, 0.7),
                                    glOptions='translucent')
        center_mesh.translate(origin[0], origin[1], origin[2])
        self.ultrasound_probe_origin = center_mesh
        self.addItem(self.ultrasound_probe_origin)

    def add_translations(self):
        """Add a plot of the translation of the tube centers throughout the
        registration."""
        if self.translations:
            self.removeItem(self.translations)
        translations = np.array(self.logic.progression[:]['Parameters'][:, 3:]) + \
            self.logic.tubes_center
        colors = self.logic.get_iterations_colors()
        translation_points = gl.GLScatterPlotItem(pos=translations,
                                                  color=colors)
        translation_points.setGLOptions('translucent')
        self.translations = translation_points
        self.addItem(self.translations)


class ImageDisplayControlsDock(pyqtgraph.dockarea.Dock):
    """Controls the image planes visualized."""

    def __init__(self, image_display_logic, *args, **kwargs):
        super(ImageDisplayControlsDock, self).__init__(*args, **kwargs)
        self.image_display_logic = image_display_logic
        self.initializeUI()
        planes_changed = image_display_logic.number_of_planes_changed
        planes_changed.connect(self.set_number_of_planes)

    def initializeUI(self):
        self.checkboxes = []
        self.sliders = []

        logic = self.image_display_logic

        x_check = QtGui.QCheckBox("X Plane")
        x_check.setChecked(False)
        self.checkboxes.append(x_check)
        QtCore.QObject.connect(x_check,
                               QtCore.SIGNAL('stateChanged(int)'),
                               logic.set_plane_visibility_x)
        self.addWidget(x_check, row=0, col=0)

        x_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.sliders.append(x_slider)
        QtCore.QObject.connect(x_slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               logic.set_plane_index_x)
        self.addWidget(x_slider, row=0, col=1)

        y_check = QtGui.QCheckBox("Y Plane")
        y_check.setChecked(True)
        self.checkboxes.append(y_check)
        QtCore.QObject.connect(y_check,
                               QtCore.SIGNAL('stateChanged(int)'),
                               logic.set_plane_visibility_y)
        self.addWidget(y_check, row=1, col=0)

        y_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.sliders.append(y_slider)
        QtCore.QObject.connect(y_slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               logic.set_plane_index_y)
        self.addWidget(y_slider, row=1, col=1)

        z_check = QtGui.QCheckBox("Z Plane")
        z_check.setChecked(True)
        self.checkboxes.append(z_check)
        QtCore.QObject.connect(z_check,
                               QtCore.SIGNAL('stateChanged(int)'),
                               logic.set_plane_visibility_z)
        self.addWidget(z_check, row=2, col=0)

        z_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.sliders.append(z_slider)
        QtCore.QObject.connect(z_slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               logic.set_plane_index_z)
        self.addWidget(z_slider, row=2, col=1)

    def set_number_of_planes(self, direction, size):
        if self.sliders[direction].maximum() != size:
            self.sliders[direction].setMaximum(size - 1)


class MetricValueDock(pyqtgraph.dockarea.Dock):
    """Displays the metric value over the iterations."""

    def __init__(self, logic, *args, **kwargs):
        super(MetricValueDock, self).__init__(*args, **kwargs)
        self.logic = logic
        self.initializeUI()

    def initializeUI(self):
        self.metric_values_plot = pg.PlotWidget()
        self.metric_values_plot.setLabel('bottom', text='Iteration')
        self.metric_values_plot.setLabel('left', text='Metric Value Inverse')
        self.metric_values_plot.setLogMode(False, True)
        self.addWidget(self.metric_values_plot)

    def plot_metric_values(self):
        iterations = self.logic.progression[:]['Iteration']
        metric_value_inverse = 1. / \
            self.logic.progression[:]['CostFunctionValue']
        min_metric = np.min(metric_value_inverse)
        self.metric_values_plot.plot({'x': iterations,
                                     'y': metric_value_inverse},
                                     pen={'color': (100, 100, 200, 200),
                                          'width': 2.5},
                                     fillLevel=min_metric,
                                     fillBrush=(255, 255, 255, 30),
                                     antialias=True,
                                     clear=True)
        self.metric_values_plot.plot({'x': [self.logic.iteration],
                                     'y': [metric_value_inverse[self.logic.iteration]]},
                                     symbolBrush='r',
                                     symbol='o',
                                     symbolSize=8.0)


class MetricSpaceLogic(QtCore.QObject):
    """Controls the metric space display logic."""

    _u_direction = 0
    _v_direction = 1

    def __init__(self, logic, parent=None):
        super(MetricSpaceLogic, self).__init__(parent)
        self.logic = logic

    u_direction_changed = QtCore.pyqtSignal(int,
                                            name='uDirectionChanged')

    def get_u_direction(self):
        return self._u_direction

    def set_u_direction(self, direction):
        if self._u_direction != direction:
            self._u_direction = direction
            self.u_direction_changed.emit(direction)

    u_direction = property(get_u_direction,
                           set_u_direction,
                           doc='Metric parameter to display in the u direction.')

    v_direction_changed = QtCore.pyqtSignal(int,
                                            name='vDirectionChanged')

    def get_v_direction(self):
        return self._v_direction

    def set_v_direction(self, direction):
        if self._v_direction != direction:
            self._v_direction = direction
            self.v_direction_changed.emit(direction)

    v_direction = property(get_v_direction,
                           set_v_direction,
                           doc='Metric parameter to display in the v direction.')

    def get_metric_range(self, direction):
        config = self.logic.config
        lower_bound = config['MetricSampler']['LowerBound']
        upper_bound = config['MetricSampler']['UpperBound']
        return lower_bound[direction], upper_bound[direction]

    def set_u_metric_range(self, bounds):
        self.set_metric_range('u', bounds)

    def set_v_metric_range(self, bounds):
        self.set_metric_range('v', bounds)

    def set_metric_range(self, u_or_v, bounds):
        if u_or_v == 'u':
            direction = self.u_direction
        else:
            direction = self.v_direction
        config = self.logic.config
        lower_bound = config['MetricSampler']['LowerBound']
        upper_bound = config['MetricSampler']['UpperBound']
        lower_bound[direction] = bounds[0]
        upper_bound[direction] = bounds[1]
        self.logic.config['MetricSampler']['LowerBound'] = lower_bound
        self.logic.config['MetricSampler']['UpperBound'] = upper_bound

    def set_sampling_size(self, u_or_v, size):
        if u_or_v == 'u':
            direction = self.u_direction
        else:
            direction = self.v_direction
        config = self.logic.config
        config['MetricSampler']['Size'][direction] = size

    def set_u_sampling_size(self, size):
        self.set_sampling_size('u', size)

    def set_v_sampling_size(self, size):
        self.set_sampling_size('v', size)


class MetricSpaceControlsDock(pyqtgraph.dockarea.Dock):
    """Controls what is visualized in the metric space display."""

    _parameters_descriptions = ['X rotation',
                                'Y rotation',
                                'Z rotation',
                                'X translation',
                                'Y translation',
                                'Z translation']

    def __init__(self, logic, metric_space_logic, *args, **kwargs):
        super(MetricSpaceControlsDock, self).__init__(*args, **kwargs)
        self.metric_space_logic = metric_space_logic
        self.logic = logic
        self.initializeUI()

    def initializeUI(self):
        self.u_buttons = []
        self.v_buttons = []

        metric_space_logic = self.metric_space_logic
        logic = self.logic

        u_buttons_layout = QtGui.QGridLayout()
        for ii in range(6):
            button = QtGui.QRadioButton(self._parameters_descriptions[ii])
            QtCore.QObject.connect(button,
                                   QtCore.SIGNAL('clicked(bool)'),
                                   self.set_logic_u_direction)
            u_buttons_layout.addWidget(button, ii % 3, ii / 3)
            self.u_buttons.append(button)
        self.u_buttons[0].setChecked(True)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('uDirectionChanged(int)'),
                               self.set_u_direction)
        u_buttons_widget = QtGui.QWidget()
        u_buttons_widget.setLayout(u_buttons_layout)

        u_controls_layout = QtGui.QVBoxLayout()
        u_controls_layout.addWidget(u_buttons_widget)
        u_controls = QtGui.QGroupBox('U Direction')
        u_controls.setLayout(u_controls_layout)
        self.addWidget(u_controls, row=0, col=0, rowspan=3)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('uDirectionChanged(int)'),
                               self.set_u_direction)

        u_histogram_widget = HistogramWidget()
        self.u_histogram_widget = u_histogram_widget
        self.addWidget(u_histogram_widget, row=0, col=1)
        self.u_image_item = ImageItem()
        u_histogram_widget.setImageItem(self.u_image_item)
        logic.analysis_run.connect(self.update_u_histogram)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('uDirectionChanged(int)'),
                               self.update_u_histogram)
        bounds = self.metric_space_logic.get_metric_range(0)
        u_histogram_widget.setLevels(bounds[0], bounds[1])
        u_histogram_widget.sigLevelChangeFinished.connect(self.set_u_metric_range)

        self.addWidget(QtGui.QLabel('U Samples'), row=1, col=1)
        u_size_spinbox = QtGui.QSpinBox()
        u_size_spinbox.setMinimum(3)
        u_size_spinbox.setSingleStep(2) # keep it odd so samples in the center
        u_size_spinbox.setValue(self.logic.config['MetricSampler']['Size'][0])
        QtCore.QObject.connect(u_size_spinbox,
                               QtCore.SIGNAL('valueChanged(int)'),
                               metric_space_logic.set_u_sampling_size)
        self.addWidget(u_size_spinbox, row=2, col=1)
        self.u_size_spinbox = u_size_spinbox

        v_buttons_layout = QtGui.QGridLayout()
        for ii in range(6):
            button = QtGui.QRadioButton(self._parameters_descriptions[ii])
            QtCore.QObject.connect(button,
                                   QtCore.SIGNAL('clicked(bool)'),
                                   self.set_logic_v_direction)
            v_buttons_layout.addWidget(button, ii % 3, ii / 3)
            self.v_buttons.append(button)
        self.v_buttons[1].setChecked(True)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('vDirectionChanged(int)'),
                               self.set_v_direction)
        v_buttons_widget = QtGui.QWidget()
        v_buttons_widget.setLayout(v_buttons_layout)

        v_controls_layout = QtGui.QVBoxLayout()
        v_controls_layout.addWidget(v_buttons_widget)
        v_controls = QtGui.QGroupBox('V Direction')
        v_controls.setLayout(v_controls_layout)
        self.addWidget(v_controls, row=3, col=0, rowspan=3)

        v_histogram_widget = HistogramWidget()
        self.v_histogram_widget = v_histogram_widget
        self.addWidget(v_histogram_widget, row=3, col=1)
        self.v_image_item = ImageItem()
        v_histogram_widget.setImageItem(self.v_image_item)
        logic.analysis_run.connect(self.update_v_histogram)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('vDirectionChanged(int)'),
                               self.update_v_histogram)
        bounds = self.metric_space_logic.get_metric_range(1)
        v_histogram_widget.setLevels(bounds[0], bounds[1])
        v_histogram_widget.sigLevelChangeFinished.connect(self.set_v_metric_range)

        self.addWidget(QtGui.QLabel('V Samples'), row=4, col=1)
        v_size_spinbox = QtGui.QSpinBox()
        v_size_spinbox.setMinimum(3)
        v_size_spinbox.setSingleStep(2) # keep it odd so samples in the center
        v_size_spinbox.setValue(self.logic.config['MetricSampler']['Size'][1])
        QtCore.QObject.connect(v_size_spinbox,
                               QtCore.SIGNAL('valueChanged(int)'),
                               metric_space_logic.set_v_sampling_size)
        self.addWidget(v_size_spinbox, row=5, col=1)
        self.v_size_spinbox = v_size_spinbox

    def set_u_direction(self, direction):
        """Only set the given direction radio button as checked."""
        for ii in range(6):
            if ii != direction:
                self.u_buttons[ii].setChecked(False)
            else:
                self.metric_space_logic.set_u_direction(ii)
                samples = self.logic.config['MetricSampler']['Size'][ii]
                self.u_size_spinbox.setValue(samples)

    def get_u_direction(self):
        """Get the currently checked u direction."""
        for ii in range(6):
            if self.u_buttons[ii].isChecked():
                return ii
        return 0

    def set_logic_u_direction(self):
        direction = self.get_u_direction()
        self.metric_space_logic.set_u_direction(direction)

    def set_v_direction(self, direction):
        """Only set the given direction radio button as checked."""
        for ii in range(6):
            if ii != direction:
                self.v_buttons[ii].setChecked(False)
            else:
                self.metric_space_logic.set_v_direction(ii)
                samples = self.logic.config['MetricSampler']['Size'][ii]
                self.v_size_spinbox.setValue(samples)

    def get_v_direction(self):
        """Get the currently checked u direction."""
        for ii in range(6):
            if self.v_buttons[ii].isChecked():
                return ii
        return 0

    def set_logic_v_direction(self):
        direction = self.get_v_direction()
        self.metric_space_logic.set_v_direction(direction)

    def update_histogram(self, u_or_v):
        """Update the histogram content and displayed levels for the given
        direction."""
        if u_or_v == 'u':
            direction = self.get_u_direction()
        else:
            direction = self.get_v_direction()

        bounds = self.metric_space_logic.get_metric_range(direction)
        if u_or_v == 'u':
            self.u_histogram_widget.setLevels(bounds[0], bounds[1])
        else:
            self.v_histogram_widget.setLevels(bounds[0], bounds[1])

        progression = self.logic.progression
        if progression is None:
            return
        values = [params[direction] for params in progression[:]['Parameters']]
        values = np.array(values)
        values.shape = (len(values), 1)
        print('parameter values', values)
        if u_or_v == 'u':
            self.u_image_item.setImage(values)
        else:
            self.v_image_item.setImage(values)

    def update_u_histogram(self):
        self.update_histogram('u')

    def update_v_histogram(self):
        self.update_histogram('v')

    def set_u_metric_range(self):
        bounds = self.u_histogram_widget.getLevels()
        self.metric_space_logic.set_u_metric_range(bounds)

    def set_v_metric_range(self):
        bounds = self.v_histogram_widget.getLevels()
        self.metric_space_logic.set_v_metric_range(bounds)


class MetricSpaceDock(pyqtgraph.dockarea.Dock):
    """Displays the metric function over the iterations."""

    _scales = [1000, 1000, 1000, 100, 100, 100]

    def __init__(self, logic, metric_space_logic, *args, **kwargs):
        super(MetricSpaceDock, self).__init__(*args, **kwargs)
        self.logic = logic
        self.metric_space_logic = metric_space_logic
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('uDirectionChanged(int)'),
                               self.plot_metric_values)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('vDirectionChanged(int)'),
                               self.plot_metric_values)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('uDirectionChanged(int)'),
                               self.plot_metric_surface)
        QtCore.QObject.connect(metric_space_logic,
                               QtCore.SIGNAL('vDirectionChanged(int)'),
                               self.plot_metric_surface)
        logic.metric_sampled.connect(self.plot_metric_surface)
        logic.iteration_changed.connect(self.plot_metric_surface)
        self.initializeUI()

    def initializeUI(self):
        self.view_widget = pg.opengl.GLViewWidget()
        grid = pg.opengl.GLGridItem()
        factor = 20
        grid.scale(factor, factor, factor)
        # does not overlap with the axis
        grid.translate(0, 0, -1)
        self.view_widget.addItem(grid)
        self.scatter_plot_item = None
        self.surface_plot_item = None
        self.line_plot_item = None
        axis = gl.GLAxisItem()
        factor = 40
        axis.scale(factor, factor, factor)
        self.view_widget.addItem(axis)
        self.addWidget(self.view_widget)
        self.view_widget.setCameraPosition(distance=500)

    def plot_metric_values(self):
        if not hasattr(self.logic.progression, '__getitem__'):
            return
        metric_value = self.logic.progression[:]['CostFunctionValue']
        parameters = self.logic.progression[:]['Parameters']
        u_direction = self.metric_space_logic.u_direction
        v_direction = self.metric_space_logic.v_direction

        u_value = []
        v_value = []
        for ii in range(len(metric_value)):
            u_value.append(
                parameters[ii][u_direction])
            v_value.append(
                parameters[ii][v_direction])
        u_value = np.array(u_value) * self._scales[u_direction]
        v_value = np.array(v_value) * self._scales[v_direction]
        print('metric_value', metric_value)
        #print('u_value', u_value)
        #print('v_value', v_value)
        pos = np.vstack((u_value, v_value, metric_value)).transpose()
        colors = self.logic.get_iterations_colors()
        if self.scatter_plot_item != None:
            self.scatter_plot_item.setData(pos=pos, color=colors)
        else:
            self.scatter_plot_item = pg.opengl.GLScatterPlotItem(pos=pos,
                                                                 color=colors)
            self.view_widget.addItem(self.scatter_plot_item)
        if self.line_plot_item != None:
            self.line_plot_item.setData(pos=pos)
        else:
            self.line_plot_item = pg.opengl.GLLinePlotItem(pos=pos,
                                                           color=(1.0, 0.8, 0.5,
                                                                  1.0),
                                                           width=2.0,
                                                           antialias=True)
            self.view_widget.addItem(self.line_plot_item)


    def plot_metric_surface(self):
        metric_image = self.logic.metric_image
        statinfo = os.stat(metric_image)
        if statinfo.st_size == 0:
            return
        config = self.logic.config

        TemporaryFile = tempfile.NamedTemporaryFile
        slice_fp = TemporaryFile(suffix='MetricSlice.mha',
                                 delete=False)
        slice_filename = slice_fp.name

        exe = config['Executables']['ExtractMetricSlice']
        metric_image = self.logic.metric_image
        u_direction = self.metric_space_logic.u_direction
        v_direction = self.metric_space_logic.v_direction
        lower_direction = min(u_direction, v_direction)
        higher_direction = max(u_direction, v_direction)
        directions = str(lower_direction)
        directions += ','
        directions += str(higher_direction)

        if self.logic.progression != None:
            iteration = self.logic.iteration
            parameters = self.logic.progression[iteration]['Parameters']
        else:
            parameters = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        lower_bound = config['MetricSampler']['LowerBound']
        upper_bound = config['MetricSampler']['UpperBound']
        size = config['MetricSampler']['Size']
        def get_closest_index(direction, value):
            if value <= lower_bound[direction]:
                return 0
            if value > upper_bound[direction]:
                index = size[direction] - 1
            else:
                step = (upper_bound[direction] - lower_bound[direction]) / \
                    size[direction]
                index = (value - lower_bound[direction]) / step
            return int(index)
        indices = ''
        for ii in range(6):
            if u_direction != ii and v_direction != ii:
                indices += str(get_closest_index(ii, parameters[ii]))
                indices += ','
        indices = indices[:-1]

        subprocess.check_call([exe,
                              '--sliceDirections', directions,
                              '--indices', indices,
                              metric_image, slice_filename])

        slice_image = itk.ImageFileReader(FileName=slice_filename)
        os.remove(slice_filename)

        z = itk.PyBuffer[slicer_image.__class__].GetArrayFromImage(slice_image)
        z = z.transpose()
        if lower_direction == v_direction:
            # correct?
            z = z.transpose()
        u = np.linspace(lower_bound[u_direction], upper_bound[u_direction],
                        z.shape[0])
        u *= self._scales[u_direction]
        v = np.linspace(lower_bound[v_direction], upper_bound[v_direction],
                        z.shape[1])
        v *= self._scales[v_direction]
        max_expected_metric = 100
        z_normalized = z / max_expected_metric
        colors = matplotlib.cm.autumn(z_normalized)
        colors[:, :, 3] = 0.75
        if self.surface_plot_item:
            self.surface_plot_item.setData(x=u, y=v, z=z, colors=colors)
        else:
            self.surface_plot_item = pg.opengl.GLSurfacePlotItem(x=u, y=v, z=z,
                                                                 edgeColor=(0.1, 0.1, 0.1, 0.8),
                                                                 drawEdges=True,
                                                                 shader='shaded',
                                                                 glOptions='translucent',
                                                                 colors=colors)
            self.view_widget.addItem(self.surface_plot_item)


class IterationDockAreaWidget(QtGui.QWidget):
    """Widgets related to iteration analysis. A DockArea containing Dock's that
    give different views on the analysis and a widget to control the current
    iteration."""

    def __init__(self, dock_area, logic):
        super(IterationDockAreaWidget, self).__init__()
        self.dock_area = dock_area
        self.logic = logic
        self.iteration_spinbox = None

        layout = QtGui.QVBoxLayout()
        self.setLayout(layout)

        layout.addWidget(self.dock_area, stretch=1)

        iteration_widget = IterationWidget()
        self.logic.iteration_changed.connect(iteration_widget.set_iteration)
        layout.addWidget(iteration_widget)
        self.logic.number_of_iterations_changed.connect(
            iteration_widget.set_number_of_iterations)
        QtCore.QObject.connect(iteration_widget.spinbox,
                               QtCore.SIGNAL('valueChanged(int)'),
                               logic.set_iteration)
        QtCore.QObject.connect(iteration_widget.slider,
                               QtCore.SIGNAL('valueChanged(int)'),
                               logic.set_iteration)
        self.iteration_widget = iteration_widget


class RegistrationTunerMainWindow(QtGui.QMainWindow):

    def __init__(self, config, logic, parent=None):
        super(RegistrationTunerMainWindow, self).__init__(parent)

        self.config = config
        self.logic = logic
        self.image_display_logic = self.logic.image_display_logic
        metric_space_logic = MetricSpaceLogic(logic)
        self.metric_space_logic = metric_space_logic

        self.initializeUI()
        self.logic.analysis_run.connect(self.image_tubes.add_image_planes)
        self.logic.iteration_changed.connect(self.image_tubes.add_tubes)
        self.logic.iteration_changed.connect(self.metric_value_dock.plot_metric_values)
        self.logic.iteration_changed.connect(self.metric_space_dock.plot_metric_values)
        if 'UltrasoundProbeGeometryFile' in self.config:
            self.logic.analysis_run.connect(self.image_tubes.add_ultrasound_probe_origin)
        self.logic.analysis_run.connect(self.image_tubes.add_translations)

    def initializeUI(self):
        self.resize(1024, 768)
        self.setWindowTitle('Registration Tuner')

        exit_action = QtGui.QAction('&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(QtGui.QApplication.instance().quit)
        self.addAction(exit_action)

        self.dock_area = pg.dockarea.DockArea()
        self.iteration_dock_area = IterationDockAreaWidget(self.dock_area,
                                                           self.logic)

        Dock = pyqtgraph.dockarea.Dock

        run_console = RunConsoleWidget(self.config, self)
        run_console_dock = Dock("Console", size=(400, 300))
        run_console_dock.addWidget(run_console)
        self.dock_area.addDock(run_console_dock, 'right')

        image_display_logic = self.image_display_logic
        dock_title = "Image Display Controls"
        # Note: must be initialized before ImageTubesWidget to get the correct
        # number of planes set.
        image_controls_dock = ImageDisplayControlsDock(image_display_logic,
                                                       dock_title,
                                                       size=(640, 60))
        self.image_tubes = ImageTubesWidget(self.logic,
                                            image_display_logic,
                                            self.config)
        image_tubes_dock = Dock("Image and Tubes", size=(640, 480))
        image_tubes_dock.addWidget(self.image_tubes)
        self.dock_area.addDock(image_tubes_dock, 'left')
        self.dock_area.addDock(image_controls_dock, 'bottom', image_tubes_dock)


        self.metric_value_dock = MetricValueDock(self.logic,
                                                 "Metric Value Inverse",
                                                 size=(500, 300))
        self.dock_area.addDock(self.metric_value_dock, 'left', run_console_dock)

        self.metric_space_dock = MetricSpaceDock(self.logic,
                                                 self.metric_space_logic,
                                                 "Metric Space",
                                                 size=(500, 300))
        self.dock_area.addDock(self.metric_space_dock, 'above',
                               self.metric_value_dock)
        self.metric_space_controls_dock = \
            MetricSpaceControlsDock(self.logic,
                                    self.metric_space_logic,
                                    "Metric Space Controls",
                                    size=(500, 60))
        self.dock_area.addDock(self.metric_space_controls_dock, 'bottom',
                               self.metric_space_dock)

        self.setCentralWidget(self.iteration_dock_area)
        self.show()

    def make_video(self):
        for it in range(self.logic.number_of_iterations + 1):
            filename = os.path.join(self.video_frame_dir,
                                    'frame_{0:04d}.png'.format(it))
            self.logic.iteration = it
            self.update()
            QtCore.QCoreApplication.processEvents()
            pixmap = QtGui.QPixmap.grabWindow(self.winId())
            print('Saving ' + filename)
            pixmap.save(filename, 'png')


class RegistrationTuner(object):
    """Class to drive registration tuning analysis.  This is the class that in
    instantiated and executed."""

    def __init__(self, config,
                 shared_iteration_logic=None,
                 image_display_logic=None, parent=None):
        self.config = config

        self.logic = RegistrationTunerLogic(config,
                                            shared_iteration_logic,
                                            image_display_logic)
        self.logic.initialize()
        self.view = RegistrationTunerMainWindow(config, self.logic)

    def run_analysis(self):
        self.logic.run_analysis()

    def __enter__(self):
        return self

    def close(self):
        self.logic.close()

    def __exit__(self, type, value, traceback):
        self.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=argparse.FileType('r'),
                        help='Configuration file for tuning analysis')
    parser.add_argument('--non-interactive', '-n', action='store_true',
                        help='Run the analysis non-interactively.')
    args = parser.parse_args()
    config_file = args.configuration

    config = json.load(config_file)
    config_file.close()

    app = pg.mkQApp()
    with RegistrationTuner(config) as tuner:
        if(args.non_interactive):
            tuner.run_analysis()
        else:
            sys.exit(app.exec_())
