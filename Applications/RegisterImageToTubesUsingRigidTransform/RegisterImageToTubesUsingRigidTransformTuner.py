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
from pyqtgraph.Qt import QtCore, QtGui


class RegistrationTuner(QtGui.QMainWindow):

    def __init__(self, config, config_filename):
        super(RegistrationTuner, self).__init__()

        self.config = config
        self.config_filename = config_filename
        self.initializeUI()

    def initializeUI(self):
        self.resize(1024, 768)
        self.setWindowTitle('Registration Tuner')

        exit_action = QtGui.QAction('&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(QtGui.QApplication.instance().quit)
        self.addAction(exit_action)

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

        self.setCentralWidget(run_console)

        self.show()

    def run_analysis(self):
        io_params = self.config['ParameterGroups'][0]['Parameters']
        print(io_params[2]['Value'])
        subprocess.check_call([config['Executables']['Analysis'],
                              '--parameterstorestore', self.config_filename,
                              io_params[0]['Value'],
                              io_params[1]['Value'],
                              io_params[2]['Value']])

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
