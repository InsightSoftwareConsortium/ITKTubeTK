# -*- coding: utf-8 -*-
from __future__ import print_function
from os import sys

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

from pathlib import Path
this_directory = Path(__file__).parent
setup_readme_text = (this_directory / "setup-readme.md").read_text()

    #include_dirs=[np.get_include()],
setup(
    name='itk-tubetk',
    version='1.2',
    author='Stephen R. Aylward',
    author_email='stephen.aylward@kitware.com',
    packages=['itk'],
    package_dir={'itk': 'itk'},
    download_url=r'https://github.com/InsightSoftwareConsortium/ITKTubeTK',
    description=r'An open-source toolkit, led by Kitware, Inc., for the segmentation, registration, and analysis of tubes and surfaces in images.',
    long_description=setup_readme_text,
    long_description_content_type='text/markdown',
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS"
        ],
    license='Apache',
    keywords='ITK InsightToolkit Tubes Vessels Nerves Ultrasound MRI CT Medical',
    url=r'https://github.com/InsightSoftwareConsortium/ITKTubeTK/',
    project_urls={
        'Dashboard': 'https://open.cdash.org/index.php?project=TubeTK',
        'Issue Tracker': 'https://github.com/InsightSoftwareConsortium/ITKTubeTK/issues',
        'Testing Data': 'https://data.kitware.com/#collection/5888b7d38d777f4f3f3085a8/folder/58a3abf08d777f0721a65b16',
        'ITK': 'https://itk.org',
        },
    install_requires=[
        r'numpy',
        r'itk>=5.3rc04',
        r'itk-minimalpathextraction>=1.2.0'
    ]
    )
