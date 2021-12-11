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

    #include_dirs=[np.get_include()],
setup(
    name='itk-tubetk',
    version='1.0.0',
    author='Stephen R. Aylward',
    author_email='stephen.aylward@kitware.com',
    packages=['itk'],
    package_dir={'itk': 'itk'},
    download_url=r'https://github.com/InsightSoftwareConsortium/ITKTubeTK',
    description=r'An open-source toolkit, led by Kitware, Inc., for the segmentation, registration, and analysis of tubes and surfaces in images.',
    long_description='TubeTK is an open-source toolkit for the segmentation, registration, and analysis of tubes and surfaces in images.',
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
        "Operating System :: Android",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS"
        ],
    license='Apache',
    keywords='ITK InsightToolkit',
    url=r'https://itk.org/',
    install_requires=[
        r'numpy',
        r'itk>=5.3rc2.post1',
        r'itk-minimalpathextraction>=1.2.0'
    ]
    )
