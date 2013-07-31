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
#!/usr/bin/env python

import Image

import cv2
import cv2.cv as cv
import numpy as np

def cv2array(im):
  '''
    Copied from http://opencv.willowgarage.com/wiki/PythonInterface
  '''
  depth2dtype = {
        cv.IPL_DEPTH_8U: 'uint8',
        cv.IPL_DEPTH_8S: 'int8',
        cv.IPL_DEPTH_16U: 'uint16',
        cv.IPL_DEPTH_16S: 'int16',
        cv.IPL_DEPTH_32S: 'int32',
        cv.IPL_DEPTH_32F: 'float32',
        cv.IPL_DEPTH_64F: 'float64',
    }

  arrdtype=im.depth
  a = np.fromstring(
         im.tostring(),
         dtype=depth2dtype[im.depth],
         count=im.width*im.height*im.nChannels)
  if im.nChannels is not 1:
    a.shape = (im.height,im.width,im.nChannels)
  else:
    a.shape = (im.height,im.width)
  return a

def array2cv(a):
  '''
    Copied from http://opencv.willowgarage.com/wiki/PythonInterface
  '''
  dtype2depth = {
        'bool':   cv.IPL_DEPTH_8U,
        'uint8':   cv.IPL_DEPTH_8U,
        'int8':    cv.IPL_DEPTH_8S,
        'uint16':  cv.IPL_DEPTH_16U,
        'int16':   cv.IPL_DEPTH_16S,
        'int32':   cv.IPL_DEPTH_32S,
        'float32': cv.IPL_DEPTH_32F,
        'float64': cv.IPL_DEPTH_64F,
    }
  try:
    nChannels = a.shape[2]
  except:
    nChannels = 1
  cv_im = cv.CreateImageHeader((a.shape[1],a.shape[0]),
          dtype2depth[str(a.dtype)],
          nChannels)
  cv.SetData(cv_im, a.tostring(),
             a.dtype.itemsize*nChannels*a.shape[1])
  return cv_im

def cv2Image(cv_im):
  pi = Image.fromstring( "L", cv.GetSize(cv_im), cv_im.tostring() )
  return pi

def Image2cv(pi):
  cv_im = cv.CreateImageHeader( pi.size, cv.IPL_DEPTH_8U, 1 )
  cv.SetData( cv_im, pi.tostring() )
  return cv_im


def gaussianFilter( inputImage, filterSize ):
    """
        Apply Gaussian filter of OpenCV to a given array
    """
    # Convert from numpy array  to CvMat
    outputImage = cv.CreateImage( cv.GetSize( inputImage ), cv.IPL_DEPTH_32F, 1 )

    cv.Smooth( inputImage, outputImage, cv.CV_GAUSSIAN, filterSize, filterSize )

    return outputImage


def adaptiveThresholding( inputImage, neighborhoodWidth=71, offsetFromMean=15 ):
    """
        Apply adaptive thresholding to a given image.  Uses a
         neighborhoodWidth x neighborhoodWidth kernel.   Threshold is set at
         mean intensity within kernel + offsetFromMean.
    """
    outputImage = cv.CreateImage( cv.GetSize( inputImage ), cv.IPL_DEPTH_8U, 1 )

    cv.AdaptiveThreshold( inputImage, outputImage, 255, cv.CV_THRESH_BINARY,
                          cv.CV_ADAPTIVE_THRESH_MEAN_C, neighborhoodWidth,
                          offsetFromMean )

    return outputImage
