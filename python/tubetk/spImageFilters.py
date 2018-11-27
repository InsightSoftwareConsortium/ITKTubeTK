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

from ctypes import *
from collections import namedtuple

import Image

import numpy as np

import scipy as sp
import scipy.ndimage
import scipy.stats.mstats

import math

class POINTF( Structure ):
     _fields_ = [( "x", c_float ),
                 ( "y", c_float )]

class POINT( Structure ):
     _fields_ = [( "x", c_int ),
                 ( "y", c_int )]

class POINT3D( Structure ):
     _fields_ = [( "x", c_int ),
                 ( "y", c_int ),
                 ( "z", c_int )]

class LINE( Structure ):
     _fields_ = [( "p1", POINT ),
                 ( "p2", POINT )]

class LINEF( Structure ):
     _fields_ = [( "p1", POINTF ),
                 ( "p2", POINTF )]

Point = namedtuple( 'Point', ['x', 'y'] )
Line = namedtuple( 'Line', ['p1', 'p2'] )

def distanceBtwTwoPoints( point1, point2 ):
    """
        Find distance between two 2D points
    """
    xDistance = point2.x - point1.x
    yDistance = point2.y - point1.y
    xyDistance = math.sqrt( xDistance * xDistance + yDistance * yDistance )
    return xyDistance, xDistance, yDistance

def wienerFilter( array1, maskSize ):
    """
        Apply wiener filter to a given array
    """
    wienerFilteredImage = scipy.signal.wiener( array1, ( maskSize, maskSize ) )

    return wienerFilteredImage


def gaussianFilter( array1, sigma ):
    """
        Apply gaussian filter to a given array
    """
    gauss_lowpass = sp.ndimage.gaussian_filter( array1, sigma )
    gauss_highpass = array1 - gauss_lowpass

    return gauss_lowpass, gauss_highpass


def meanFilter( array1 ):
    """
        Apply mean filtering usign a fixed kernel for a given array1
    """
    kernel = np.array( [[1./9.0, 1./9.0, 1./9.0],
                        [1./9.0, 1./9.0, 1./9.0],
                        [1./9.0, 1./9.0, 1./9.0]] )

    lowpass_3x3 = sp.ndimage.convolve( array1, kernel )


def highPassFilter( array1 ):
    """
        Apply a high pass filter to a given array
    """
    # A slightly "wider", but sill very simple highpass filter
    kernel = np.array( [[-1, -1, -1, -1, -1],
                        [-1,  1,  2,  1, -1],
                        [-1,  2,  4,  2, -1],
                        [-1,  1,  2,  1, -1],
                        [-1, -1, -1, -1, -1]] )

    highpass_5x5 = sp.ndimage.convolve( array1, kernel )


def threshold( inputImage, edge_threshold, lowVal=0, highVal=255 ):

   thresholdedImage1 = sp.stats.mstats.threshold( inputImage, None,
     edge_threshold,  highVal )

   thresholdedImage2 = sp.stats.mstats.threshold( thresholdedImage1,
     edge_threshold+0.00000000001, None, lowVal )

   return thresholdedImage2

def thresholdBinary( x, T, delta=1e-3, minval=1, maxval=0 ):
    x = np.asarray( x ).copy()
    Told = T + 1
    while np.abs( T - Told ) > delta:
        Told = T
        mask = ( x<=T )
        T = ( x[mask].mean() + x[~mask].mean() )/2.

    x[mask] = minval
    x[~mask] = maxval
    return x, T

_4_connected = np.array( [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]], dtype=bool )

_8_connected = np.array( [[1, 1, 1],
                         [1, 1, 1],
                         [1, 1, 1]], dtype=bool )

def removeSmallRegionBlobs( inputImage, regionThreshold ):

    prunedImage1 = removeSmallSizes( inputImage, regionThreshold,
      structure=_8_connected )

    prunedImage2= removeSmallSizes( prunedImage1, regionThreshold,
      structure=_4_connected )

    return prunedImage2

def removeSmallSizes( inputImage, minSize, structure=_4_connected ):

   labelImage, num_objects = sp.ndimage.label( inputImage, structure )

   objectLabels = np.arange( 1, num_objects+1 )

   if( inputImage.max() == 0 ):
    return inputImage

   areas = np.array( sp.ndimage.sum( inputImage / inputImage.max(), labelImage,
                                  objectLabels ) )

   bigObjects = objectLabels[areas >= minSize]

   bigObjectImage = np.zeros( inputImage.shape, dtype=bool )

   for bo in bigObjects:
     bigObjectImage |= labelImage == bo

   return bigObjectImage

def getLargestRegionBlob( inputImage ):

   labelImage, num_objects = sp.ndimage.label( inputImage )
   objectLabels = np.arange( 1, num_objects+1 )

   if( inputImage.max() == 0 ):
    return inputImage

   areas = np.array( sp.ndimage.sum( inputImage / inputImage.max(),
                                     labelImage, objectLabels ) )
   max_size = sp.ndimage.maximum( areas )

   bigObjects = objectLabels[areas == max_size]

   bigObjectImage = np.zeros( inputImage.shape, dtype=bool )
   for bo in bigObjects:
     bigObjectImage |= labelImage == bo

   return bigObjectImage

def findConnectedComponents( inputImage ):
    """
        Finds the connected components in an image
    """
    objects, num_objects = sp.ndimage.label( inputImage )
    object_slices =  sp.ndimage.find_objects( objects )
    return object_slices

def scale( array1 ):
    """
        Scales the values in array to be between 0 and 255
    """
    array2 = normalize( array1 )
    return array2 * 255

def scaleUsingFullWidthHalfMax( array1, stdDevScale=10 ):
    """
        Scales the values in array to be between 0 and 255
    """
    array2 = normalizeUsingFullWidthHalfMax( array1, False )

    array2 = (array2 * stdDevScale) + ( 128 * np.ones( array2.shape ) )
    array2 = array2.clip( min=0, max=255 )

    return array2

def getNMax( arr, n ):
    """

    """
    indices = arr.ravel().argsort()[-n:]
    return indices

def normalize( array1 ):
    """
        Normalizes the values in an array between 0 and 1
    """
    minimum = float( array1.min() )
    array1 = array1-minimum
    maximum = float( array1.max() )
    return array1 / maximum

def normalizeUsingFullWidthHalfMax( array, isArrayOfInts ):
    nBins = 50
    binMin = array.min()
    binMax = array.max()

    loop = 0
    while loop < 4:
        loop += 1
        meanV = 0
        stdDevV = 1

        if binMax - binMin < nBins and isArrayOfInts:
            binMid = (binMax + binMin ) / 2.0
            binStep = ( nBins - 1 ) / 2.0
            binMin = binMid - binStep
            binMax = binMin + (nBins-1)

        bins, binEdges = np.histogram( array, bins=nBins, range=(binMin, binMax) )
        bins = sp.ndimage.gaussian_filter( bins, 0.5 )

        maxBinV = bins.max()
        maxBinList = sp.transpose( ( bins == maxBinV ).nonzero() )[0]
        maxBin = maxBinList[ len(maxBinList)/2 ]

        fwhm = maxBinV / 2

        binFWHMMin = maxBin
        while binFWHMMin > 0 and bins[ binFWHMMin ] >= fwhm:
            binFWHMMin -= 1
        denom = ( bins[ binFWHMMin+1 ] - bins[ binFWHMMin ] )
        if denom != 0:
            binFWHMMin += ( ( fwhm - bins[ binFWHMMin ] ) / denom )

        binFWHMMax = maxBin
        while binFWHMMax < nBins-1 and bins[ binFWHMMax ] >= fwhm:
            binFWHMMax += 1
        denom = ( bins[ binFWHMMax-1 ] - bins[ binFWHMMax ] )
        if denom != 0:
            binFWHMMax -= ( ( fwhm - bins[ binFWHMMax ] ) / denom )

        if binFWHMMax <= binFWHMMin:
            binFWHMMin = maxBin-1
            binFWHMMax = maxBin+1

        minV = ( ( ( binFWHMMin + 0.5 ) / ( nBins - 1.0 ) )
            * ( binMax-binMin ) + binMin )
        maxV = ( ( ( binFWHMMax + 0.5 ) / ( nBins - 1.0 ) )
            * ( binMax-binMin ) + binMin )

        meanV = ( maxV + minV ) / 2.0

        # FWHM to StdDev relationship from
        #   http://mathworld.wolfram.com/GaussianFunction.html
        stdDevV = ( maxV - minV ) / 2.3548

        binMin = meanV - 3 * stdDevV
        binMax = meanV + 3 * stdDevV

    return ( array.astype(float) - meanV ) / stdDevV

def extractSubImage( inputImage, center, mask ):
    """
        Extracts a sub region around a center point in a given image
        Extracted region size is mask by mask
    """
    inputImageHeight = inputImage.shape[0]
    inputImageWidth = inputImage.shape[1]
    corner = np.zeros( ( 4,2 ), dtype =np.int64 )
    # upperleft
    corner[0,0] = center.x - mask.x
    corner[0,1] = center.y - mask.y
    #lowerleft
    corner[1,0] = center.x - mask.x
    corner[1,1] = center.y + mask.y
    #lowerright
    corner[2,0] = center.x + mask.x
    corner[2,1] = center.y + mask.y
    #upperright
    corner[3,0] = center.x + mask.x
    corner[3,1] = center.y - mask.y

    for i in range( 0,4 ):
        for j in range( 0,2 ):
            if ( corner[i][j]< 0 ):
                corner[i][j] =  0
        if( corner[i][1]>=inputImageHeight ):
            corner[i][1] = inputImageHeight-1
        if( corner[i][0]>=inputImageWidth ):
            corner[i][0] = inputImageWidth-1

    return inputImage[corner[0,1]:corner[1,1],corner[1,0]:corner[2,0]], corner


def cropLargestBinaryRegion( inputImage ):

    B = np.argwhere( inputImage )
    ( ystart, xstart ), ( ystop, xstop ) = B.min( 0 ), B.max( 0 ) + 1
    Atrim = inputImage[ystart:ystop, xstart:xstop]

    cropMin = POINT( xstart, ystart )
    cropMax = POINT( xstop, ystop )

    return Atrim, cropMin, cropMax

def cropImageByRegion( inputImage, minColumn, maxColumn, minRow, maxRow ):
    croppedImage = inputImage[minColumn:maxColumn,minRow:maxRow]
    return croppedImage

def floodfill(inputImage, seed, threshMin, threshMax, return_region = False):
    """
    Fill bounded region

    "value" should be a value that otherwise does not occur within the image.
    """

    x, y = seed

    mask = np.zeros( inputImage.shape, dtype=bool )

    if return_region:
        roi = [ (x, y) ]

    edge = [ (x, y) ]
    while edge:
        newedge = []
        for (x, y) in edge:
            for (s, t) in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
                try:
                    item = inputImage[s, t]
                except IndexError:
                    pass
                else:
                    itemB = mask[s, t]
                    if ( item >= threshMin and item <= threshMax and
                         itemB == False ):
                        mask[s, t] = True
                        newedge.append((s, t))
                        if return_region and not (s, t) in roi:
                            roi.append( (s, t) )
        edge = newedge

    if return_region:
        return mask, roi

    return mask

def cropImage( inputImage, roiMinThreshold, roiMaxThreshold ):
    '''
       Hack method to find a dark box within a bright boarder
    '''
    imageT = (inputImage >= roiMinThreshold) & (inputImage <= roiMaxThreshold)
    imageC = getLargestRegionBlob( imageT )
    imageT, cropMin, cropMax = cropLargestBinaryRegion( imageC )

    minRow = cropMin.y
    maxRow = cropMax.y
    minColumn = cropMin.x
    maxColumn = cropMax.x

    croppedImage = inputImage[ minColumn:maxColumn, minRow:maxRow ]

    minCoord = POINT( minColumn, minRow )
    maxCoord = POINT( maxColumn, maxRow )

    return croppedImage, minCoord, maxCoord

def removeNoise( inputImage, noiseSize ):
    #opSize = int( noiseSize / 2 )
    #deNoised1 = sp.ndimage.grey_opening( inputImage,
        #size=(opSize*2+1, opSize*2+1) )
    #deNoised2 = sp.ndimage.grey_closing( deNoised1,
        #size=(opSize*2+1, opSize*2+1) )
    #smoothedImage, tmp = gaussianFilter( deNoised2, noiseSize/3 )
    smoothedImage, tmp = gaussianFilter( inputImage, noiseSize/3 )
    return smoothedImage

def findBackground( inputImage, backgroundScale ):
    background, tmp = gaussianFilter( inputImage, backgroundScale )
    return background

def divideImages( image1, image2 ):
    return np.divide( image1, image2 )

def removeBackground( inputImage, backgroundScale, noiseSize=1 ):
    smoothedImage = removeNoise( inputImage, noiseSize/2 ).astype( float )
    backgroundImage = findBackground( smoothedImage,
        backgroundScale ).astype( float )
    inputImageCleaned = divideImages( smoothedImage, backgroundImage )
    inputImageCleaned = scaleUsingFullWidthHalfMax( inputImageCleaned, 3 )
    return inputImageCleaned
