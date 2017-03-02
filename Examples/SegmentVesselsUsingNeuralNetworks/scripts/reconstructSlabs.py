#!/usr/bin/python

###########################################################################
# reconstructSlabs.py:
#
# Retrieve all processed slabs for the specified animal and save them as a
# single .mha meta-image.
#
###########################################################################

import json
import os
import sys

## Append ITK libs
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Modules/ThirdParty/VNL/src/vxl/lib'))
import itk


def reconstructSlabs(animalName, directory):

    PixelType = itk.UC
    Dimension = 3
    ImageType=itk.Image[PixelType,Dimension]

    SeriesReaderType = itk.ImageSeriesReader[ImageType]
    WriterType = itk.ImageFileWriter[ImageType]

    seriesReader = SeriesReaderType.New()
    writer = WriterType.New()

    NameGeneratorType = itk.NumericSeriesFileNames
    nameGenerator = NameGeneratorType.New()

    nameGenerator.SetStartIndex( 0 )
    nameGenerator.SetEndIndex( 9 )
    nameGenerator.SetIncrementIndex( 1 )
    nameGenerator.SetSeriesFormat(str(os.path.join(directory, "%d_" + animalName)))

    #  seriesReader.SetImageIO( itk.PNGImageIO.New() )
    seriesReader.SetFileNames( nameGenerator.GetFileNames() )
    writer.SetFileName(str(os.path.join(directory, animalName + ".mha")))
    writer.SetInput(seriesReader.GetOutput())
    writer.Update()

# Define paths
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']

outputDir = os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks/output/")
outputAnimal = "pp07_A36_left.png"  # WARNING Hardcoded

reconstructSlabs(outputAnimal, outputDir)
