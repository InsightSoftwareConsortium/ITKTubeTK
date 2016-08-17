#!/usr/bin/python

###########################################################################
# reconstructSlabs.py:
#
# Retrieve all processed slabs for the specified animal and save them as a
# single .mha meta-image.
#
###########################################################################

import sys
## Append ITK libs
sys.path.append('/home/lucas/Projects/ITK-Release/Wrapping/Generators/Python')
sys.path.append('/home/lucas/Projects/ITK-Release/Modules/ThirdParty/VNL/src/vxl/lib')
import itk

def reconstructSlabs( animalName, directory ):
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
  nameGenerator.SetSeriesFormat( directory + "%d_" + animalName )

  #seriesReader.SetImageIO( itk.PNGImageIO.New() )
  seriesReader.SetFileNames( nameGenerator.GetFileNames() )
  writer.SetFileName( directory + animalName + ".mha" )
  writer.SetInput( seriesReader.GetOutput() )
  writer.Update()



hardDrive_root = "/media/lucas/krs0014/"

outputDir = hardDrive_root + "SegmentVesselsUsingNeuralNetworks/output/"
outputAnimal = "pp07_A34_left.png" #WARNING Hardcoded

reconstructSlabs( outputAnimal, outputDir )
