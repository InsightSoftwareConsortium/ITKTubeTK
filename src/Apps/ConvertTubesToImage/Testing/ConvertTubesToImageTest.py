import os
import sys
from ctk_cli import CLIArgumentParser

# Append ITK libs
sys.path.append(os.path.join(os.environ['ITK_BUILD_DIR'],
                             'Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['ITK_BUILD_DIR'], 'lib')

# Append TubeTK libs :
sys.path.append(os.environ['TUBETK_BUILD_DIR'], 'TubeTK-build/lib/TubeTK')

import itk
from itk import TubeTKITK as itktube

def run(args):

  PixelType = itk.UC
  Dimension = 3

  # Read tre file
  TubeFileReaderType = itk.SpatialObjectReader[Dimension]

  tubeFileReader = TubeFileReaderType.New()
  tubeFileReader.SetFileName(args.inputTubeFile)
  tubeFileReader.Update()

  # Read template image
  TemplateImageType = itk.Image[PixelType, Dimension]
  TemplateImageReaderType = itk.ImageFileReader[TemplateImageType]

  templateImageReader = TemplateImageReaderType.New()
  templateImageReader.SetFileName(args.inputTemplateImage)
  templateImageReader.Update()

  # call ConvertTubesToImage
  TubesToImageFilterType = itktube.ConvertTubesToImage[Dimension, PixelType]

  tubesToImageFilter = TubesToImageFilterType.New()
  tubesToImageFilter.SetUseRadius(args.useRadii)
  tubesToImageFilter.SetTemplateImage(templateImageReader.GetOutput())
  tubesToImageFilter.SetInput(tubeFileReader.GetOutput())

  # write output image
  TubeImageWriterType = itk.ImageFileWriter[TemplateImageType]
  tubeImageWriter = TubeImageWriterType.New()
  tubeImageWriter.SetInput(tubesToImageFilter.GetOutput())
  tubeImageWriter.SetFileName(args.outputImageFile)
  tubeImageWriter.Update()

if __name__ == "__main__":
    run(CLIArgumentParser().parse_args())
