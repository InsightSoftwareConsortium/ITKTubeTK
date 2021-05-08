#!/usr/bin/python

###########################################################################
# UnShrink.py :
#
# Map all slabs pixel back into 3D image space
#
###########################################################################

import json
import os
import sys

## Append ITK lib dirs to PYTHONPATH and import itk
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Wrapping/Generators/Python'))
sys.path.append(os.path.join(os.environ['TubeTK_BUILD_DIR'], 'ITK-build/Modules/ThirdParty/VNL/src/vxl/lib'))
import itk

# Define paths
script_params = json.load(open('params.json'))
caffe_root = script_params['CAFFE_SRC_ROOT']
hardDrive_root = script_params['CNN_DATA_ROOT']

hardDrive_vseg_root = os.path.join(hardDrive_root, "SegmentVesselsUsingNeuralNetworks")

# Animal name. WARNING: This script will look for the animal point map generated
#  by ConvertData.py and saved in hardDrive_root + "SegmentVesselsUsingNeuralNetworks/testing/images/
#  The animal name has to come from the testing set, otherwise re-run the ConvertData.py
#  script with the appropriate location
animal = "pp07_A36_left"

pointImage = str(os.path.join(hardDrive_vseg_root, "testing/points", animal + "_points.mha"))
vesselImage = str(os.path.join(hardDrive_vseg_root, "output", animal + ".png.mha"))

########
# Main #
########
PixelType = itk.F
Dimension = 3

# Read segmented vessels slabs
ImageType=itk.Image[PixelType, Dimension]
ReaderType = itk.ImageFileReader[ImageType]
reader = ReaderType.New()
reader.SetFileName(vesselImage)
reader.Update()
vesselImg = reader.GetOutput()
vesselBuf = itk.PyBuffer[ImageType].GetArrayFromImage(vesselImg)

# Read point map image
VectorImageType=itk.VectorImage[PixelType, Dimension]
VectorReaderType = itk.ImageFileReader[VectorImageType]
vectorReader = VectorReaderType.New()
vectorReader.SetFileName(pointImage)
vectorReader.Update()
img = vectorReader.GetOutput()
buf = itk.PyBuffer[VectorImageType].GetArrayFromImage(img)

# Create output 3D image
RegionType=itk.ImageRegion[Dimension]
region = itk.ImageRegion[Dimension]([0,0,0], [512,512,393])
img3d = ImageType.New()
img3d.CopyInformation(img)
img3d.SetRegions(region)
img3d.Allocate()
img3d.SetSpacing([0.05,0.05,0.05])  # WARNING: Spacing and Origin should be read from the input .mhd ultrasound scan
img3d.SetOrigin([0,0,0])

# Iterate through slabs
for slabIdx in range(0,vesselBuf.shape[0]):
    print "Reconstructing slab",slabIdx,"/",vesselBuf.shape[0]-1,"..."
    for i in range(vesselBuf.shape[1]):
        for j in range(vesselBuf.shape[2]):
            if vesselBuf[slabIdx,i,j] > 249: # WARNING: Threshold value can be changed
                # Get current bright pixel 3D position
                x=float(buf[slabIdx,i,j,0])
                y=float(buf[slabIdx,i,j,1])
                z=float(buf[slabIdx,i,j,2])
                # Paint pixel back into the 3D image
                index = img3d.TransformPhysicalPointToIndex([x,y,z])
                img3d.SetPixel(index, 255)

# Blur image to fill gaps between painted pixels.
BlurFilterType=itk.DiscreteGaussianImageFilter[ImageType,ImageType]
blurFilter = BlurFilterType.New()
blurFilter.SetInput(img3d)
blurFilter.SetVariance(1)
blurFilter.SetMaximumKernelWidth(1)

print("Writing output ..")
WriterType = itk.ImageFileWriter[ImageType]
writer = WriterType.New()
writer.UseCompressionOn()
writer.SetFileName(str(os.path.join(hardDrive_vseg_root, "output/3d.mha")))
writer.SetInput( blurFilter.GetOutput() )
writer.Update()
