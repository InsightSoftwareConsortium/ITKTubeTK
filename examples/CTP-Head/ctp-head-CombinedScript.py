#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os 
import sys
import glob
import numpy as np
import itk
from itk import TubeTK as ttk
#from itkwidgets import view


# In[2]:

if len(sys.argv) != 2:
  print("ctp-head-CombinedScript.py <ctp-directory>")
  print("   <ctp-directory> format = C:/Data/unc/HighRes-005-ctp")
  sys.exit()

# NRRD Study Name
print("Directory = ", sys.argv[1])
studyname = sys.argv[1] #'C:/Users/steph/Desktop/Data/unc/HighRes-005-ctp'


# NRRD Files
directory = (studyname + '/')

# Saved NRRD Files 
directory2 = (studyname + '-Reg/')
try:
    os.mkdir(directory2)
except OSError as error:
    print(error)

# Mask Creation and Location
directory3 = (studyname + '-MinMax/')
try:
    os.mkdir(directory3)
except OSError as error:
    print(error)

pic_folder = os.listdir(directory)
pic_folder = [pic_folder for pic_folder in pic_folder if ".nii" in pic_folder]
pic_folder.sort()
print(pic_folder)
num_images = len(pic_folder)

im0Tmp = itk.imread(directory + pic_folder[int(num_images/2)], itk.F)

resample = ttk.ResampleImage.New(Input=im0Tmp,MakeIsotropic=True)
resample.Update()
im0 = resample.GetOutput()
immath = ttk.ImageMath.New(Input=im0)
immath.Blur(1)
im0Blur = immath.GetOutput()

immath.Threshold(150, 800, 1, 0)
immath.Dilate(10, 1, 0)
mask0 = immath.GetOutputUChar()
mask0Tmp = itk.GetArrayViewFromImage(mask0)
mask0Tmp[0:4,:,:] = 0
sizeZ = mask0Tmp.shape[0]
mask0Tmp[sizeZ-4:sizeZ,:,:] = 0   #No need to update mask0 since mask0Tmp is a view of mask0 (shared memory)

itk.imwrite(mask0, directory3 + 'mask.mha', compression=True)
maskObj = itk.ImageMaskSpatialObject[3].New()
maskObj.SetImage(mask0)
maskObj.Update()


# In[3]:


#view(mask0)


# In[ ]:


Dimension = 3
PixelType = itk.ctype('float')
ImageType = itk.Image[PixelType, Dimension]

imdatamax = itk.GetArrayFromImage(im0)
imdatamin = imdatamax
imdatamax2 = imdatamax
imdatamin2 = imdatamax
imdatamax3 = imdatamax
imdatamin3 = imdatamax

imFixedBlur = im0Blur

for imNum in range(num_images):
    imMoving = itk.imread( directory + pic_folder[imNum], itk.F )
    
    immath.SetInput(imMoving)
    immath.Blur(1)
    imMovingBlur = immath.GetOutput()
    
    imreg = ttk.RegisterImages[ImageType].New()
    imreg.SetFixedImage(imFixedBlur)
    imreg.SetMovingImage(imMovingBlur)
    
    imreg.SetRigidMaxIterations(3000)
    imreg.SetRegistration("RIGID")
    imreg.SetExpectedOffsetMagnitude(20)
    imreg.SetExpectedRotationMagnitude(0.3)
    imreg.SetMetric("MEAN_SQUARED_ERROR_METRIC")
    
    imreg.SetFixedImageMaskObject(maskObj)
    #imreg.SetSampleFromOverlap(True)

    imreg.SetReportProgress(True)
    imreg.Update()
    
    tfm = imreg.GetCurrentMatrixTransform()
    #imFixedBlur = imreg.GetFinalMovingImage("LINEAR_INTERPOLATION", -1024)
    imMovingReg = imreg.ResampleImage("LINEAR_INTERPOLATION", imMoving, tfm, -1024)
    
    itk.imwrite( imMovingReg, directory2 + pic_folder[imNum], compression=True )
    
    print(tfm)
    
    imdataTmp = itk.GetArrayFromImage(imMovingReg)
    
    imdatamax = np.maximum(imdatamax,imdataTmp)
    imdatamin = np.minimum(imdatamin,imdataTmp)
    imdataTmp[np.where(imdataTmp==imdatamax)] = 0
    imdataTmp[np.where(imdataTmp==imdatamin)] = 0
    imdatamax2 = np.maximum(imdatamax2,imdataTmp)
    imdatamin2 = np.minimum(imdatamin2,imdataTmp)
    imdataTmp[np.where(imdataTmp==imdatamax)] = 0
    imdataTmp[np.where(imdataTmp==imdatamin)] = 0
    imdatamax3 = np.maximum(imdatamax3,imdataTmp)
    imdatamin3 = np.minimum(imdatamin3,imdataTmp)
    
    #out = itk.GetImageFromArray(imdatamax)
    #out.CopyInformation(im0)
    #itk.imwrite(out, (directory3 + 'max_' + str(imNum) + '.nrrd'))
    
    #out = itk.GetImageFromArray(imdatamax3)
    #out.CopyInformation(im0)
    #itk.imwrite(out, (directory3 + 'max3_' + str(imNum) + '.nrrd'))
    
    percent = (imNum + 1) / num_images * 100
    print(str(round(percent)) + '% : ' + pic_folder[imNum])
    
print('Done')    


# In[ ]:


out = itk.GetImageFromArray(imdatamax3)
out.CopyInformation(im0)
itk.imwrite(out, (directory3 + 'max3.nrrd'), compression=True)

out = itk.GetImageFromArray(imdatamin3)
out.CopyInformation(im0)
itk.imwrite(out, (directory3 + 'min3.nrrd'), compression=True)

out = itk.GetImageFromArray(imdatamax3 - imdatamin3)
out.CopyInformation(im0)
itk.imwrite(out, (directory3 + 'diff3.nrrd'), compression=True)

print('Done3')


# In[ ]:


out = itk.GetImageFromArray(imdatamax)
out.CopyInformation(im0)
itk.imwrite(out, (directory3 + 'max.nrrd'), compression=True)

out = itk.GetImageFromArray(imdatamin)
out.CopyInformation(im0)
itk.imwrite(out, (directory3 + 'min.nrrd'), compression=True)

out = itk.GetImageFromArray(imdatamax - imdatamin)
out.CopyInformation(im0)
itk.imwrite(out, (directory3 + 'diff.nrrd'), compression=True)

print('Done')


# In[ ]:


#out = itk.GetImageFromArray(imdatamax)
#view(out)


# In[ ]:


#!/usr/bin/env python
# coding: utf-8

# This notebook is intended to demonstrate how select registration, segmentation, and image mathematical methods of ITKTubeTK can be combined to perform multi-channel brain extraction (aka. skull stripping for patient data containing multiple MRI sequences).
# 
# There are many other (probably more effective) brain extraction methods available as open-source software such as BET and BET2 in the FSL package (albeit such methods are only for single channel data).   If you need to perform brain extraction for a large collection of scans that do not contain major pathologies, please use one of those packages.   This notebook is meant to show off the capabilities of specific ITKTubeTK methods, not to demonstration how to "solve" brain extraction.

# In[1]:


#import itk
#from itk import TubeTK as ttk

#from itkwidgets import view

#import numpy as np


# In[2]:


ImageType = itk.Image[itk.F, 3]

InputBaseName = studyname + "-MinMax/max3"

filename = InputBaseName + ".nrrd"
im1iso = itk.imread(filename, itk.F)


# In[3]:


N = 8
readerList = ["003", "010", "026", "034", "045", "056", "063", "071"]

imBase = []
imBaseB = []
for i in range(0,N):
    name = "../Data/Normal"+readerList[i]+"-FLASH.mha"
    nameB = "../Data/Normal"+readerList[i]+"-FLASH-Brain.mha"
    imBaseTmp = itk.imread(name, itk.F)
    imBaseBTmp = itk.imread(nameB, itk.F)
    imBase.append(imBaseTmp)
    imBaseB.append(imBaseBTmp)


# In[4]:


#view(im1iso)


# In[5]:


#view(imBase[0])


# In[6]:


imMath = ttk.ImageMath.New(Input=im1iso)
#imMath.Threshold(-4000,-500,1,0)
#headMask = imMath.GetOutput()
imMath.SetInput(im1iso)
#imMath.IntensityWindow(0,1000,1000,0)
#imMath.ReplaceValuesOutsideMaskRange(headMask,-0.5,0.5,-500)
imMath.Blur(1)
imMath.NormalizeMeanStdDev()
imMath.IntensityWindow(-5,5,-500,500)
im1isoBlur = imMath.GetOutput()
#view(im1isoBlur)


# In[7]:


RegisterImagesType = ttk.RegisterImages[ImageType]
regB = []
regBB = []
for i in range(0,N):
    imMath.SetInput(imBase[i])
    imMath.Blur(1)
    imMath.NormalizeMeanStdDev()
    imMath.IntensityWindow(-5,5,-500,500)
    imBaseBlur = imMath.GetOutput()
    
    #regBTo1 = RegisterImagesType.New(FixedImage=im1isoBlur, MovingImage=imBaseBlur)
    regBTo1 = RegisterImagesType.New(FixedImage=imBaseBlur, MovingImage=im1isoBlur)
    regBTo1.SetReportProgress(True)
    
    regBTo1.SetRigidMaxIterations(3000)
    regBTo1.SetAffineMaxIterations(3000)
    
    regBTo1.SetExpectedRotationMagnitude(0.2)
    regBTo1.SetExpectedScaleMagnitude(0.25)
    regBTo1.SetExpectedSkewMagnitude(0.01)
    regBTo1.SetExpectedOffsetMagnitude(40) 

    regBTo1.SetRigidSamplingRatio(0.1)
    regBTo1.SetAffineSamplingRatio(0.1)
    
    regBTo1.SetSampleFromOverlap(True)
    
    regBTo1.SetInitialMethodEnum("INIT_WITH_IMAGE_CENTERS")
    regBTo1.SetRegistration("PIPELINE_AFFINE")
    regBTo1.SetMetric("MATTES_MI_METRIC")
    
    regBTo1.Update()
    
    tfm = regBTo1.GetCurrentMatrixTransform()
    tfmInv = tfm.GetInverseTransform()
    print(tfm)
    
    resm = ttk.ResampleImage.New(Input=imBase[i])
    resm.SetMatchImage(im1iso)
    resm.SetTransform(tfmInv)
    resm.SetLoadTransform(True)
    resm.Update()
    img = resm.GetOutput()
    regB.append( img )

    resm = ttk.ResampleImage.New(Input=imBaseB[i])
    resm.SetMatchImage(im1iso)
    resm.SetTransform(tfmInv)
    resm.SetLoadTransform(True)
    resm.Update()
    img = resm.GetOutput()
    regBB.append( img )


# In[8]:


imMath.SetInput(regB[1])
imMath.AddImages(im1iso,20,1)
img = imMath.GetOutput()
#view( img )


# In[9]:


regBBT = []
for i in range(0,N):
    imMath.SetInput(regBB[i])
    imMath.Threshold(0,1,0,1)
    img = imMath.GetOutput()
    if i==0:
        imMath.SetInput( img )
        imMath.AddImages( img, 1.0/N, 0 )
        sumBBT = imMath.GetOutput()
    else:
        imMath.SetInput( sumBBT )
        imMath.AddImages( img, 1, 1.0/N )
        sumBBT = imMath.GetOutput()
        
#view(sumBBT)


# In[10]:


imMath.SetInput(sumBBT)
imMath.Threshold(0.85,1.1,1,0)
imMath.Dilate(5,1,0)
imMath.Erode(25,1,0)
brainInside = imMath.GetOutput()

imMath.SetInput( sumBBT )
imMath.Threshold(0,0,1,0)
imMath.Erode(1,1,0)
brainOutsideAll = imMath.GetOutput()
imMath.Erode(20,1,0)
imMath.AddImages(brainOutsideAll, -1, 1)
brainOutside = imMath.GetOutput()

imMath.AddImages(brainInside,1,2)
brainCombinedMask = imMath.GetOutputUChar()
brainCombinedMaskF = imMath.GetOutput()


# In[11]:


imMath.SetInput(brainCombinedMaskF)
imMath.AddImages(im1iso, 100, 1)
brainCombinedMaskView = imMath.GetOutput()
#view(brainCombinedMaskView)


# In[12]:


LabelMapType = itk.Image[itk.UC,3]

segmenter = ttk.SegmentConnectedComponentsUsingParzenPDFs[ImageType,LabelMapType].New()
segmenter.SetFeatureImage( im1iso )
segmenter.SetInputLabelMap( brainCombinedMask )
segmenter.SetObjectId( 2 )
segmenter.AddObjectId( 1 )
segmenter.SetVoidId( 0 )
segmenter.SetErodeDilateRadius( 10 )
segmenter.SetHoleFillIterations( 40 )
segmenter.Update()
segmenter.ClassifyImages()
brainCombinedMaskClassified = segmenter.GetOutputLabelMap()


# In[13]:


#view(brainCombinedMaskClassified)


# In[14]:


cast = itk.CastImageFilter[LabelMapType, ImageType].New()
cast.SetInput(brainCombinedMaskClassified)
cast.Update()
brainMaskF = cast.GetOutput()

brainMath = ttk.ImageMath[ImageType,ImageType].New(Input = brainMaskF)
brainMath.Threshold(2,2,1,0)
brainMath.Erode(1,1,0)
brainMaskD = brainMath.GetOutput()
brainMath.SetInput( im1iso )
brainMath.ReplaceValuesOutsideMaskRange( brainMaskD, 1, 1, 0)
brain = brainMath.GetOutput()


# In[15]:


#view(brain)


# In[16]:


writer = itk.ImageFileWriter[ImageType].New(Input = brain)
filename = InputBaseName + "-Brain.nrrd"
writer.SetFileName(filename)
writer.SetUseCompression(True)
writer.Update()


# In[ ]:



#!/usr/bin/env python
# coding: utf-8

# This notebook is intended to demonstrate how select registration, segmentation, and image mathematical methods of ITKTubeTK can be combined to perform multi-channel brain extraction (aka. skull stripping for patient data containing multiple MRI sequences).
# 
# There are many other (probably more effective) brain extraction methods available as open-source software such as BET and BET2 in the FSL package (albeit such methods are only for single channel data).   If you need to perform brain extraction for a large collection of scans that do not contain major pathologies, please use one of those packages.   This notebook is meant to show off the capabilities of specific ITKTubeTK methods, not to demonstration how to "solve" brain extraction.

# In[1]:


#import itk
#from itk import TubeTK as ttk

#from itkwidgets import view

#import numpy as np

# In[2]:

InputBaseDir = studyname 

CTPMaxFilename = InputBaseDir + "-MinMax/max.nrrd"
CTPMinFilename = InputBaseDir + "-MinMax/min.nrrd"
CTPBrainFilename = InputBaseDir + "-MinMax/max3-Brain.nrrd"

imMax = itk.imread(CTPMaxFilename, itk.F)
imMin = itk.imread(CTPMinFilename, itk.F)
imBrain = itk.imread(CTPBrainFilename, itk.F)


# In[3]:


#view(imBrain)


# In[4]:


ImageType = itk.Image[itk.F, 3]

imMath = ttk.ImageMath.New(Input=imBrain)
imMath.Threshold( 0.00001, 2000, 1, 0)
imMath.Erode(10,1,0)
imBrainMaskErode = imMath.GetOutput()

imMath.SetInput(imMax)
imMath.AddImages(imMin,1,-1)
imDiff = imMath.GetOutput()
imMath.ReplaceValuesOutsideMaskRange(imBrain, 0.0001, 2000, 0)
imDiffBrain = imMath.GetOutput()
imMath.ReplaceValuesOutsideMaskRange(imBrainMaskErode, 0.5, 1.5, 0)
imDiffBrainErode = imMath.GetOutput()


# In[5]:


tmpA = itk.GetArrayViewFromImage(imDiffBrain)
tmpAE = itk.GetArrayViewFromImage(imDiffBrainErode)
zMax = tmpA.shape[0]
clip = 0
while((np.amax(tmpA[clip:clip+1,:,:])>1000) | (np.amax(tmpA[clip:clip+1,:,:])==0)):
    clip += 1
if(clip>0):
    tmpA[0:clip,:,:]=0
    tmpAE[0:clip,:,:]=0
clip = 1
while((np.amax(tmpA[zMax-clip:zMax-clip+1,:,:])>1000) | (np.amax(tmpA[zMax-clip:zMax-clip+1,:,:])==0)):
    clip += 1
print(clip, np.amax(tmpA[zMax-clip:zMax-clip+1,:,:]))
clip = clip - 1
if(clip>0):
    tmpA[zMax-clip:zMax,:,:]=0  #Happens to imDiffBrain since this array is a view of an itk image
    tmpAE[zMax-clip:zMax,:,:]=0  #Happens to imDiffBrain since this array is a view of an itk image


# In[6]:


#view(imDiffBrain)


# In[7]:


imMath = ttk.ImageMath[ImageType,ImageType].New()
imMath.SetInput(imDiffBrainErode)
imMath.Blur(1.5)
imBlur = imMath.GetOutput()
imBlurArray = itk.GetArrayViewFromImage(imBlur)

numSeeds = 15
seedCoverage = 20
seedCoord = np.zeros([numSeeds,3])
for i in range(numSeeds):
    seedCoord[i] = np.unravel_index(np.argmax(imBlurArray, axis=None), imBlurArray.shape)
    indx = [int(seedCoord[i][0]),int(seedCoord[i][1]),int(seedCoord[i][2])]
    minX = max(indx[0]-seedCoverage,0)
    maxX = max(indx[0]+seedCoverage,imBlurArray.shape[0])
    minY = max(indx[1]-seedCoverage,0)
    maxY = max(indx[1]+seedCoverage,imBlurArray.shape[1])
    minZ = max(indx[2]-seedCoverage,0)
    maxZ = max(indx[2]+seedCoverage,imBlurArray.shape[2])
    imBlurArray[minX:maxX,minY:maxY,minZ:maxZ]=0
    indx.reverse()
    seedCoord[:][i] = imDiffBrain.TransformIndexToPhysicalPoint(indx)
print(seedCoord)


# In[8]:


# Manually extract a few vessels to form an image-specific training set
vSeg = ttk.SegmentTubes.New(Input=imDiffBrain)
vSeg.SetVerbose(True)
vSeg.SetMinRoundness(0.4)
vSeg.SetMinCurvature(0.002)
vSeg.SetRadiusInObjectSpace( 1 )
for i in range(numSeeds):
    print("**** Processing seed " + str(i) + " : " + str(seedCoord[i]))
    vSeg.ExtractTubeInObjectSpace( seedCoord[i], i )
    
tubeMaskImage = vSeg.GetTubeMaskImage()


# In[9]:


imMath.SetInput(tubeMaskImage)
imMath.AddImages(imDiffBrain, 200, 1)
blendIm = imMath.GetOutput()
#view(blendIm)


# In[10]:


LabelMapType = itk.Image[itk.UC,3]

trMask = ttk.ComputeTrainingMask[ImageType,LabelMapType].New()
trMask.SetInput( tubeMaskImage )
trMask.SetGap( 4 )
trMask.SetObjectWidth( 1 )
trMask.SetNotObjectWidth( 1 )
trMask.Update()
fgMask = trMask.GetOutput()


# In[11]:


#view(fgMask)


# In[12]:


enhancer = ttk.EnhanceTubesUsingDiscriminantAnalysis[ImageType,LabelMapType].New()
enhancer.AddInput( imDiff )
enhancer.SetLabelMap( fgMask )
enhancer.SetRidgeId( 255 )
enhancer.SetBackgroundId( 128 )
enhancer.SetUnknownId( 0 )
enhancer.SetTrainClassifier(True)
enhancer.SetUseIntensityOnly(True)
enhancer.SetScales([0.43,1.29,3.01])
enhancer.Update()
enhancer.ClassifyImages()


# In[13]:


im1vess = itk.SubtractImageFilter( Input1=enhancer.GetClassProbabilityImage(0), Input2=enhancer.GetClassProbabilityImage(1))

imMath.SetInput(imDiffBrain)
imMath.Threshold(0.0001,2000,1,0)
imMath.Erode(2,1,0)
imBrainE = imMath.GetOutput()

imMath.SetInput(im1vess)
imMath.ReplaceValuesOutsideMaskRange(imBrainE, 1, 1, -0.001)
im1vessBrain = imMath.GetOutput()
#view(enhancer.GetClassProbabilityImage(0))
#view(im1vessBrain)


# In[14]:


itk.imwrite( im1vess, InputBaseDir + "-MinMax/diff-VesselEnhanced.nrrd", compression=True)

itk.imwrite( im1vessBrain, InputBaseDir + "-MinMax/Brain-VesselEnhanced.nrrd", compression=True)


# In[ ]:


#!/usr/bin/env python
# coding: utf-8

# This notebook is intended to demonstrate how vessel segmentation methods of ITKTubeTK can be applied to multi-channel MRI (MRA + T1, T2, etc).

# In[1]:


#import itk
#from itk import TubeTK as ttk

#from itkwidgets import view

#import numpy as np


# In[2]:


ImageType = itk.Image[itk.F, 3]

imDir = studyname + "-MinMax/"

im1iso = itk.imread(imDir + "diff3.nrrd")
im1BrainVess = itk.imread(imDir + "Brain-VesselEnhanced.nrrd")


# In[3]:


imMath = ttk.ImageMath.New(im1BrainVess)
imMath.MedianFilter(1)
imMath.Threshold(0.000001, 1, 1, 0)
im1VessMask = imMath.GetOutputShort()

ccSeg = ttk.SegmentConnectedComponents.New(im1VessMask)
ccSeg.SetMinimumVolume(10)
ccSeg.Update()
im1VessMaskCC = ccSeg.GetOutput()


# In[4]:


#view(im1VessMaskCC)


# In[5]:


imMathSS = ttk.ImageMath.New(im1VessMaskCC)
imMathSS.Threshold(0,0,1,0)
im1VessMaskInv = imMathSS.GetOutputFloat()

distFilter = itk.DanielssonDistanceMapImageFilter.New(im1VessMaskInv)
distFilter.Update()
dist = distFilter.GetOutput()

imMath.SetInput(dist)
imMath.Blur(0.4)
tmp = imMath.GetOutput()
imMath.ReplaceValuesOutsideMaskRange(tmp, 0.1, 10, 0)
im1SeedRadius = imMath.GetOutput()

itk.imwrite(im1SeedRadius, imDir+"VesselsSeedRadius.mha")


# In[6]:


#view(im1SeedRadius)


# In[7]:


imMath.SetInput(im1iso)
imMath.ReplaceValuesOutsideMaskRange(im1BrainVess, 0, 1000, 0)
imMath.Blur(0.4)
imMath.IntensityWindow(0.5,1000,0,1000)
im1Input = imMath.GetOutput()

itk.imwrite(im1iso, imDir+"VesselsInput.mha")

#view(im1Input)


# In[8]:


numSeeds = 40

vSeg = ttk.SegmentTubes.New(Input=im1Input)
#vSeg.SetVerbose(True)
vSeg.SetMinCurvature(0)#.0001)
vSeg.SetMinRoundness(0.02)
vSeg.SetMinRidgeness(0.5)
vSeg.SetMinLevelness(0.0)
vSeg.SetRadiusInObjectSpace( 0.8 )
vSeg.SetBorderInIndexSpace(3)
vSeg.SetSeedMask( im1SeedRadius )
vSeg.SetSeedRadiusMask( im1SeedRadius )
vSeg.SetOptimizeRadius(False)
vSeg.SetUseSeedMaskAsProbabilities(True)
vSeg.SetSeedExtractionMinimumProbability(0.4)
#vSeg.SetSeedMaskMaximumNumberOfPoints( numSeeds )
vSeg.ProcessSeeds()


# In[9]:


tubeMaskImage = vSeg.GetTubeMaskImage()
#view(tubeMaskImage)


# In[10]:


SOWriter = itk.SpatialObjectWriter[3].New()
SOWriter.SetInput(vSeg.GetTubeGroup())
SOWriter.SetBinaryPoints(True)
SOWriter.SetFileName( imDir+"Vessels.tre" )
SOWriter.Update()


# In[ ]:









