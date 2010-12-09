
#
# InteractivePDFSegmenter
#

class InteractivePDFSegmenter:
  def __init__(self, parent):
    parent.title = "Interactive PDF Segmenter"
    parent.category = "Segmentation"
    parent.contributor = ""
    parent.helpText = """
    TODO
    """
    parent.acknowledgementText = """
    This work is part of the TubeTK project at Kitware.
    Module implemented by Danielle Pace.  PDF Segmenter implemented by
    Stephen Aylward.
    """
    self.parent = parent

#
# qSlicerPythonModuleExampleWidget
#

class InteractivePDFSegmenterWidget:
  def __init__(self, parent):
    self.parent = parent
    self.layout = parent.layout()

    self.labelMapNode = None
    self.labelMap = None

    self.goalSegmentation = None

    # TODO don't save anything as attributes unless it is needed
    self.inputNode1 = None
    self.inputVolume1 = None
    self.inputNode2 = None
    self.inputVolume2 = None
    self.inputNode3 = None
    self.inputVolume3 = None
    self.outputNode = None
    self.outputVolume = None
    self.outputProbabilityNode1 = None
    self.outputProbabilityVolume1 = None
    self.outputProbabilityNode2 = None
    self.outputProbabilityVolume2 = None
    self.outputProbabilityNode3 = None
    self.outputProbabilityVolume3 = None

    self.goalButtonTexts = ["J", "S/J", "S", "U"]
    self.goalButtonDefault = 2
    self.goalButtonUserDefined = 3
    self.precomputedErosionRadii = [1, 3, 5] # TODO do actual computation
    self.precomputedHoleFillIterations = [2, 4, 6] # TODO do actual computation

    self.erosionRadius = 0
    self.holeFillIterations = 0

    self.CLINode = None

  def setup(self):
    # LABEL MAP COLLAPSIBLE BUTTON
    labelMapCollapsibleButton = ctk.ctkCollapsibleButton()
    labelMapCollapsibleButton.text = "Label Maps"
    self.layout.addWidget(labelMapCollapsibleButton)

    # Layout within the labelMap collapsible button
    labelMapFormLayout = qt.QFormLayout(labelMapCollapsibleButton)

    # TODO integrate Steve Pieper's editor here
    #    insideLabelComboBox = slicer.qMRMLLabelComboBox() # TODO maybe?

    # labelMap node selector
    labelMapNodeSelector = slicer.qMRMLNodeComboBox()
    labelMapNodeSelector.objectName = 'labelMapNodeSelector'
    labelMapNodeSelector.toolTip = "Select the label map roughly outlining the structure to be segmented and its background."
    labelMapNodeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    # TODO select label maps only, using addAttribute (but not wrapped)
    labelMapNodeSelector.noneEnabled = False
    labelMapNodeSelector.addEnabled = False
    labelMapNodeSelector.removeEnabled = False
    labelMapNodeSelector.editEnabled = True
    labelMapNodeSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.setLabelMapNode)
    labelMapFormLayout.addRow("Label Map:", labelMapNodeSelector)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        labelMapNodeSelector, 'setMRMLScene(vtkMRMLScene*)')
    self.labelMapNodeSelector = labelMapNodeSelector

    # inside label spin box # TODO should be list instead
    insideLabelSpinBox = qt.QSpinBox()
    insideLabelSpinBox.objectName = 'insideLabelSpinBox'
    insideLabelSpinBox.toolTip = "Set the label for the structure to be segmented."
    insideLabelSpinBox.setMinimum(0)
    insideLabelSpinBox.setMaximum(1000)
    labelMapFormLayout.addRow("Inside Label:", insideLabelSpinBox)
    self.insideLabelSpinBox = insideLabelSpinBox

    # outside label spin box # TODO should be list instead
    outsideLabelSpinBox = qt.QSpinBox()
    outsideLabelSpinBox.objectName = 'outsideLabelSpinBox'
    outsideLabelSpinBox.toolTip = "Set the label for the structure to be segmented."
    outsideLabelSpinBox.setMinimum(0)
    outsideLabelSpinBox.setMaximum(1000)
    labelMapFormLayout.addRow("Outside Label:", outsideLabelSpinBox)
    self.outsideLabelSpinBox = outsideLabelSpinBox

    # SEGMENTATION PARAMETERS COLLAPSIBLE BUTTON
    segmentationCollapsibleButton = ctk.ctkCollapsibleButton()
    segmentationCollapsibleButton.text = "Segmentation Parameters"
    self.layout.addWidget(segmentationCollapsibleButton)

    # Layout within the parameters collapsible button
    segmentationFormLayout = qt.QFormLayout(segmentationCollapsibleButton)

    # segmentation "goal" buttons
    self.goalButtonList = []
    goalButtonGroup = qt.QButtonGroup()
    goalGroupBox = qt.QGroupBox()
    goalGroupBox.objectName = 'goalGroupBox'
    goalGroupBox.toolTip = "Select what the goal segmentation looks like"
    goalGroupBoxLayout = qt.QHBoxLayout()

    for i in range(0, len(self.goalButtonTexts)):
      button = qt.QToolButton()
      button.setText(self.goalButtonTexts[i]) # TODO replace with icons
      button.setCheckable(True)
      goalButtonGroup.addButton(button, i)
      goalGroupBoxLayout.addWidget(button)
      self.goalButtonList.append(button)

    self.goalButtonList[self.goalButtonDefault].setChecked(True)

    goalButtonGroup.setExclusive(True)
    goalButtonGroup.connect('buttonClicked(int)', self.setGoalSegmentationType)
    goalGroupBox.setLayout(goalGroupBoxLayout)
    goalGroupBox.setFlat(True)

    segmentationFormLayout.addRow("Goal Segmentation:", goalGroupBox)
    self.goalButtonGroup = goalButtonGroup

    # ADVANCED PARAMETERS COLLAPSIBLE BUTTON
    advancedCollapsibleButton = ctk.ctkCollapsibleButton()
    advancedCollapsibleButton.text = "Advanced Parameters"
    self.layout.addWidget(advancedCollapsibleButton)

    # Layout within the parameters collapsible button
    advancedFormLayout = qt.QFormLayout(advancedCollapsibleButton)

    # Erosion radius spin box
    erosionSpinBox = qt.QSpinBox()
    erosionSpinBox.objectName = 'erosionSpinBox'
    erosionSpinBox.toolTip = "Set the erosion radius."
    erosionSpinBox.setMinimum(0)
    erosionSpinBox.connect('valueChanged(int)', self.setErosionRadius)
    advancedFormLayout.addRow("Erosion Radius:", erosionSpinBox)
    self.erosionSpinBox = erosionSpinBox

    # Hole fill iterations spin box
    holeFillSpinBox = qt.QSpinBox()
    holeFillSpinBox.objectName = 'holeFillSpinBox'
    holeFillSpinBox.toolTip = "Set the number of hole filling iterations."
    holeFillSpinBox.setMinimum(0)
    holeFillSpinBox.connect('valueChanged(int)', self.setHoleFillIterations)
    advancedFormLayout.addRow("Hole Fill Iterations:", holeFillSpinBox)
    self.holeFillSpinBox = holeFillSpinBox

    # useTexture check box
    useTextureCheckBox = qt.QCheckBox()
    useTextureCheckBox.objectName = 'useTextureCheckBox'
    useTextureCheckBox.toolTip = "Consider image texture during segmentation?"
    advancedFormLayout.addRow("Use Texture:", useTextureCheckBox)
    self.useTextureCheckBox = useTextureCheckBox

    # falsePositiveRatio spin box
    falsePositiveRatioSpinBox = qt.QDoubleSpinBox()
    falsePositiveRatioSpinBox.objectName = 'falsePositiveRatioSpinBox'
    falsePositiveRatioSpinBox.toolTip = "Relative Cost of False Positive vs. false negative."
    falsePositiveRatioSpinBox.setMinimum(0.0)
    falsePositiveRatioSpinBox.setValue(1.0) # Default
    falsePositiveRatioSpinBox.setSingleStep(0.1)
    advancedFormLayout.addRow("False Positive Ratio:", falsePositiveRatioSpinBox)
    self.falsePositiveRatioSpinBox = falsePositiveRatioSpinBox

    # probabilitySmoothingStandardDeviation spin box
    probabilitySmoothingStdDevSpinBox = qt.QDoubleSpinBox()
    probabilitySmoothingStdDevSpinBox.objectName = 'probabilitySmoothingStdDevSpinBox'
    probabilitySmoothingStdDevSpinBox.toolTip = "Standard deviation of blur applied to probability images prior to computing maximum likelihood of each class at each pixel."
    probabilitySmoothingStdDevSpinBox.setMinimum(0.0)
    probabilitySmoothingStdDevSpinBox.setValue(3.0) # Default
    probabilitySmoothingStdDevSpinBox.setSingleStep(0.1)
    advancedFormLayout.addRow("Probability Smoothing Standard Deviation:", probabilitySmoothingStdDevSpinBox)
    self.probabilitySmoothingStdDevSpinBox = probabilitySmoothingStdDevSpinBox

    # draft check box
    draftCheckBox = qt.QCheckBox()
    draftCheckBox.objectName = 'draftCheckBox'
    draftCheckBox.toolTip = "Generate draft results?"
    advancedFormLayout.addRow("Draft Mode:", draftCheckBox)
    self.draftCheckBox = draftCheckBox

    # reclassifyObjectMask check box
    reclassifyObjectMaskCheckBox = qt.QCheckBox()
    reclassifyObjectMaskCheckBox.objectName = 'reclassifyObjectMaskCheckBox'
    reclassifyObjectMaskCheckBox.toolTip = "Perform classification on voxels within the object mask?"
    advancedFormLayout.addRow("Reclassify Object Mask:", reclassifyObjectMaskCheckBox)
    self.reclassifyObjectMaskCheckBox = reclassifyObjectMaskCheckBox

    # reclassifyNotObjectMask check box
    reclassifyNotObjectMaskCheckBox = qt.QCheckBox()
    reclassifyNotObjectMaskCheckBox.objectName = 'reclassifyNotObjectMaskCheckBox'
    reclassifyNotObjectMaskCheckBox.toolTip = "Perform classification on all non-void voxels?"
    advancedFormLayout.addRow("Reclassify Not Object Mask:", reclassifyNotObjectMaskCheckBox)
    self.reclassifyNotObjectMaskCheckBox = reclassifyNotObjectMaskCheckBox

    # IO COLLAPSIBLE BUTTON
    ioCollapsibleButton = ctk.ctkCollapsibleButton()
    ioCollapsibleButton.text = "IO"
    self.layout.addWidget(ioCollapsibleButton)

    # Layout within the io collapsible button
    ioFormLayout = qt.QFormLayout(ioCollapsibleButton)

    # inputVolume node selector 1
    inputNodeSelector1 = slicer.qMRMLNodeComboBox()
    inputNodeSelector1.objectName = 'inputNodeSelector1'
    inputNodeSelector1.toolTip = "Select the 1st input volume to be segmented."
    inputNodeSelector1.nodeTypes = ['vtkMRMLVolumeNode']
    inputNodeSelector1.noneEnabled = False
    inputNodeSelector1.addEnabled = False
    inputNodeSelector1.removeEnabled = False
    inputNodeSelector1.editEnabled = True
    inputNodeSelector1.connect('currentNodeChanged(vtkMRMLNode*)', self.setInputNode1)
    ioFormLayout.addRow("Input Volume 1:", inputNodeSelector1)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        inputNodeSelector1, 'setMRMLScene(vtkMRMLScene*)')
    self.inputNodeSelector1 = inputNodeSelector1

    # inputVolume node selector 2
    inputNodeSelector2 = slicer.qMRMLNodeComboBox()
    inputNodeSelector2.objectName = 'inputNodeSelector2'
    inputNodeSelector2.toolTip = "Select the 2nd input volume to be segmented."
    inputNodeSelector2.nodeTypes = ['vtkMRMLVolumeNode']
    inputNodeSelector2.noneEnabled = True
    inputNodeSelector2.addEnabled = False
    inputNodeSelector2.removeEnabled = False
    inputNodeSelector2.editEnabled = True
    inputNodeSelector2.connect('currentNodeChanged(vtkMRMLNode*)', self.setInputNode2)
    ioFormLayout.addRow("Input Volume 2 (optional):", inputNodeSelector2)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        inputNodeSelector2, 'setMRMLScene(vtkMRMLScene*)')
    self.inputNodeSelector2 = inputNodeSelector2

    # inputVolume node selector 3
    inputNodeSelector3 = slicer.qMRMLNodeComboBox()
    inputNodeSelector3.objectName = 'inputNodeSelector3'
    inputNodeSelector3.toolTip = "Select the 3rd input volume to be segmented."
    inputNodeSelector3.nodeTypes = ['vtkMRMLVolumeNode']
    inputNodeSelector3.noneEnabled = True
    inputNodeSelector3.addEnabled = False
    inputNodeSelector3.removeEnabled = False
    inputNodeSelector3.editEnabled = True
    inputNodeSelector3.connect('currentNodeChanged(vtkMRMLNode*)', self.setInputNode3)
    ioFormLayout.addRow("Input Volume 3 (optional):", inputNodeSelector3)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        inputNodeSelector3, 'setMRMLScene(vtkMRMLScene*)')
    self.inputNodeSelector3 = inputNodeSelector3

    # outputVolume node selector
    outputNodeSelector = slicer.qMRMLNodeComboBox()
    outputNodeSelector.objectName = 'outputNodeSelector'
    outputNodeSelector.toolTip = "Select the output volume to be segmented."
    outputNodeSelector.nodeTypes = ['vtkMRMLVolumeNode']
    outputNodeSelector.noneEnabled = False
    outputNodeSelector.addEnabled = True
    outputNodeSelector.removeEnabled = False
    outputNodeSelector.editEnabled = True
    outputNodeSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.setOutputNode)
    ioFormLayout.addRow("Output Volume:", outputNodeSelector)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        outputNodeSelector, 'setMRMLScene(vtkMRMLScene*)')
    self.outputNodeSelector = outputNodeSelector

    # outputProbabilityVolume1 node selector
    outputProbabilityNodeSelector1 = slicer.qMRMLNodeComboBox()
    outputProbabilityNodeSelector1.objectName = 'outputProbabilityNodeSelector1'
    outputProbabilityNodeSelector1.toolTip = "Probability-of-being-1st-object estimate for each voxel"
    outputProbabilityNodeSelector1.nodeTypes = ['vtkMRMLVolumeNode']
    outputProbabilityNodeSelector1.noneEnabled = True
    outputProbabilityNodeSelector1.addEnabled = True
    outputProbabilityNodeSelector1.removeEnabled = False
    outputProbabilityNodeSelector1.editEnabled = True
    outputProbabilityNodeSelector1.connect('currentNodeChanged(vtkMRMLNode*)', self.setOutputProbabilityNode1)
    ioFormLayout.addRow("Probability Map for Object 1 (optional):", outputProbabilityNodeSelector1)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        outputProbabilityNodeSelector1, 'setMRMLScene(vtkMRMLScene*)')
    self.outputProbabilityNodeSelector1 = outputProbabilityNodeSelector1

    # outputProbabilityVolume2 node selector
    outputProbabilityNodeSelector2 = slicer.qMRMLNodeComboBox()
    outputProbabilityNodeSelector2.objectName = 'outputProbabilityNodeSelector2'
    outputProbabilityNodeSelector2.toolTip = "Probability-of-being-2st-object estimate for each voxel"
    outputProbabilityNodeSelector2.nodeTypes = ['vtkMRMLVolumeNode']
    outputProbabilityNodeSelector2.noneEnabled = True
    outputProbabilityNodeSelector2.addEnabled = True
    outputProbabilityNodeSelector2.removeEnabled = False
    outputProbabilityNodeSelector2.editEnabled = True
    outputProbabilityNodeSelector2.connect('currentNodeChanged(vtkMRMLNode*)', self.setOutputProbabilityNode2)
    ioFormLayout.addRow("Probability Map for Object 2 (optional):", outputProbabilityNodeSelector2)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        outputProbabilityNodeSelector2, 'setMRMLScene(vtkMRMLScene*)')
    self.outputProbabilityNodeSelector2 = outputProbabilityNodeSelector2

    # outputProbabilityVolume3 node selector
    outputProbabilityNodeSelector3 = slicer.qMRMLNodeComboBox()
    outputProbabilityNodeSelector3.objectName = 'outputProbabilityNodeSelector3'
    outputProbabilityNodeSelector3.toolTip = "Probability-of-being-3st-object estimate for each voxel"
    outputProbabilityNodeSelector3.nodeTypes = ['vtkMRMLVolumeNode']
    outputProbabilityNodeSelector3.noneEnabled = True
    outputProbabilityNodeSelector3.addEnabled = True
    outputProbabilityNodeSelector3.removeEnabled = False
    outputProbabilityNodeSelector3.editEnabled = True
    outputProbabilityNodeSelector3.connect('currentNodeChanged(vtkMRMLNode*)', self.setOutputProbabilityNode3)
    ioFormLayout.addRow("Probability Map for Object 3 (optional):", outputProbabilityNodeSelector3)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        outputProbabilityNodeSelector3, 'setMRMLScene(vtkMRMLScene*)')
    self.outputProbabilityNodeSelector3 = outputProbabilityNodeSelector3

    # SEGMENTATION BUTTON
    segmentationButton = qt.QPushButton("Segment")
    segmentationButton.toolTip = "Perform PDF Segmentation."
    ioFormLayout.addRow(segmentationButton)
    segmentationButton.connect('clicked()', self.onSegmentationButtonClicked)

    # Now that we've created all UI elements, apply the default goal segmentation type
    self.setGoalSegmentationType(self.goalButtonDefault)

    print "DONE"

  def setLabelMapNode(self, newLabelMapNode):
    """Sets the current node for the 'labelMap' label map
    Connected to signal 'currentNodeChanged()' emitted from the labelMapNodeSelector."""

    newLabelMap = None
    if newLabelMapNode:
      newLabelMap = newLabelMapNode.GetImageData()

    self.labelMapNode = newLabelMapNode
    self.labelMap = newLabelMap

  def setInputNode1(self, newInputNode1):
    """Sets the current node for the 1st input volume
    Connected to signal 'currentNodeChanged()' emitted from the inputNodeSelector1."""

    newInputVolume1 = None
    if newInputNode1:
      newInputVolume1 = newInputNode1.GetImageData()

    self.inputNode1 = newInputNode1
    self.inputVolume1 = newInputVolume1

  def setInputNode2(self, newInputNode2):
    """Sets the current node for the 2nd input volume
    Connected to signal 'currentNodeChanged()' emitted from the inputNodeSelector2."""

    newInputVolume2 = None
    if newInputNode2:
      newInputVolume2 = newInputNode2.GetImageData()

    self.inputNode2 = newInputNode2
    self.inputVolume2 = newInputVolume2

  def setInputNode3(self, newInputNode3):
    """Sets the current node for the 3rd input volume
    Connected to signal 'currentNodeChanged()' emitted from the inputNodeSelector3."""

    newInputVolume3 = None
    if newInputNode3:
      newInputVolume3 = newInputNode3.GetImageData()

    self.inputNode3 = newInputNode3
    self.inputVolume3 = newInputVolume3

  def setOutputNode(self, newOutputNode):
    """Sets the current node for the output volume
    Connected to signal 'currentNodeChanged()' emitted from the outputNodeSelector."""

    newOutputVolume = None
    if newOutputNode:
      newOutputVolume = newOutputNode.GetImageData()

    self.outputNode = newOutputNode
    self.outputVolume = newOutputVolume

  def setOutputProbabilityNode1(self, newOutputProbabilityNode1):
    """Sets the current node for the 1st output probability node
    Connected to signal 'currentNodeChanged()' emitted from the outputProbabilityNodeSelector1."""

    newOutputProbabilityVolume1 = None
    if newOutputProbabilityNode1:
      newOutputProbabilityVolume1 = newOutputProbabilityNode1.GetImageData()

    self.outputProbabilityNode1 = newOutputProbabilityNode1
    self.outputProbabilityVolume1 = newOutputProbabilityVolume1

  def setOutputProbabilityNode2(self, newOutputProbabilityNode2):
    """Sets the current node for the 2nd output probability node
    Connected to signal 'currentNodeChanged()' emitted from the outputProbabilityNodeSelector2."""

    newOutputProbabilityVolume2 = None
    if newOutputProbabilityNode2:
      newOutputProbabilityVolume2 = newOutputProbabilityNode2.GetImageData()

    self.outputProbabilityNode2 = newOutputProbabilityNode2
    self.outputProbabilityVolume2 = newOutputProbabilityVolume2

  def setOutputProbabilityNode3(self, newOutputProbabilityNode3):
    """Sets the current node for the 3rd output probability node
    Connected to signal 'currentNodeChanged()' emitted from the outputProbabilityNodeSelector3."""

    newOutputProbabilityVolume3 = None
    if newOutputProbabilityNode3:
      newOutputProbabilityVolume3 = newOutputProbabilityNode3.GetImageData()

    self.outputProbabilityNode3 = newOutputProbabilityNode3
    self.outputProbabilityVolume3 = newOutputProbabilityVolume3

  def setGoalSegmentationType(self, goalId):
    """Sets the goal segmentation 'type': jagged, semi-jagged, smooth or user defined
    Connected to signal 'buttonClicked()' emitted from the goalButtonGroup.
    Computes appropriate values for erosion radius and hole fill iterations"""
    self.computeSegmentationParameters(goalId)
    print self.goalButtonTexts[goalId]

  def computeSegmentationParameters(self, goalId):
    """Actually computes the appropriate values for the erosion radius and hole fill iterations"""
    if (goalId < len(self.precomputedErosionRadii)):
        newErosionRadius = self.precomputedErosionRadii[goalId]
        self.erosionRadius = newErosionRadius # TODO sketchy
        self.erosionSpinBox.setValue(newErosionRadius)
    if (goalId < len(self.precomputedHoleFillIterations)):
        newHoleFillIterations = self.precomputedHoleFillIterations[goalId]
        self.holeFillIterations = newHoleFillIterations # TODO sketchy
        self.holeFillSpinBox.setValue(newHoleFillIterations)

  def setErosionRadius(self, newErosionRadius):
    if (self.erosionRadius != newErosionRadius):
      self.erosionRadius = newErosionRadius
      self.goalButtonList[self.goalButtonUserDefined].click()
    print self.erosionRadius

  def setHoleFillIterations(self, newHoleFillIterations):
    if (self.holeFillIterations != newHoleFillIterations):
      self.holeFillIterations = newHoleFillIterations
      self.goalButtonList[self.goalButtonUserDefined].click()
    print self.holeFillIterations

  def onSegmentationButtonClicked(self):
    parameters = {}
    parameters['inputVolume1'] = self.inputNode1
    parameters['inputVolume2'] = self.inputNode2
    parameters['inputVolume3'] = self.inputNode3
    parameters['objectId'] = self.insideLabelSpinBox.value
    parameters['voidId'] = self.outsideLabelSpinBox.value
    parameters['labelmap'] = self.labelMapNode
    parameters['outputVolume'] = self.outputNode
    parameters['useTexture'] = self.useTextureCheckBox.checked
    parameters['erodeRadius'] = self.erosionRadius
    parameters['holeFillIterations'] = self.holeFillIterations
    parameters['fprWeight'] = self.falsePositiveRatioSpinBox.value
    parameters['probSmoothingStdDev'] = self.probabilitySmoothingStdDevSpinBox.value
    parameters['draft'] = self.draftCheckBox.checked
    parameters['reclassifyObjectMask'] = self.reclassifyObjectMaskCheckBox.checked
    parameters['reclassifyNotObjectMask'] = self.reclassifyNotObjectMaskCheckBox.checked
    parameters['probabilityVolume0'] = self.outputProbabilityNode1
    parameters['probabilityVolume1'] = self.outputProbabilityNode2
    parameters['probabilityVolume2'] = self.outputProbabilityNode3

    pdfSegmenter = slicer.modules.pdfsegmenter
    self.CLINode = slicer.cli.run(pdfSegmenter, self.CLINode, parameters)

    print "SEGMENTED"
