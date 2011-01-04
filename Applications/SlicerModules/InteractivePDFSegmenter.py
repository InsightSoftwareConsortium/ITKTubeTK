import __main__
import qt

#
# InteractivePDFSegmenter
#

class InteractivePDFSegmenter:
  def __init__(self, parent):
    parent.title = "Interactive PDF Segmenter (TubeTK)"
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
  def __init__(self, parent=None):

    self.labelMapNode = None
    self.labelMap = None

    self.goalSegmentation = None

    self.inputNode1 = None
    self.inputNode2 = None
    self.inputNode3 = None
    self.outputNode = None
    self.outputProbabilityNode1 = None
    self.outputProbabilityNode2 = None
    self.outputProbabilityNode3 = None
    self.goalButtonTexts = ["J", "S/J", "S", "U"]
    self.goalButtonDefault = 2
    self.goalButtonUserDefined = 3
    self.precomputedErosionRadii = [1, 3, 5] #->> TODO do actual computation
    self.precomputedHoleFillIterations = [2, 4, 6] #->> TODO do actual computation

    self.erosionRadius = 0
    self.holeFillIterations = 0

    self.CLINode = None

    self.editorEffects = ["DefaultTool", "EraseLabel", "Paint", "Threshold"]
    self.editorWidget = None

    if not parent:
      self.parent = qt.QFrame()
      self.parent.setLayout( qt.QVBoxLayout() )
      self.layout = self.parent.layout()
      self.parent.show()
    else:
      self.parent = parent
      self.layout = parent.layout()

  def exit(self):
    if self.editorWidget:
      self.editorWidget.pauseEffect()

  def setup(self):

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
    inputNodeSelector1.nodeTypes = ['vtkMRMLScalarVolumeNode']
    inputNodeSelector1.noneEnabled = False
    inputNodeSelector1.addEnabled = False
    inputNodeSelector1.removeEnabled = False
    inputNodeSelector1.editEnabled = True
    inputNodeSelector1.connect('currentNodeChanged(vtkMRMLNode*)', self.setInputNode1)
    ioFormLayout.addRow("Input Volume 1:", inputNodeSelector1)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        inputNodeSelector1, 'setMRMLScene(vtkMRMLScene*)')
    self.inputNodeSelector1 = inputNodeSelector1

    #->> for all parameters, provide slots to set them, and use internal values to eventually pass them into the PDF segmenter - for using the interactive PDF segmenter somewhere else, ex in editor module

    # inputVolume node selector 2
    inputNodeSelector2 = slicer.qMRMLNodeComboBox()
    inputNodeSelector2.objectName = 'inputNodeSelector2'
    inputNodeSelector2.toolTip = "Select the 2nd input volume to be segmented."
    inputNodeSelector2.nodeTypes = ['vtkMRMLScalarVolumeNode']
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
    inputNodeSelector3.nodeTypes = ['vtkMRMLScalarVolumeNode']
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
    outputNodeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
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
    outputProbabilityNodeSelector1.nodeTypes = ['vtkMRMLScalarVolumeNode']
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
    outputProbabilityNodeSelector2.nodeTypes = ['vtkMRMLScalarVolumeNode']
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
    outputProbabilityNodeSelector3.nodeTypes = ['vtkMRMLScalarVolumeNode']
    outputProbabilityNodeSelector3.noneEnabled = True
    outputProbabilityNodeSelector3.addEnabled = True
    outputProbabilityNodeSelector3.removeEnabled = False
    outputProbabilityNodeSelector3.editEnabled = True
    outputProbabilityNodeSelector3.connect('currentNodeChanged(vtkMRMLNode*)', self.setOutputProbabilityNode3)
    ioFormLayout.addRow("Probability Map for Object 3 (optional):", outputProbabilityNodeSelector3)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        outputProbabilityNodeSelector3, 'setMRMLScene(vtkMRMLScene*)')
    self.outputProbabilityNodeSelector3 = outputProbabilityNodeSelector3

    # LABEL MAP COLLAPSIBLE BUTTON

    labelMapCollapsibleButton = ctk.ctkCollapsibleButton()
    labelMapCollapsibleButton.text = "Label Maps"
    self.layout.addWidget(labelMapCollapsibleButton)

    # Layout within the labelMap collapsible button
    labelMapFormLayout = qt.QFormLayout(labelMapCollapsibleButton)

    # Create frame editor widget
    editorFrame = qt.QFrame()
    editorFrame.setLayout(qt.QHBoxLayout())
    self.layout.addWidget(editorFrame)
    self.editorFrame = editorFrame

    # initialize editor widget: using parent frame, embedded is true and list of effects
    self.editorWidget = __main__.EditorWidget(parent=editorFrame, embedded=True, suppliedEffects=self.editorEffects, showVolumesFrame=False)

    #->> create another selector that picks colors from filled in label maps, to set background color (others will represent objects)

    # labelMap node selector
    labelMapNodeSelector = slicer.qMRMLNodeComboBox()
    labelMapNodeSelector.objectName = 'labelMapNodeSelector'
    labelMapNodeSelector.toolTip = "Select the label map roughly outlining the structure to be segmented and its background."
    labelMapNodeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    labelMapNodeSelector.addAttribute("vtkMRMLScalarVolumeNode", "LabelMap", True);
    labelMapNodeSelector.noneEnabled = False
    labelMapNodeSelector.addEnabled = False
    labelMapNodeSelector.removeEnabled = False
    labelMapNodeSelector.editEnabled = True
    labelMapNodeSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.setLabelMapNode)
    labelMapFormLayout.addRow("Label Map:", labelMapNodeSelector)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
                        labelMapNodeSelector, 'setMRMLScene(vtkMRMLScene*)')
    self.labelMapNodeSelector = labelMapNodeSelector

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
      button.setText(self.goalButtonTexts[i]) #->>TODO replace with icons
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

    # SEGMENTATION BUTTON
    segmentCollapsibleButton = ctk.ctkCollapsibleButton()
    segmentCollapsibleButton.text = "Run Segmentation"
    self.layout.addWidget(segmentCollapsibleButton)

    # Layout within the parameters collapsible button
    segmentFormLayout = qt.QFormLayout(segmentCollapsibleButton)

    # segmentation button
    segmentationButton = qt.QPushButton("Segment")
    segmentationButton.toolTip = "Perform PDF Segmentation."
    segmentFormLayout.addRow(segmentationButton)
    segmentationButton.connect('clicked()', self.onSegmentationButtonClicked)

    # Now that we've created all UI elements, apply the default goal segmentation type
    self.setGoalSegmentationType(self.goalButtonDefault)

  def getAllCompositeNodes(self):
    nodes = []
    count = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLSliceCompositeNode')
    for n in xrange(count):
      nodes.append(slicer.mrmlScene.GetNthNodeByClass(n, 'vtkMRMLSliceCompositeNode'))
    return nodes

  def setLabelMapNode(self, newLabelMapNode):
    """Sets the current node for the 'labelMap' label map
    Connected to signal 'currentNodeChanged()' emitted from the labelMapNodeSelector."""

    newLabelMap = None
    if newLabelMapNode:
      newLabelMap = newLabelMapNode.GetImageData()

      # the editor widget pulls the label map from the red slice's composite node,
      # so set the slice label maps to the new label map node
      #->> problem when adding a new volume, and it automatically switches
      #->> problem when toggling between two label maps, colors don't match
      if self.editorWidget:
        compositeNodes = self.getAllCompositeNodes()
        for node in compositeNodes:
          node.SetReferenceLabelVolumeID(newLabelMapNode.GetID())

    self.labelMapNode = newLabelMapNode
    self.labelMap = newLabelMap

  def setInputNode1(self, newInputNode1):
    """Sets the current node for the 1st input volume
    Connected to signal 'currentNodeChanged()' emitted from the inputNodeSelector1."""

    if newInputNode1:
      # the editor widget pulls the "master node" from the red slice's composite node,
      # so set the slice  background volumes to the new input node 1
      #->> problem when adding a new volume, and it automatically switches
      if self.editorWidget:
        compositeNodes = self.getAllCompositeNodes()
        for node in compositeNodes:
          node.SetReferenceBackgroundVolumeID(newInputNode1.GetID())

    self.inputNode1 = newInputNode1

  def setInputNode2(self, newInputNode2):
    """Sets the current node for the 2nd input volume
    Connected to signal 'currentNodeChanged()' emitted from the inputNodeSelector2."""

    self.inputNode2 = newInputNode2

  def setInputNode3(self, newInputNode3):
    """Sets the current node for the 3rd input volume
    Connected to signal 'currentNodeChanged()' emitted from the inputNodeSelector3."""

    self.inputNode3 = newInputNode3

  def setOutputNode(self, newOutputNode):
    """Sets the current node for the output volume
    Connected to signal 'currentNodeChanged()' emitted from the outputNodeSelector."""

    self.outputNode = newOutputNode

  def setOutputProbabilityNode1(self, newOutputProbabilityNode1):
    """Sets the current node for the 1st output probability node
    Connected to signal 'currentNodeChanged()' emitted from the outputProbabilityNodeSelector1."""

    self.outputProbabilityNode1 = newOutputProbabilityNode1

  def setOutputProbabilityNode2(self, newOutputProbabilityNode2):
    """Sets the current node for the 2nd output probability node
    Connected to signal 'currentNodeChanged()' emitted from the outputProbabilityNodeSelector2."""

    self.outputProbabilityNode2 = newOutputProbabilityNode2

  def setOutputProbabilityNode3(self, newOutputProbabilityNode3):
    """Sets the current node for the 3rd output probability node
    Connected to signal 'currentNodeChanged()' emitted from the outputProbabilityNodeSelector3."""

    self.outputProbabilityNode3 = newOutputProbabilityNode3

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
        self.erosionRadius = newErosionRadius
        self.erosionSpinBox.setValue(newErosionRadius)
    if (goalId < len(self.precomputedHoleFillIterations)):
        newHoleFillIterations = self.precomputedHoleFillIterations[goalId]
        self.holeFillIterations = newHoleFillIterations
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

    #->> voidID pulled from selector widget (to be added)
    #->> objectID all other labels in merge volume that are not the voidID
    parameters['objectId'] = [127, 255]
    parameters['voidId'] = 0

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

    tubepdfSegmenter = slicer.modules.tubepdfsegmenter
    self.CLINode = slicer.cli.run(tubepdfSegmenter, self.CLINode, parameters)

  def getLabelsFromLabelMap(self, labelMapNode):
    if not labelMapNode:
      return
