import __main__
import qt, vtk
import volumesLogic

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

    self.goalSegmentation = None

    self.inputNode1 = None
    self.inputNode2 = None
    self.inputNode3 = None
    self.outputNode = None
    self.goalButtonTexts = ["Jagged", "Intermediate", "Smooth", "Custom"]
    self.goalButtonDefault = 2
    self.goalButtonUserDefined = 3
    self.precomputedErosionRadii = [1, 3, 5] #->> TODO do actual computation
    self.precomputedHoleFillIterations = [2, 4, 6] #->> TODO do actual computation

    self.erosionRadius = 0
    self.holeFillIterations = 0

    self.CLINode = None

    self.editorEffects = ["DefaultTool", "EraseLabel", "Paint", "Threshold"]
    self.editorWidget = None

    self.labelMapNodeSelector = None
    self.voidLabelSpinBox = None
##    self.labelsColorNode = None

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

    #->> TODO could also specify with Qt Designer instead in future (QtUiTools)

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

    #->> TODO for all parameters, provide slots to set them, and use internal values to eventually pass them into the PDF segmenter - for using the interactive PDF segmenter somewhere else, ex in editor module

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

    # LABEL MAP COLLAPSIBLE BUTTON

    labelMapCollapsibleButton = ctk.ctkCollapsibleButton()
    labelMapCollapsibleButton.text = "Label Maps"
    self.layout.addWidget(labelMapCollapsibleButton)

    # Layout within the labelMap collapsible button
    labelMapFormLayout = qt.QFormLayout(labelMapCollapsibleButton)

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

    # Create frame editor widget
    editorFrame = qt.QFrame()
    editorFrame.setLayout(qt.QVBoxLayout())
    palette = editorFrame.palette
    bgColor = 240
    palette.setColor(qt.QPalette.Background, qt.QColor(bgColor, bgColor, bgColor))
    editorFrame.setPalette(palette)
    editorFrame.setAutoFillBackground(True);
    labelMapFormLayout.addRow(editorFrame)
    self.editorFrame = editorFrame

    # initialize editor widget: using parent frame, embedded is true and list of effects
    self.editorWidget = __main__.EditorWidget(parent=self.editorFrame, embedded=True, suppliedEffects=self.editorEffects, showVolumesFrame=False)

    # voidLabel selector
    # The voidLabel selector selects which label corresponds to the void label
    # All other labels in the label map will be extracted and set to object labels
    voidLabelSpinBox = qt.QSpinBox()
    voidLabelSpinBox.objectName = 'voidLabelSpinBox'
    voidLabelSpinBox.toolTip = "Value that represents nothing in the label map.  All other labels represent objects."
    voidLabelSpinBox.setMinimum(0)
    voidLabelSpinBox.setMaximum(255) # temporary value to start
    voidLabelSpinBox.enabled = False
    labelMapFormLayout.addRow("Void Id:", voidLabelSpinBox)
    self.voidLabelSpinBox = voidLabelSpinBox

    #->> TODO: later on, would like a label combo box that shows only those labels
    # that are included in the label map
    # The following code starts in that direction, but does not work

    # BEGIN hacking
   ##  voidLabelSelector = slicer.qMRMLLabelComboBox()
   ##  voidLabelSelector.maximumColorCount = 256 #->> TODO
   ##  labelMapFormLayout.addRow("Void Label:", voidLabelSelector)
   ##  self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)',
   ##                      voidLabelSelector, 'setMRMLScene(vtkMRMLScene*)')
   ##  # create a new vtkMRMLColorTableNode to hold the labels in the label map
   ##  colorLogic = slicer.vtkSlicerColorLogic()
   ##  defaultID = colorLogic.GetDefaultEditorColorNodeID()
   ##  defaultNode = slicer.mrmlScene.GetNodeByID(defaultID)
   ##  if defaultNode:
   ##    # create the node based on the default editor color node
   ##    self.labelsColorNode = slicer.vtkMRMLColorTableNode()
   ##    self.labelsColorNode.Copy(defaultNode)
   ##    # substitute in a new lookup table that we will manipulate
   ##    lookupTable = vtk.vtkLookupTable()
   ##    lookupTable.DeepCopy(defaultNode.GetLookupTable())
   ##    defaultLookupTable = defaultNode.GetLookupTable()
   ##    list = [3,5,7,9]
   ##    lookupTable.SetNumberOfTableValues(len(list))
   ##    for i in range(0, len(list)):
   ##      orig = []
   ##      defaultLookupTable.GetTableValue(list[i], orig)
   ##      lookupTable.SetTableValue(i, defaultNode.GetLookupTable().
   ##    self.labelsColorNode.SetLookupTable(lookupTable)
   ##    # set the new color node to the selector
   ##    # voidLabelSelector.setMRMLColorNode(self.labelsColorNode)
   ##    print "lut:", self.labelsColorNode.GetLookupTable()
   ## self.voidLabelSelector = voidLabelSelector
    # END hacking

    #->> TODO: another alternative is to use an EditColor - but it is heavily coupled
    # to the editor logic, using an editor parameter node that ties it with the editor
    # widget's EditColor
    # Create a frame to give the EditColor a suitable parent
    ## voidLabelFrame = qt.QFrame()
    ## voidLabelFrame.setLayout(qt.QHBoxLayout())
    ## voidLabelFrame.toolTip = "Value that represents nothing in the label map.  All labels not equal to the void id represent objects."
    ## voidLabelSelector = EditorLib.EditColor(parent=voidLabelFrame)
    ## labelMapFormLayout.addRow("Void Id", voidLabelFrame)
    ## self.voidLabelSelector = voidLabelSelector

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
      button.setText(self.goalButtonTexts[i])
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
    reclassifyObjectMaskCheckBox.setChecked(True)
    advancedFormLayout.addRow("Reclassify Object Mask:", reclassifyObjectMaskCheckBox)
    self.reclassifyObjectMaskCheckBox = reclassifyObjectMaskCheckBox

    # reclassifyNotObjectMask check box
    reclassifyNotObjectMaskCheckBox = qt.QCheckBox()
    reclassifyNotObjectMaskCheckBox.objectName = 'reclassifyNotObjectMaskCheckBox'
    reclassifyNotObjectMaskCheckBox.toolTip = "Perform classification on all non-void voxels?"
    reclassifyNotObjectMaskCheckBox.setChecked(True)
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

  def setSliceLabelMaps(self, newLabelMapNode):
    compositeNodes = self.getAllCompositeNodes()
    for node in compositeNodes:
      node.SetReferenceLabelVolumeID(newLabelMapNode.GetID())

  def setLabelMapNode(self, newLabelMapNode):
    """Sets the current node for the 'labelMap' label map
    Connected to signal 'currentNodeChanged()' emitted from the labelMapNodeSelector."""

    if newLabelMapNode:

      if self.labelMapNode == newLabelMapNode:
        return

      # if we don't have a display node (i.e. creating node), add one here
      if (not newLabelMapNode.GetDisplayNode()):
        self.onLabelMapAddedByUser(newLabelMapNode)

      # the editor widget pulls the label map from the red slice's composite node,
      # so set the slice label maps to the new label map node
      #->> problem when toggling between two label maps, colors don't match
      if self.editorWidget:
        self.setSliceLabelMaps(newLabelMapNode)

      # enable the void label spin box only when there is a label map, and set its range
      # to the extent of the label image
      if self.voidLabelSpinBox:
        if (self.voidLabelSpinBox.enabled == (not newLabelMapNode)):
          self.voidLabelSpinBox.enabled = not self.voidLabelSpinBox.enabled
        #->> TODO change limits depending on the values in the image, but must be updated
        # whenever the user adds a new label with the editor
        ## image = newLabelMapNode.GetImageData()

    self.labelMapNode = newLabelMapNode

    if self.labelMapNodeSelector:
      if self.labelMapNodeSelector.currentNode() != newLabelMapNode:
        self.labelMapNodeSelector.setCurrentNode(newLabelMapNode)

  def onLabelMapAddedByUser(self, newLabelMapNode):
    """Creating a new label map volume does not instantiate the image data, so we will
    do that here"""
    logic = volumesLogic.vtkSlicerVolumesLogic()
    logic.SetMRMLScene(slicer.mrmlScene)
    logic.FillLabelVolumeFromTemplate(slicer.mrmlScene, newLabelMapNode, self.inputNode1)

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

      # toggle the status of the label map node selector on whether or not there is an
      # inputNode1 - necessary because we use inputNode1 as the template to create the
      # label map
      if self.labelMapNodeSelector:
        if (self.labelMapNodeSelector.addEnabled == (not newInputNode1)):
          self.labelMapNodeSelector.addEnabled = not self.labelMapNodeSelector.addEnabled

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

    # Calculate the void ID and the object IDs
    # The void ID is provided by the voidLabelSpinBox (for now) and the object IDs
    # are any other labels in the label map that are not the void ID
    voidId = self.voidLabelSpinBox.value
    objectIds = self.getObjectIds(self.labelMapNode, voidId)
    if len(objectIds) == 0:
      print "Error - no valid object Ids"
      return

    parameters = {}
    parameters['inputVolume1'] = self.inputNode1
    parameters['inputVolume2'] = self.inputNode2
    parameters['inputVolume3'] = self.inputNode3
    parameters['voidId'] = voidId
    parameters['objectId'] = objectIds
    parameters['labelmap'] = self.labelMapNode
    parameters['outputVolume'] = self.outputNode
    parameters['erodeRadius'] = self.erosionRadius
    parameters['holeFillIterations'] = self.holeFillIterations
    parameters['fprWeight'] = self.falsePositiveRatioSpinBox.value
    parameters['probSmoothingStdDev'] = self.probabilitySmoothingStdDevSpinBox.value
    parameters['draft'] = self.draftCheckBox.checked
    parameters['reclassifyObjectMask'] = self.reclassifyObjectMaskCheckBox.checked
    parameters['reclassifyNotObjectMask'] = self.reclassifyNotObjectMaskCheckBox.checked

    #->> TODO additional processing here
    #->> cropping
    #->> calculate values for erosion radius and hole fill iterations, based on goal
    # segmentation type and image properties

    # get the pdf segmenter module and run the cli
    tubepdfSegmenter = slicer.modules.tubepdfsegmenter
    self.CLINode = slicer.cli.run(tubepdfSegmenter, self.CLINode, parameters)

    # For a nice display in Slicer, make the output node a label map with the same
    # coloring as the input label map
    if self.outputNode:
      self.outputNode.LabelMapOn()
      if self.labelMapNode:
        #->>TODO Setting the slice label maps is a hack to get a display node
        # We will need to set the slice label maps back to the label map node so that
        # the user can keep editing on it - keeps synchrony with the label map selector
        self.setSliceLabelMaps(self.outputNode)
        labelMapDisplayNode = self.labelMapNode.GetDisplayNode()
        outputDisplayNode = self.outputNode.GetDisplayNode()
        if labelMapDisplayNode and outputDisplayNode:
          outputDisplayNode.SetAndObserveColorNodeID(labelMapDisplayNode.GetColorNodeID())
        self.setSliceLabelMaps(self.labelMapNode)

  def getObjectIds(self, labelMapNode, voidId):
    if not labelMapNode:
      return []
    nonZeroLabels = self.getLabelsFromLabelMap(labelMapNode)
    if len(nonZeroLabels) == 0:
      print "Error - no labels within the label map"
      return []
    if not voidId in nonZeroLabels:
      print "Error - void Id is not represented in the label map"
      print "voidID: ", voidId
      print "label map labels: ", nonZeroLabels
      return []
    nonZeroLabels.remove(voidId)
    return nonZeroLabels

  def getLabelsFromLabelMap(self, labelMapNode):
    if not labelMapNode:
      return
    accum = vtk.vtkImageAccumulate()
    accum.SetInput(labelMapNode.GetImageData())
    accum.UpdateWholeExtent()
    data = accum.GetOutput()
    data.Update()
    numBins = accum.GetComponentExtent()[1]
    nonZeroLabels = []
    for i in range(0, numBins + 1):
      numVoxels = data.GetScalarComponentAsDouble(i,0,0,0)
      if (numVoxels != 0):
        nonZeroLabels.append(i)
    return nonZeroLabels
