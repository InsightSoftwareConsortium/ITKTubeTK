import os
from __main__ import vtk, qt, ctk, slicer
import EditorLib
from EditorLib.EditOptions import HelpButton
from EditorLib.EditOptions import EditOptions
from EditorLib import EditUtil
from EditorLib import LabelEffect


class InteractiveConnectedComponentsUsingParzenPDFsOptions(EditorLib.LabelEffectOptions):
  """ Editor Effect gui
  """

  def __init__(self, parent=0):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).__init__(parent)

    self.attributes = ('PaintTool')
    self.displayName = 'Interactive Connected Components using Parzen PDFs'
    self.undoRedo = EditorLib.EditUtil.UndoRedo()

    # Create the normal PDF segmenter node cli if it doesn't exists yet.
    # This is because we want the normal module's cli to be selected
    # when opening the cli module.
    module = slicer.modules.segmentconnectedcomponentsusingparzenpdfs
    self.logic.getCLINode(module, module.title)

  def __del__(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).__del__()

  def create(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).create()

    ioCollapsibleButton = ctk.ctkCollapsibleGroupBox()
    ioCollapsibleButton.title = "IO"
    ioCollapsibleButton.collapsed = 0
    self.frame.layout().addWidget(ioCollapsibleButton)

    # Layout within the io collapsible button
    ioFormLayout = qt.QFormLayout(ioCollapsibleButton)
    self.additionalInputNodeSelectors = []
    for i in range(0,2):
      self.additionalInputNodeSelectors.append(self.addInputNodeSelector(i, ioFormLayout))
    self.additionalInputNodeSelectors[0].toolTip = "Select the 1st additional input volume to be segmented"
    self.additionalInputNodeSelectors[1].toolTip = "Select the 2nd additional input volume to be segmented"

    # Objects
    objectCollapsibleGroupBox = ctk.ctkCollapsibleGroupBox()
    objectCollapsibleGroupBox.title = "Objects"
    self.frame.layout().addWidget(objectCollapsibleGroupBox)

    # Layout within the io collapsible button
    objectFormLayout = qt.QFormLayout(objectCollapsibleGroupBox)
    foregroundLayout = qt.QHBoxLayout()
    foregroundLabel = slicer.qMRMLLabelComboBox()
    foregroundLabel.objectName = 'Foreground label'
    foregroundLabel.setMRMLScene(slicer.mrmlScene)
    foregroundLabel.setMRMLColorNode(self.editUtil.getColorNode())
    foregroundLabel.labelValueVisible = True
    foregroundLabel.currentColor = 1
    self.foregroundLabel = foregroundLabel
    self.connections.append( (self.foregroundLabel, 'currentColorChanged(int)', self.updateMRMLFromGUI ) )
    foregroundWeightSpinBox = qt.QDoubleSpinBox(foregroundLabel)
    self.foregroundWeightSpinBox = foregroundWeightSpinBox
    foregroundWeightSpinBox.setRange(0.0, 1.0)
    foregroundWeightSpinBox.setSingleStep(0.1)
    foregroundWeightSpinBox.value = 1.0
    foregroundPopup = ctk.ctkPopupWidget( foregroundWeightSpinBox )
    foregroundPopupLayout = qt.QHBoxLayout( foregroundPopup )
    foregroundPopupSlider = ctk.ctkDoubleSlider( foregroundPopup )
    self.foregroundPopupSlider = foregroundPopupSlider
    foregroundPopupSlider.maximum = 1.0
    foregroundPopupSlider.minimum = 0.0
    foregroundPopupSlider.singleStep = 0.1
    foregroundPopupSlider.connect('valueChanged(double)', self.foregroundWeightSpinBox.setValue)
    foregroundWeightSpinBox.connect('valueChanged(double)', self.foregroundPopupSlider.setValue)
    self.connections.append( (self.foregroundWeightSpinBox, 'valueChanged(double)', self.updateMRMLFromGUI ) )
    foregroundLayout.addWidget( foregroundLabel )
    foregroundLayout.addWidget( foregroundWeightSpinBox )
    foregroundPopupLayout.addWidget( foregroundPopupSlider )
    objectFormLayout.addRow("Foreground Label:", foregroundLayout )
    self.objectLabel = foregroundLabel
    # http://qt-project.org/doc/qt-4.7/qt.html
    foregroundPopup.alignment = 0x0082 # Qt::AlignVCenter | Qt::AlignRight
    foregroundPopup.horizontalDirection = 0 # Qt::LeftToRight
    foregroundPopup.verticalDirection = 0 # Qt::TopToBottom
    foregroundPopup.animationEffect = 1 # Qt::ScrollEffect

    backgroundLayout = qt.QHBoxLayout()
    backgroundLabel = slicer.qMRMLLabelComboBox()
    backgroundLabel.objectName = 'Background label'
    backgroundLabel.setMRMLScene(slicer.mrmlScene)
    backgroundLabel.setMRMLColorNode(self.editUtil.getColorNode())
    backgroundLabel.labelValueVisible = True
    backgroundLabel.currentColor = 2
    self.backgroundLabel = backgroundLabel
    self.connections.append( (self.backgroundLabel, 'currentColorChanged(int)', self.updateMRMLFromGUI ) )
    backgroundWeightSpinBox = qt.QDoubleSpinBox(backgroundLabel)
    self.backgroundWeightSpinBox = backgroundWeightSpinBox
    backgroundWeightSpinBox.setRange(0.0, 1.0)
    backgroundWeightSpinBox.setSingleStep(0.1)
    backgroundWeightSpinBox.value = 1.0
    backgroundPopup = ctk.ctkPopupWidget( backgroundWeightSpinBox )
    backgroundPopupLayout = qt.QHBoxLayout( backgroundPopup )
    backgroundPopupSlider = ctk.ctkDoubleSlider( backgroundPopup )
    self.backgroundPopupSlider = backgroundPopupSlider
    backgroundPopupSlider.maximum = 1.0
    backgroundPopupSlider.minimum = 0.0
    backgroundPopupSlider.singleStep = 0.1
    backgroundPopupSlider.connect('valueChanged(double)', self.backgroundWeightSpinBox.setValue)
    backgroundWeightSpinBox.connect('valueChanged(double)', self.backgroundPopupSlider.setValue)
    self.connections.append( (self.backgroundWeightSpinBox, 'valueChanged(double)', self.updateMRMLFromGUI ) )
    backgroundLayout.addWidget( backgroundLabel )
    backgroundLayout.addWidget( backgroundWeightSpinBox )
    backgroundPopupLayout.addWidget( backgroundPopupSlider )
    objectFormLayout.addRow("Background Label:", backgroundLayout)
    self.backgroundLabel = backgroundLabel
    # http://qt-project.org/doc/qt-4.7/qt.html
    backgroundPopup.alignment = 0x0082 # Qt::AlignVCenter | Qt::AlignRight
    backgroundPopup.horizontalDirection = 0 # Qt::LeftToRight
    backgroundPopup.verticalDirection = 0 # Qt::TopToBottom
    backgroundPopup.animationEffect = 1 # Qt::ScrollEffect

    # Presets
    # Placeholder
    presetsCollapsibleGroupBox = ctk.ctkCollapsibleGroupBox()
    presetsCollapsibleGroupBox.title = "Preset"
    self.frame.layout().addWidget(presetsCollapsibleGroupBox)
    presetComboBox = slicer.qSlicerPresetComboBox()

    # Advanced Parameters
    paramsCollapsibleGroupBox = ctk.ctkCollapsibleGroupBox()
    paramsCollapsibleGroupBox.title = "Advanced Parameters"
    paramsCollapsibleGroupBox.collapsed = 1
    self.frame.layout().addWidget(paramsCollapsibleGroupBox)

    # Layout within the io collapsible button
    paramsFormLayout = qt.QFormLayout(paramsCollapsibleGroupBox)

    erosionSpinBox = qt.QSpinBox()
    erosionSpinBox.objectName = 'erosionSpinBox'
    erosionSpinBox.toolTip = "Set the erosion radius."
    erosionSpinBox.setMinimum(0)
    erosionSpinBox.setValue(5) # Default
    paramsFormLayout.addRow("Erosion Radius:", erosionSpinBox)
    self.erosionSpinBox = erosionSpinBox
    self.connections.append( (self.erosionSpinBox, "valueChanged(int)", self.updateMRMLFromGUI ) )

    holeFillSpinBox = qt.QSpinBox()
    holeFillSpinBox.objectName = 'holeFillSpinBox'
    holeFillSpinBox.toolTip = "Set the hole fill iterations."
    holeFillSpinBox.setMinimum(0)
    holeFillSpinBox.setValue(5) #Default
    paramsFormLayout.addRow("Hole Fill Iterations:", holeFillSpinBox)
    self.holeFillSpinBox = holeFillSpinBox
    self.connections.append( (self.holeFillSpinBox, "valueChanged(int)", self.updateMRMLFromGUI ) )

    # probabilitySmoothingStandardDeviation spin box
    probabilitySmoothingStdDevSpinBox = qt.QDoubleSpinBox()
    probabilitySmoothingStdDevSpinBox.objectName = 'probabilitySmoothingStdDevSpinBox'
    probabilitySmoothingStdDevSpinBox.toolTip = "Standard deviation of blur applied to probability images prior to computing maximum likelihood of each class at each pixel."
    probabilitySmoothingStdDevSpinBox.setMinimum(0.0)
    probabilitySmoothingStdDevSpinBox.setValue(1.0) # Default
    probabilitySmoothingStdDevSpinBox.setSingleStep(0.5)
    paramsFormLayout.addRow("Probability Smoothing Standard Deviation:", probabilitySmoothingStdDevSpinBox)
    self.probabilitySmoothingStdDevSpinBox = probabilitySmoothingStdDevSpinBox
    self.connections.append( (self.probabilitySmoothingStdDevSpinBox, "valueChanged(double)", self.updateMRMLFromGUI ) )

    # histogramSmoothingStandardDeviation spin box
    histogramSmoothingStdDevSpinBox = qt.QDoubleSpinBox()
    histogramSmoothingStdDevSpinBox.objectName = 'histogramSmoothingStdDevSpinBox'
    histogramSmoothingStdDevSpinBox.toolTip = "Standard deviation of blur applied to histograms to convert them to probability density function estimates."
    histogramSmoothingStdDevSpinBox.setMinimum(0.0)
    histogramSmoothingStdDevSpinBox.setValue(5.0) # Default
    histogramSmoothingStdDevSpinBox.setSingleStep(0.5)
    paramsFormLayout.addRow("Probability Smoothing Standard Deviation:", histogramSmoothingStdDevSpinBox)
    self.histogramSmoothingStdDevSpinBox = histogramSmoothingStdDevSpinBox
    self.connections.append( (self.histogramSmoothingStdDevSpinBox, "valueChanged(double)", self.updateMRMLFromGUI ) )

    # draft check box
    draftCheckBox = qt.QCheckBox()
    draftCheckBox.objectName = 'draftCheckBox'
    draftCheckBox.toolTip = "Downsamples results by a factor of 4."
    paramsFormLayout.addRow("Draft Mode:", draftCheckBox)
    self.draftCheckBox = draftCheckBox
    self.connections.append( (self.draftCheckBox, "stateChanged(int)", self.updateMRMLFromGUI ) )

    # force classification check box
    forceClassificationCheckBox = qt.QCheckBox()
    forceClassificationCheckBox.objectName = 'forceClassificationCheckBox'
    forceClassificationCheckBox.toolTip = "Perform the classification of all voxels?"
    forceClassificationCheckBox.setChecked(False)
    paramsFormLayout.addRow("Classify all voxels: ", forceClassificationCheckBox)
    self.forceClassificationCheckBox = forceClassificationCheckBox
    self.connections.append( (self.forceClassificationCheckBox, "stateChanged(int)", self.updateMRMLFromGUI ) )

    # dilate first check box
    dilateFirstCheckBox = qt.QCheckBox()
    dilateFirstCheckBox.objectName = 'dilateFirstCheckBox'
    dilateFirstCheckBox.toolTip = "Dilate and then erode so as to fill-in holes?"
    dilateFirstCheckBox.setChecked(False)
    paramsFormLayout.addRow("Dilate First: ", dilateFirstCheckBox)
    self.dilateFirstCheckBox = dilateFirstCheckBox
    self.connections.append( (self.dilateFirstCheckBox, "stateChanged(int)", self.updateMRMLFromGUI ) )

    self.helpLabel = qt.QLabel("Run the PDF Segmentation on the current label map.", self.frame)
    self.frame.layout().addWidget(self.helpLabel)

    self.apply = qt.QPushButton("Apply", self.frame)
    self.apply.setToolTip("Apply to run segmentation.\nCreates a new label volume using the current volume as input")
    self.frame.layout().addWidget(self.apply)
    self.widgets.append(self.apply)

    EditorLib.HelpButton(self.frame, "Use this tool to apply segmentation using Parzen windowed PDFs.\n\n Select different label colors and paint on the foreground or background voxels using the paint effect.\nTo run segmentation correctly, you need to supply a minimum or two class labels.")

    self.connections.append( (self.apply, 'clicked()', self.onApply) )

  def destroy(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).destroy()

  # note: this method needs to be implemented exactly as-is
  # in each leaf subclass so that "self" in the observer
  # is of the correct type
  def updateParameterNode(self, caller, event):
    node = EditUtil.EditUtil().getParameterNode()
    if node != self.parameterNode:
      if self.parameterNode:
        node.RemoveObserver(self.parameterNodeTag)
      self.parameterNode = node
      self.parameterNodeTag = node.AddObserver(vtk.vtkCommand.ModifiedEvent, self.updateGUIFromMRML)

  def setMRMLDefaults(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).setMRMLDefaults()
    disableState = self.parameterNode.GetDisableModifiedEvent()
    self.parameterNode.SetDisableModifiedEvent(1)

    defaults = [
      ("outputVolume", "0"),
      ("labelmap", "0"),
      ("objectId", "1,2"),
      ("erodeRadius", "5"),
      ("holeFillIterations", "5"),
      ("objectPDFWeight", "1.0,1.0"),
      ("probImageSmoothingStdDev", "1.0"),
      ("histogramSmoothingStdDev", "5.0"),
      ("draft", "0"),
      ("forceClassification", "0"),
      ("dilateFirst", "0"),
      ]
    for i in range(0, 2):
      defaults.append(("additionalInputVolumeID" + str(i), "0"))

    # Set logic here because this function is called before the end
    # of the superclass constructor
    self.logic = InteractiveConnectedComponentsUsingParzenPDFsLogic(None)

    for default in defaults:
      pvalue = self.getParameter(default[0])
      if pvalue == "":
        self.setParameter(default[0], default[1])
    self.parameterNode.SetDisableModifiedEvent(disableState)

  def updateGUIFromMRML(self,caller,event):
    parameters = ["objectId",
                  "erodeRadius",
                  "holeFillIterations",
                  "objectPDFWeight",
                  "probImageSmoothingStdDev",
                  "histogramSmoothingStdDev",
                  "draft",
                  "forceClassification",
                  "dilateFirst",
                  ]
    for i in range(0, 2):
      parameters.append("additionalInputVolumeID" + str(i))

    for parameter in parameters:
      if self.getParameter(parameter) == "":
        # don't update if the parameter node has not got all values yet
        return
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).updateGUIFromMRML(caller,event)
    self.disconnectWidgets()

    # Additional inputs
    for i in range(0, 2):
      self.additionalInputNodeSelectors[i].currentNodeID = self.getParameter("additionalInputVolumeID" + str(i))

    # labels
    objectIds = self.logic.listFromStringList(self.getParameter("objectId"))
    self.foregroundLabel.currentColor = objectIds[0]
    self.backgroundLabel.currentColor = objectIds[1]

    # Parameters
    self.erosionSpinBox.value = int(self.getParameter("erodeRadius"))
    self.holeFillSpinBox.value = int(self.getParameter("holeFillIterations"))
    self.probabilitySmoothingStdDevSpinBox.value = float(self.getParameter("probImageSmoothingStdDev"))
    self.histogramSmoothingStdDevSpinBox.value = float(self.getParameter("histogramSmoothingStdDev"))
    self.draftCheckBox.setChecked(int(self.getParameter("draft")))
    self.forceClassificationCheckBox.setChecked(int(self.getParameter("forceClassification")))
    self.dilateFirstCheckBox.setChecked(int(self.getParameter("dilateFirst")))

    self.connectWidgets()

  def onApply(self):
    self.undoRedo.saveState()
    self.logic.applyPDFSegmenter()

  def updateMRMLFromGUI(self):
    if self.updatingGUI:
      return

    disableState = self.parameterNode.GetDisableModifiedEvent()
    self.parameterNode.SetDisableModifiedEvent(1)
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).updateMRMLFromGUI()

    # Input
    for i in range(0, 2):
      self.setParameter("additionalInputVolumeID" + str(i), self.additionalInputNodeSelectors[i].currentNodeID)

    # Labels
    objectIds = (str(self.foregroundLabel.currentColor) + ","
                 + str(self.backgroundLabel.currentColor)
                )
    self.setParameter("objectId", objectIds)

    # Parameters
    self.setParameter("erodeRadius", self.erosionSpinBox.text)
    self.setParameter("holeFillIterations", self.holeFillSpinBox.text)
    self.setParameter("probImageSmoothingStdDev", self.probabilitySmoothingStdDevSpinBox.text)
    self.setParameter("histogramSmoothingStdDev", self.histogramSmoothingStdDevSpinBox.text)
    self.setParameter("draft", str(int(self.draftCheckBox.isChecked())))
    self.setParameter("forceClassification", str(int(self.forceClassificationCheckBox.isChecked())))
    self.setParameter("dilateFirst", str(int(self.dilateFirstCheckBox.isChecked())))

    self.parameterNode.SetDisableModifiedEvent(disableState)
    if not disableState:
      self.parameterNode.InvokePendingModifiedEvent()

  def addInputNodeSelector(self, index, layout):
    inputNodeSelector =  slicer.qMRMLNodeComboBox()
    inputNodeSelector.objectName = 'additionalInputNodeSelector'+str(index+1)
    inputNodeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    inputNodeSelector.noneEnabled = True
    inputNodeSelector.addEnabled = False
    inputNodeSelector.removeEnabled = False
    inputNodeSelector.editEnabled = True
    inputNodeSelector.enabled = 1
    inputNodeSelector.setMRMLScene(slicer.mrmlScene)
    layout.addRow("Additional Input Volume "+str(index+1)+":", inputNodeSelector)
    self.connections.append( (inputNodeSelector, "currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI ) )
    return inputNodeSelector

  def setParameter(self, parameterName, value):
    self.logic.setParameter(parameterName, value)

  def getParameter(self, parameterName):
    return self.logic.getParameter(parameterName)

#
# EditorEffectTemplateTool
#

class InteractiveConnectedComponentsUsingParzenPDFsTool(LabelEffect.LabelEffectTool):
  """
  One instance of this will be created per-view when the effect
  is selected.  It is responsible for implementing feedback and
  label map changes in response to user input.
  This class observes the editor parameter node to configure itself
  and queries the current view for background and label volume
  nodes to operate on.
  """

  def __init__(self, sliceWidget):
    super(InteractiveConnectedComponentsUsingParzenPDFsTool,self).__init__(sliceWidget)
    # create a logic instance to do the non-gui work
    self.logic = InteractiveConnectedComponentsUsingParzenPDFsLogic(self.sliceWidget.sliceLogic())

  def cleanup(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsTool,self).cleanup()

  def processEvent(self, caller=None, event=None):
    """
    handle events from the render window interactor
    """

    # let the superclass deal with the event if it wants to
    if super(InteractiveConnectedComponentsUsingParzenPDFsTool,self).processEvent(caller,event):
      return

    if event == "LeftButtonPressEvent":
      xy = self.interactor.GetEventPosition()
      sliceLogic = self.sliceWidget.sliceLogic()
      logic = InteractiveConnectedComponentsUsingParzenPDFsLogic(sliceLogic)
      logic.apply(xy)
      print("Got a %s at %s in %s", (event,str(xy),self.sliceWidget.sliceLogic().GetSliceNode().GetName()))
      self.abortEvent(event)
    else:
      pass
    #if event == "LeftButtonPressEvent"
    #  self.actionState = "painting"
    #  if not self.pixelMode:
    #    self.cursorOff()
    #  xy = self.interactor.GetEventPosition()
    #elif event == "LeftButtonReleaseEvent":
    #  self.paintApply

    # events from the slice node
    #if caller and caller.IsA('vtkMRMLSliceNode'):
      # here you can respond to pan/zoom or other changes
      # to the view
    #  pass


#
# EditorEffectTemplateLogic
#

class InteractiveConnectedComponentsUsingParzenPDFsLogic(LabelEffect.LabelEffectLogic):
  """
  This class contains helper methods for a given effect
  type.  It can be instanced as needed by an EditorEffectTemplateTool
  or EditorEffectTemplateOptions instance in order to compute intermediate
  results (say, for user feedback) or to implement the final
  segmentation editing operation.  This class is split
  from the EditorEffectTemplateTool so that the operations can be used
  by other code without the need for a view context.
  """

  def __init__(self,sliceLogic):
    super(InteractiveConnectedComponentsUsingParzenPDFsLogic,self).__init__(sliceLogic)
    self.effectName = 'InteractiveConnectedComponentsUsingParzenPDFsOptions'
    self.parameterNode = self.editUtil.getParameterNode()

  def getCLINode(self, module, nodeName):
    cliNode = slicer.mrmlScene.GetFirstNodeByName(nodeName)
    # Also check path to make sure the CLI isn't a scripted module
    if (cliNode == None) and ("qt-scripted-modules" not in module.path):
      cliNode = slicer.cli.createNode(module)
      cliNode.SetName(nodeName)
    return cliNode

  def setParameter(self, parameterName, value):
    self.parameterNode.SetParameter(self.getFullParameterName(parameterName), value)

  def getParameter(self, parameterName):
    return self.parameterNode.GetParameter(self.getFullParameterName(parameterName))

  def getFullParameterName(self, parameterName):
    return self.effectName + ',' + parameterName

  def listFromStringList(self, stringlist):
    '''Convert a stringlist of the format '1.0, 2.0, 3.0' to a list
       of the format [1.0, 2.0, 3.0].'''
    list = []
    for string in stringlist.split(","):
      try:
        list.append(int(string))
      except ValueError:
        list.append(float(string))
    return list

  def apply(self,xy):
    pass

  def applyPDFSegmenter(self):
    #
    # Apply PDF segmenter based on the parameter node
    #
    if not self.sliceLogic:
      self.sliceLogic = self.editUtil.getSliceLogic()
    cliParameters = {}

    # IO
    cliParameters["inputVolume1"] = self.editUtil.getBackgroundVolume()
    for i in range(0,2):
      # Get input nodes by their IDs
      nodeID = self.getParameter("additionalInputVolumeID" + str(i))
      cliParameters["inputVolume"+str(i+2)] = slicer.mrmlScene.GetNodeByID(nodeID)

    cliParameters["labelmap"] = self.editUtil.getLabelVolume()
    cliParameters["outputVolume"] = self.editUtil.getLabelVolume()

    # Labels
    cliParameters["objectId"] = self.listFromStringList(self.getParameter("objectId"))

    # Parameters
    cliParameters["erodeRadius"] = int(self.getParameter( "erodeRadius"))
    cliParameters["holeFillIterations"] = int(self.getParameter("holeFillIterations"))
    cliParameters["objectPDFWeight"] = self.listFromStringList(self.getParameter("objectPDFWeight"))
    cliParameters["probImageSmoothingStdDev"] = float(self.getParameter("probImageSmoothingStdDev"))
    cliParameters["histogramSmoothingStdDev"] = float(self.getParameter("histogramSmoothingStdDev"))
    cliParameters["draft"] = int(self.getParameter("draft"))
    cliParameters["forceClassification"] = int(self.getParameter("forceClassification"))
    cliParameters["dilateFirst"] = int(self.getParameter("dilateFirst"))

    module = slicer.modules.segmentconnectedcomponentsusingparzenpdfs
    cliNode = self.getCLINode(module, "PDFSegmenterEditorEffect")
    slicer.cli.run(module, cliNode, cliParameters)

#
# The InteractiveConnectedComponentsUsingParzenPDFs Template class definition
#

class InteractiveConnectedComponentsUsingParzenPDFsExtension(LabelEffect.LabelEffect):
  """Organizes the Options, Tool, and Logic classes into a single instance
  that can be managed by the EditBox
  """

  def __init__(self):
    # name is used to define the name of the icon image resource (e.g. EditorEffectTemplate.png)
    self.name = "InteractiveConnectedComponentsUsingParzenPDFs"
    self.title = "InteractiveConnectedComponentsUsingParzenPDFs"
    # tool tip is displayed on mouse hover
    self.toolTip = "Perform PDF Segmentation"

    self.options = InteractiveConnectedComponentsUsingParzenPDFsOptions
    self.tool = InteractiveConnectedComponentsUsingParzenPDFsTool
    self.logic = InteractiveConnectedComponentsUsingParzenPDFsLogic


#
# EditorEffectTemplate
#

class InteractiveConnectedComponentsUsingParzenPDFs:
  """
  This class is the 'hook' for slicer to detect and recognize the extension
  as a loadable scripted module
  """
  def __init__(self, parent):
    parent.title = "Editor InteractiveConnectedComponentsUsingParzenPDFs Effect"
    parent.categories = ["Developer Tools.Editor Extensions"]
    parent.contributors = ["Danielle Pace (Kitware)",
                           "Christopher Mullins (Kitware)",
                           "Stephen Aylward (Kitware)",
                           "Johan Andruejol (Kitware)",]
    parent.helpText = """
    The PDF Segmenter is a framework for using connected components in
    conjunction with intensity histograms for classifying images in pixel space.

    This module is available as an editor tool via the editor module in Slicer.
    This module cannot be run as a standard module in Slicer.
    """
    parent.acknowledgementText = """
    This work is part of the TubeTK project at Kitware.
    Module implemented by Danielle Pace.  PDF Segmenter implemented by
    Stephen Aylward.
    """

    # TODO:
    # don't show this module - it only appears in the Editor module
    #parent.hidden = True

    # Add this extension to the editor's list for discovery when the module
    # is created.  Since this module may be discovered before the Editor itself,
    # create the list if it doesn't already exist.
    try:
      slicer.modules.editorExtensions
    except AttributeError:
      slicer.modules.editorExtensions = {}
    slicer.modules.editorExtensions['InteractiveConnectedComponentsUsingParzenPDFs'] = InteractiveConnectedComponentsUsingParzenPDFsExtension

#
# EditorEffectTemplateWidget
#

class InteractiveConnectedComponentsUsingParzenPDFsWidget:
  def __init__(self, parent = None):
    self.parent = parent

  def setup(self):
    # don't display anything for this widget - it will be hidden anyway
    pass

  def enter(self):
    pass

  def exit(self):
    pass
