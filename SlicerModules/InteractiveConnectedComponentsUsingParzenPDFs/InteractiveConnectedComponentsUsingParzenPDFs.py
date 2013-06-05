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

  def __del__(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).__del__()

  def create(self):
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).create()

    ioCollapsibleButton = ctk.ctkCollapsibleButton()
    ioCollapsibleButton.text = "IO"
    self.frame.layout().addWidget(ioCollapsibleButton)

    # Layout within the io collapsible button
    ioFormLayout = qt.QFormLayout(ioCollapsibleButton)

    inputNodeSelector1 =  slicer.qMRMLNodeComboBox()
    inputNodeSelector1.objectName = 'inputNodeSelector1'
    inputNodeSelector1.toolTip = "Select the 1st input volume to be segmented."
    inputNodeSelector1.nodeTypes = ['vtkMRMLScalarVolumeNode']
    inputNodeSelector1.noneEnabled = False
    inputNodeSelector1.addEnabled = False
    inputNodeSelector1.removeEnabled = False
    inputNodeSelector1.editEnabled = True
    inputNodeSelector1.enabled = 1
    inputNodeSelector1.setMRMLScene(slicer.mrmlScene)
    ioFormLayout.addRow("Input Volume 1:", inputNodeSelector1)
    self.inputNodeSelector1 = inputNodeSelector1

    inputNodeSelector2 = slicer.qMRMLNodeComboBox()
    inputNodeSelector2.objectName = 'inputNodeSelector2'
    inputNodeSelector2.toolTip = "Select the 2nd coregistered input volume to be segmented."
    inputNodeSelector2.nodeTypes = ['vtkMRMLScalarVolumeNode']
    inputNodeSelector2.noneEnabled = True
    inputNodeSelector2.addEnabled = False
    inputNodeSelector2.removeEnabled = False
    inputNodeSelector2.editEnabled = True
    inputNodeSelector2.enabled = 1
    inputNodeSelector2.setMRMLScene(slicer.mrmlScene)
    ioFormLayout.addRow("Input Volume 2:", inputNodeSelector2)
    self.inputNodeSelector2 = inputNodeSelector2

    inputNodeSelector3 = slicer.qMRMLNodeComboBox()
    inputNodeSelector3.objectName = 'inputNodeSelector3'
    inputNodeSelector3.toolTip = "Select the 3rd coregistered input volume to be segmented."
    inputNodeSelector3.nodeTypes = ['vtkMRMLScalarVolumeNode']
    inputNodeSelector3.noneEnabled = True
    inputNodeSelector3.addEnabled = False
    inputNodeSelector3.removeEnabled = False
    inputNodeSelector3.editEnabled = True
    inputNodeSelector3.enabled = 1
    inputNodeSelector3.setMRMLScene(slicer.mrmlScene)
    ioFormLayout.addRow("Input Volume 3:", inputNodeSelector3)
    self.inputNodeSelector3 = inputNodeSelector3

    # Objects
    objectCollapsibleButton = ctk.ctkCollapsibleButton()
    objectCollapsibleButton.text = "Objects"
    self.frame.layout().addWidget(objectCollapsibleButton)

    # Layout within the io collapsible button
    objectFormLayout = qt.QFormLayout(objectCollapsibleButton)
    objectLabel = slicer.qMRMLLabelComboBox()
    objectLabel.objectName = 'Foreground label selector'
    objectLabel.setMRMLScene(slicer.mrmlScene)
    objectLabel.setMRMLColorNode(slicer.modules.EditorWidget.editUtil.getColorNode())
    objectLabel.labelValueVisible = True
    objectLabel.currentColor = 1
    objectFormLayout.addRow("Foreground Label:", objectLabel)
    self.objectLabel = objectLabel

    backgroundLabel = slicer.qMRMLLabelComboBox()
    backgroundLabel.objectName = 'Background label selector'
    backgroundLabel.setMRMLScene(slicer.mrmlScene)
    backgroundLabel.setMRMLColorNode(slicer.modules.EditorWidget.editUtil.getColorNode())
    backgroundLabel.labelValueVisible = True
    backgroundLabel.currentColor = 2
    objectFormLayout.addRow("Background Label:", backgroundLabel)
    self.backgroundLabel = backgroundLabel

    barrierLabel = slicer.qMRMLLabelComboBox()
    barrierLabel.objectName = 'Barrier label selector'
    barrierLabel.setMRMLScene(slicer.mrmlScene)
    barrierLabel.setMRMLColorNode(slicer.modules.EditorWidget.editUtil.getColorNode())
    barrierLabel.labelValueVisible = True
    barrierLabel.currentColor = 3
    objectFormLayout.addRow("Barrier Label:", barrierLabel)
    self.barrierLabel = barrierLabel

    voidLabel = slicer.qMRMLLabelComboBox()
    voidLabel.objectName = 'Void label selector'
    voidLabel.setMRMLScene(slicer.mrmlScene)
    voidLabel.setMRMLColorNode(slicer.modules.EditorWidget.editUtil.getColorNode())
    voidLabel.labelValueVisible = True
    # HACK: Set this to 1 first to get around a bug in the labelComboBox MRML interaction
    #voidLabel.setCurrentColor(1)
    voidLabel.setCurrentColor(0)
    objectFormLayout.addRow("Void Label:", voidLabel)
    self.voidLabel = voidLabel
    self.voidLabel.currentColor = 0

    # Presets
    # Placeholder
    presetsCollapsibleButton = ctk.ctkCollapsibleButton()
    presetsCollapsibleButton.text = "Preset"
    self.frame.layout().addWidget(presetsCollapsibleButton)
    presetComboBox = slicer.qSlicerPresetComboBox()

    # Advanced Parameters
    paramsCollapsibleButton = ctk.ctkCollapsibleButton()
    paramsCollapsibleButton.text = "Advanced Parameters"
    self.frame.layout().addWidget(paramsCollapsibleButton)

    # Layout within the io collapsible button
    paramsFormLayout = qt.QFormLayout(paramsCollapsibleButton)

    erosionSpinBox = qt.QSpinBox()
    erosionSpinBox.objectName = 'erosionSpinBox'
    erosionSpinBox.toolTip = "Set the erosion radius."
    erosionSpinBox.setMinimum(0)
    paramsFormLayout.addRow("Erosion Radius:", erosionSpinBox)
    self.erosionSpinBox = erosionSpinBox

    holeFillSpinBox = qt.QSpinBox()
    holeFillSpinBox.objectName = 'holeFillSpinBox'
    holeFillSpinBox.toolTip = "Set the hole fill iterations."
    holeFillSpinBox.setMinimum(0)
    paramsFormLayout.addRow("Hole Fill Iterations:", holeFillSpinBox)
    self.holeFillSpinBox = holeFillSpinBox

    falsePositiveRatioSpinBox = qt.QDoubleSpinBox()
    falsePositiveRatioSpinBox.objectName = 'falsePositiveRatioSpinBox'
    falsePositiveRatioSpinBox.toolTip = "Relative Cost of False Positive vs. false negative."
    falsePositiveRatioSpinBox.setMinimum(0.0)
    falsePositiveRatioSpinBox.setValue(1.0) # Default
    falsePositiveRatioSpinBox.setSingleStep(0.1)
    paramsFormLayout.addRow("False Positive Ratio:", falsePositiveRatioSpinBox)
    self.falsePositiveRatioSpinBox = falsePositiveRatioSpinBox

    # probabilitySmoothingStandardDeviation spin box
    probabilitySmoothingStdDevSpinBox = qt.QDoubleSpinBox()
    probabilitySmoothingStdDevSpinBox.objectName = 'probabilitySmoothingStdDevSpinBox'
    probabilitySmoothingStdDevSpinBox.toolTip = "Standard deviation of blur applied to probability images prior to computing maximum likelihood of each class at each pixel."
    probabilitySmoothingStdDevSpinBox.setMinimum(0.0)
    probabilitySmoothingStdDevSpinBox.setValue(3.0) # Default
    probabilitySmoothingStdDevSpinBox.setSingleStep(0.1)
    paramsFormLayout.addRow("Probability Smoothing Standard Deviation:", probabilitySmoothingStdDevSpinBox)
    self.probabilitySmoothingStdDevSpinBox = probabilitySmoothingStdDevSpinBox

    # draft check box
    draftCheckBox = qt.QCheckBox()
    draftCheckBox.objectName = 'draftCheckBox'
    draftCheckBox.toolTip = "Downsamples results by a factor of 4."
    paramsFormLayout.addRow("Draft Mode:", draftCheckBox)
    self.draftCheckBox = draftCheckBox

    # reclassifyObjectMask check box
    reclassifyObjectMaskCheckBox = qt.QCheckBox()
    reclassifyObjectMaskCheckBox.objectName = 'reclassifyObjectMaskCheckBox'
    reclassifyObjectMaskCheckBox.toolTip = "Perform classification on voxels within the foreground mask?"
    reclassifyObjectMaskCheckBox.setChecked(True)
    paramsFormLayout.addRow("Reclassify Foreground Mask:", reclassifyObjectMaskCheckBox)
    self.reclassifyObjectMaskCheckBox = reclassifyObjectMaskCheckBox

    # reclassifyNotObjectMask check box
    reclassifyNotObjectMaskCheckBox = qt.QCheckBox()
    reclassifyNotObjectMaskCheckBox.objectName = 'reclassifyNotObjectMaskCheckBox'
    reclassifyNotObjectMaskCheckBox.toolTip = "Perform classification on voxels within the barrier mask?"
    reclassifyNotObjectMaskCheckBox.setChecked(True)
    paramsFormLayout.addRow("Reclassify Barrier Mask:", reclassifyNotObjectMaskCheckBox)
    self.reclassifyNotObjectMaskCheckBox = reclassifyNotObjectMaskCheckBox

    # force classification check box
    forceClassificationCheckBox = qt.QCheckBox()
    forceClassificationCheckBox.objectName = 'forceClassificationCheckBox'
    forceClassificationCheckBox.toolTip = "Perform the classification of all voxels?"
    forceClassificationCheckBox.setChecked(False)
    paramsFormLayout.addRow("Force Classification: ", forceClassificationCheckBox)
    self.forceClassificationCheckBox = forceClassificationCheckBox

    self.helpLabel = qt.QLabel("Run the PDF Segmentation on the current label map.", self.frame)
    self.frame.layout().addWidget(self.helpLabel)

    self.apply = qt.QPushButton("Apply", self.frame)
    self.apply.setToolTip("Apply to run segmentation.\nCreates a new label volume using the current volume as input")
    self.frame.layout().addWidget(self.apply)
    self.widgets.append(self.apply)

    EditorLib.HelpButton(self.frame, "Use this tool to apply segmentation using Parzen windowed PDFs.\n\n Select different label colors and paint on the foreground, background, or barrier voxels using the paint effect.\nTo run segmentation correctly, you need to supply a minimum or two class labels.")

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

  def updateGUIFromMRML(self,caller,event):
    self.disconnectWidgets()
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).updateGUIFromMRML(caller,event)
    self.connectWidgets()

  def onApply(self):

    self.undoRedo.saveState()
    objectIds = []
    objectIds.append(self.objectLabel.currentColor)
    objectIds.append(self.backgroundLabel.currentColor)
    parameters = {}
    parameters['inputVolume1'] = self.inputNodeSelector1.currentNode()
    parameters['inputVolume2'] = self.inputNodeSelector2.currentNode()
    parameters['inputVolume3'] = self.inputNodeSelector3.currentNode()
    parameters['voidId'] = self.voidLabel.currentColor
    parameters['objectId'] = objectIds
    parameters['labelmap'] = slicer.modules.EditorWidget.editUtil.getLabelVolume()
    parameters['outputVolume'] = slicer.modules.EditorWidget.editUtil.getLabelVolume()
    parameters['erodeRadius'] = self.erosionSpinBox.value
    parameters['holeFillIterations'] = self.holeFillSpinBox.value
    parameters['fprWeight'] = self.falsePositiveRatioSpinBox.value
    parameters['probSmoothingStdDev'] = self.probabilitySmoothingStdDevSpinBox.value
    parameters['draft'] = self.draftCheckBox.checked
    parameters['reclassifyObjectMask'] = self.reclassifyObjectMaskCheckBox.checked
    parameters['reclassifyNotObjectMask'] = self.reclassifyNotObjectMaskCheckBox.checked
    parameters['forceClassification'] = self.forceClassificationCheckBox.checked

    slicer.cli.run(slicer.modules.segmentconnectedcomponentsusingparzenpdfs, None, parameters)

  def updateMRMLFromGUI(self):
    if self.updatingGUI:
      return
    disableState = self.parameterNode.GetDisableModifiedEvent()
    self.parameterNode.SetDisableModifiedEvent(1)
    super(InteractiveConnectedComponentsUsingParzenPDFsOptions,self).updateMRMLFromGUI()
    self.parameterNode.SetDisableModifiedEvent(disableState)
    if not disableState:
      self.parameterNode.InvokePendingModifiedEvent()

  def setInputNode1(self):
    backgroundLogic = self.sliceLogic.GetBackgroundLayer()
    backgroundNode = backgroundLogic.GetVolumeNode()


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
    self.sliceLogic = sliceLogic

  def apply(self,xy):
    pass


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
                           "Stephen Aylward (Kitware)"]
    parent.helpText = """
    The PDF Segmenter is a framework for using connected components alonside intensity
    histograms for classifying images in pixel space.
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
