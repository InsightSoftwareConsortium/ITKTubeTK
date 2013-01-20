/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Michael Jeulin-L, Kitware Inc.

==============================================================================*/

#include "vtkObjectFactory.h"

// VTK includes
#include "vtkAssignAttribute.h"
#include "vtkCallbackCommand.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// MRML includes
#include "vtkMRMLScene.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLSpatialObjectsLineDisplayNode.h"
#include "vtkMRMLSpatialObjectsDisplayPropertiesNode.h"

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSpatialObjectsLineDisplayNode);


//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsLineDisplayNode::vtkMRMLSpatialObjectsLineDisplayNode()
{
  this->ColorMode = vtkMRMLSpatialObjectsDisplayNode::colorModeSolid;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsLineDisplayNode::~vtkMRMLSpatialObjectsLineDisplayNode()
{
  this->RemoveObservers (vtkCommand::ModifiedEvent, this->MRMLCallbackCommand);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::ReadXMLAttributes(const char** atts)
{
  Superclass::ReadXMLAttributes(atts);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::Copy(vtkMRMLNode *anode)
{
  Superclass::Copy(anode);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::PrintSelf(ostream& os,
                                                     vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::UpdatePolyDataPipeline() 
{
  if (!this->GetInputPolyData() || !this->Visibility)
    {
    return;
    }

  this->Superclass::UpdatePolyDataPipeline();

  // Set display properties according to the
  // line display properties node
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    SpatialObjectsDisplayPropertiesNode =
      this->GetSpatialObjectsDisplayPropertiesNode();

  if (SpatialObjectsDisplayPropertiesNode != NULL)
    {
    const int colorMode = this->GetColorMode();
    if (colorMode ==
          vtkMRMLSpatialObjectsDisplayNode::colorModeSolid)
      {
      this->ScalarVisibilityOff();

      vtkMRMLNode* colorNode =
        this->GetScene()->GetNodeByID("vtkMRMLColorTableNodeFullRainbow");
      if (colorNode)
        {
        this->SetAndObserveColorNodeID(colorNode->GetID());
        }

      this->AutoScalarRangeOff();
      this->SetScalarRange(0, 255);
      }
    else if (colorMode ==
               vtkMRMLSpatialObjectsDisplayNode::colorModeScalarData)
      {
      this->ScalarVisibilityOn();
      this->AssignAttribute->Update();
      }
    }
  else
    {
    this->ScalarVisibilityOff();
    }
}
