// QT includes
#include <QDebug>
#include <QColor>

// qMRML includes
#include "qSlicerSpatialObjectsBasicWidget.h"
#include "ui_qSlicerSpatialObjectsBasicWidget.h"

// MRML includes
#include <vtkMRMLStorageNode.h>
#include <vtkMRMLSpatialObjectsNode.h>
#include <vtkMRMLSpatialObjectsDisplayNode.h>

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>

//------------------------------------------------------------------------------
class qSlicerSpatialObjectsBasicWidgetPrivate :
public Ui_qSlicerSpatialObjectsBasicWidget
{
Q_DECLARE_PUBLIC(qSlicerSpatialObjectsBasicWidget);

protected:
  qSlicerSpatialObjectsBasicWidget* const q_ptr;

public:
  qSlicerSpatialObjectsBasicWidgetPrivate(
    qSlicerSpatialObjectsBasicWidget& object);
  void init();

  vtkMRMLSpatialObjectsNode* SpatialObjectsNode;
  vtkMRMLSpatialObjectsDisplayNode *LineDN, *TubeDN, *GlyphDN;
};

//------------------------------------------------------------------------------
qSlicerSpatialObjectsBasicWidgetPrivate::
qSlicerSpatialObjectsBasicWidgetPrivate
  (qSlicerSpatialObjectsBasicWidget& object)
  : q_ptr(&object)
{
  this->SpatialObjectsNode = 0;
  this->LineDN = 0;
  this->TubeDN = 0;
  this->GlyphDN = 0;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidgetPrivate::init()
{
  Q_Q(qSlicerSpatialObjectsBasicWidget);
  this->setupUi(q);

  QObject::connect(this->LineVisibility, SIGNAL(stateChanged(int)), q,
                   SLOT(setLineVisibility(int)) );
  QObject::connect(this->TubeVisibility, SIGNAL(stateChanged(int)), q,
                   SLOT(setTubeVisibility(int)) );
  QObject::connect(this->GlyphVisibility, SIGNAL(stateChanged(int)), q,
                   SLOT(setGlyphVisibility(int)) );
  }

//------------------------------------------------------------------------------
qSlicerSpatialObjectsBasicWidget::
qSlicerSpatialObjectsBasicWidget(QWidget *_parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSpatialObjectsBasicWidgetPrivate(*this))
{
  Q_D(qSlicerSpatialObjectsBasicWidget);
  d->init();
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsBasicWidget::~qSlicerSpatialObjectsBasicWidget()
{}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode* qSlicerSpatialObjectsBasicWidget::
spatialObjectsNode() const
{
  Q_D(const qSlicerSpatialObjectsBasicWidget);
  return d->SpatialObjectsNode;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidget::setSpatialObjectsNode(vtkMRMLNode* node)
{
  this->setSpatialObjectsNode(vtkMRMLSpatialObjectsNode::SafeDownCast(node));
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidget::
setSpatialObjectsNode(vtkMRMLSpatialObjectsNode* SpatialObjectsNode)
{
  Q_D(qSlicerSpatialObjectsBasicWidget);

  vtkMRMLSpatialObjectsNode *oldNode = this->spatialObjectsNode();
  vtkMRMLSpatialObjectsDisplayNode *oldLDN = NULL;
  vtkMRMLSpatialObjectsDisplayNode *oldTDN = NULL;
  vtkMRMLSpatialObjectsDisplayNode *oldGDN = NULL;

  if (d->SpatialObjectsNode)
    {
    oldLDN = d->SpatialObjectsNode->GetLineDisplayNode();
    oldTDN = d->SpatialObjectsNode->GetTubeDisplayNode();
    oldGDN = d->SpatialObjectsNode->GetGlyphDisplayNode();
    }

  if (SpatialObjectsNode)
    {
    d->LineDN = SpatialObjectsNode->GetLineDisplayNode();
    d->TubeDN = SpatialObjectsNode->GetTubeDisplayNode();
    d->GlyphDN = SpatialObjectsNode->GetGlyphDisplayNode();
    }

  d->SpatialObjectsNode = SpatialObjectsNode;

  qvtkReconnect(oldNode, d->SpatialObjectsNode,
                vtkCommand::ModifiedEvent, this,
                SLOT(updateWidgetFromMRML()));

  qvtkReconnect(oldLDN, d->LineDN,
                vtkCommand::ModifiedEvent, this,
                SLOT(updateWidgetFromMRML()) );
  qvtkReconnect(oldTDN, d->TubeDN,
                vtkCommand::ModifiedEvent, this,
                SLOT(updateWidgetFromMRML()));
  qvtkReconnect(oldGDN, d->GlyphDN,
                vtkCommand::ModifiedEvent, this,
                SLOT(updateWidgetFromMRML()));

  this->updateWidgetFromMRML();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidget::updateWidgetFromMRML()
{
  Q_D(qSlicerSpatialObjectsBasicWidget);

  if (!d->SpatialObjectsNode)
    {
    return;
    }

  if (d->SpatialObjectsNode->GetLineDisplayNode())
    {
    d->LineVisibility->setChecked(
      d->SpatialObjectsNode->GetLineDisplayNode()->GetVisibility());
    }

  if (d->SpatialObjectsNode->GetTubeDisplayNode())
    {
    d->TubeVisibility->setChecked(
      d->SpatialObjectsNode->GetTubeDisplayNode()->GetVisibility());
    }

  if (d->SpatialObjectsNode->GetGlyphDisplayNode())
    {
    d->GlyphVisibility->setChecked( 
      d->SpatialObjectsNode->GetGlyphDisplayNode()->GetVisibility());
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidget::setLineVisibility(int state)
{
  Q_D(qSlicerSpatialObjectsBasicWidget);

  if (!d->SpatialObjectsNode)
    {
    return;
    }

  vtkMRMLSpatialObjectsDisplayNode* dNode =
    d->SpatialObjectsNode->GetLineDisplayNode();
  if (dNode)
    {
    dNode->SetVisibility((state > 0));
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidget::setTubeVisibility(int state)
{
  Q_D(qSlicerSpatialObjectsBasicWidget);

  if (!d->SpatialObjectsNode)
    {
    return;
    }

  vtkMRMLSpatialObjectsDisplayNode* dNode =
    d->SpatialObjectsNode->GetTubeDisplayNode();
  if (dNode)
    {
    dNode->SetVisibility((state > 0));
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsBasicWidget::setGlyphVisibility(int state)
{
  Q_D(qSlicerSpatialObjectsBasicWidget);

  if (!d->SpatialObjectsNode)
    {
    return;
    }

  vtkMRMLSpatialObjectsDisplayNode* dNode =
    d->SpatialObjectsNode->GetGlyphDisplayNode();
  if (dNode)
    {
    dNode->SetVisibility((state > 0));
    }
}
