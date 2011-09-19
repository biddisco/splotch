/*=========================================================================

 Program:   Visualization Toolkit
 Module:    pqSplotchRaytraceDisplayPanelDecorator.cxx

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/

// .NAME pqSplotchRaytraceDisplayPanelDecorator
// .SECTION Thanks
// <verbatim>
//
//  This file is part of the SplotchRaytraces plugin developed and contributed by
//
// </verbatim>

#include "pqSplotchRaytraceDisplayPanelDecorator.h"
#include "ui_pqSplotchRaytraceDisplayPanelDecorator.h"

// Server Manager Includes.
#include "vtkCommand.h"
#include "vtkDataObject.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkSmartPointer.h"
#include "vtkSMProperty.h"
#include "vtkSMPropertyHelper.h"
#include "vtkSMPVRepresentationProxy.h"
#include "vtkPVDataInformation.h"
#include "vtkSMStringVectorProperty.h"
#include "vtkSMEnumerationDomain.h"

// Qt Includes.
#include <QVBoxLayout>
#include <QComboBox>
#include <QGroupBox>

// ParaView Includes.
#include "pqDisplayProxyEditor.h"
#include "pqRepresentation.h"
#include "pqFieldSelectionAdaptor.h"
#include "pqPropertyLinks.h"
#include "pqSMAdaptor.h"
#include "pqPipelineRepresentation.h"
#include "pqVariableType.h"
#include "pqScalarsToColors.h"
#include "pqWidgetRangeDomain.h"
#include "pqFieldSelectionAdaptor.h"
#include "pqSignalAdaptors.h"

class pqSplotchRaytraceDisplayPanelDecorator::pqInternals: public Ui::pqSplotchRaytraceDisplayPanelDecorator
{
public:
  pqPropertyLinks Links;
  vtkSMProxy                * Representation;
  vtkSMPVRepresentationProxy* RepresentationProxy;
  vtkSmartPointer<vtkEventQtSlotConnect> VTKConnect;
  pqPipelineRepresentation* PipelineRepresentation;

  pqWidgetRangeDomain* MaxPixelSizeRangeDomain;
  pqWidgetRangeDomain* OpacityRangeDomain;
  pqWidgetRangeDomain* RadiusRangeDomain;
  QWidget* Frame;

  pqInternals(QWidget* parent)
  {
    this->VTKConnect     = vtkSmartPointer<vtkEventQtSlotConnect>::New();
    this->Representation = 0;
    this->Frame          = 0;
/*
    this->RepresentationProxy = 0;
    this->TransferFunctionDialog = new pqTransferFunctionDialog(parent);
    this->MaxPixelSizeRangeDomain = NULL;
    this->OpacityRangeDomain = NULL;
    this->RadiusRangeDomain = NULL;
*/
  }
};

//-----------------------------------------------------------------------------
pqSplotchRaytraceDisplayPanelDecorator::pqSplotchRaytraceDisplayPanelDecorator(
  pqDisplayPanel* _panel):Superclass(_panel)
{
  pqDisplayProxyEditor *panel = qobject_cast<pqDisplayProxyEditor*> (_panel);
  pqRepresentation     *repr = panel->getRepresentation();
  vtkSMProxy      *reprProxy = (repr) ? repr->getProxy() : NULL;

  this->Internals = NULL;
  vtkSMProperty* prop = reprProxy->GetProperty("IntensityScalars");
  if (prop)
    {
    this->Internals = new pqInternals(this);
    this->Internals->Representation = reprProxy;

    QWidget* wid = new QWidget(panel);
    this->Internals->Frame = wid;
    this->Internals->setupUi(wid);
    QVBoxLayout* l = qobject_cast<QVBoxLayout*>(panel->layout());
    l->addWidget(wid);
    //
    this->Internals->RepresentationProxy = vtkSMPVRepresentationProxy::SafeDownCast(reprProxy);

    pqPipelineRepresentation *pr = qobject_cast<pqPipelineRepresentation*>(panel->getRepresentation());

    //
    // Intensity
    //
    pqFieldSelectionAdaptor* adaptor= new pqFieldSelectionAdaptor(
      this->Internals->IntensityArray, prop);
    this->Internals->Links.addPropertyLink(
      adaptor, "attributeMode", SIGNAL(selectionChanged()),
      reprProxy, prop, 0);
    this->Internals->Links.addPropertyLink(
      adaptor, "scalar", SIGNAL(selectionChanged()),
      reprProxy, prop, 1);
    reprProxy->GetProperty("Input")->UpdateDependentDomains();
    prop->UpdateDependentDomains();

    //
    // Radius
    //
    prop = reprProxy->GetProperty("RadiusScalars");
    adaptor = new pqFieldSelectionAdaptor(
      this->Internals->RadiusArray, prop);
    this->Internals->Links.addPropertyLink(
      adaptor, "attributeMode", SIGNAL(selectionChanged()),
      reprProxy, prop, 0);
    this->Internals->Links.addPropertyLink(
      adaptor, "scalar", SIGNAL(selectionChanged()),
      reprProxy, prop, 1);
    prop->UpdateDependentDomains();

    //
    // Type
    //
    prop = reprProxy->GetProperty("TypeScalars");
    adaptor = new pqFieldSelectionAdaptor(
      this->Internals->TypeArray, prop);
    this->Internals->Links.addPropertyLink(
      adaptor, "attributeMode", SIGNAL(selectionChanged()),
      reprProxy, prop, 0);
    this->Internals->Links.addPropertyLink(
      adaptor, "scalar", SIGNAL(selectionChanged()),
      reprProxy, prop, 1);
    prop->UpdateDependentDomains();

    //
    // Active
    //
    prop = reprProxy->GetProperty("ActiveScalars");
    adaptor = new pqFieldSelectionAdaptor(
      this->Internals->ActiveArray, prop);
    this->Internals->Links.addPropertyLink(
      adaptor, "attributeMode", SIGNAL(selectionChanged()),
      reprProxy, prop, 0);
    this->Internals->Links.addPropertyLink(
      adaptor, "scalar", SIGNAL(selectionChanged()),
      reprProxy, prop, 1);
    prop->UpdateDependentDomains();

    //
    //
    //
    this->Internals->Links.addPropertyLink(
      this->Internals->brightness, "value", SIGNAL(valueChanged(int)),
      reprProxy, reprProxy->GetProperty("Brightness"));

  }

  this->setupGUIConnections();
}
//-----------------------------------------------------------------------------
pqSplotchRaytraceDisplayPanelDecorator::~pqSplotchRaytraceDisplayPanelDecorator()
{
  delete this->Internals;
  this->Internals = 0;
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::setRepresentation(
    pqPipelineRepresentation* repr)
{
//  this->Internals->IntensityArray->setRepresentation(repr);

//  QObject::connect(this->Internals->ScaleBy, SIGNAL(modified()), this,
//      SLOT(updateEnableState()), Qt::QueuedConnection);

}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::setupGUIConnections()
{
  this->Internals->VTKConnect->Connect(
      this->Internals->RepresentationProxy->GetProperty("Representation"),
      vtkCommand::ModifiedEvent, this, SLOT(representationTypeChanged()));

//  this->connect(this->Internals->IntensityArray,
//      SIGNAL(variableChanged(pqVariableType, const QString&)), this,
//      SLOT(onRadiusArrayChanged(pqVariableType, const QString&)));
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::representationTypeChanged()
{
  if (this->Internals)
    {
    const char* reprType = vtkSMPropertyHelper
        ( this->Internals->Representation, "Representation" ).GetAsString();
    if ( strcmp(  reprType, "Splotch particles"  ) == 0 )
      {
      this->Internals->Frame->setEnabled(true);
      vtkSMPropertyHelper(this->Internals->Representation,
        "InterpolateScalarsBeforeMapping").Set(0);
      this->Internals->Representation->UpdateVTKObjects();
      }
    else
      {
      this->Internals->Frame->setEnabled(false);
      }
    }
}
//-----------------------------------------------------------------------------
/*
void pqSplotchRaytraceDisplayPanelDecorator::onIntensityArrayChanged(
    pqVariableType type, const QString& name)
{
  if (!this->Internals->PipelineRepresentation) {
    return;
  }

  vtkSMStringVectorProperty* svp = vtkSMStringVectorProperty::SafeDownCast(
      this->Internals->PipelineRepresentation->GetProperty("IntensityArray"));
  svp->SetElement(0, 0); // idx
  svp->SetElement(1, 0); //port
  svp->SetElement(2, 0); //connection
  svp->SetElement(3, (int) vtkDataObject::FIELD_ASSOCIATION_POINTS); //type
  svp->SetElement(4, name.toAscii().data()); //name

  this->Internals->IntensityArray->reloadGUI();

//  this->Internals->PipelineRepresentation->UpdateVTKObjects();
  this->updateAllViews();
}
*/
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::updateAllViews()
{
//  if (this->Internals->PipelineRepresentation)
    {
//    this->Internals->PipelineRepresentation->renderViewEventually();
    }
}

