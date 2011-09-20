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
#include "pqLookupTableManager.h"
#include "pqApplicationCore.h"
#include "pqCoreUtilities.h"

// splotch enhanced qt classes
#include "pqSplotchColorScaleEditor.h"

class pqSplotchRaytraceDisplayPanelDecorator::pqInternals: public Ui::pqSplotchRaytraceDisplayPanelDecorator
{
public:
  vtkSmartPointer<vtkEventQtSlotConnect> VTKConnect;
  pqPropertyLinks                        Links;
  vtkSMPVRepresentationProxy            *RepresentationProxy;
  pqPipelineRepresentation              *PipelineRepresentation;
  QWidget                               *Frame;
  QList<pqScalarsToColors*>              ColorTableList;
  int                                    TableIndex;
  //
  pqInternals(QWidget* parent)
  {
    this->VTKConnect     = vtkSmartPointer<vtkEventQtSlotConnect>::New();
    this->Frame          = 0;
    this->TableIndex     = 0;
    this->RepresentationProxy = 0;
  }
};

//-----------------------------------------------------------------------------
pqSplotchRaytraceDisplayPanelDecorator::pqSplotchRaytraceDisplayPanelDecorator(
  pqDisplayPanel* _panel):Superclass(_panel)
{
  pqDisplayProxyEditor *panel = qobject_cast<pqDisplayProxyEditor*> (_panel);
  pqRepresentation     *repr  = panel->getRepresentation();
  vtkSMProxy      *reprProxy  = (repr) ? repr->getProxy() : NULL;
  this->Internals             = NULL;

  //
  // If the representation doesn't have this property, then it's not our splotch representation
  //
  vtkSMProperty* prop = reprProxy->GetProperty("IntensityScalars");
  if (!prop)  {
    return;
  }

  QWidget* wid = new QWidget(panel);
  this->Internals = new pqInternals(this);
  this->Internals->Frame = wid;
  this->Internals->setupUi(wid);
  QVBoxLayout* l = qobject_cast<QVBoxLayout*>(panel->layout());
  l->addWidget(wid);
  //
  this->Internals->RepresentationProxy = vtkSMPVRepresentationProxy::SafeDownCast(reprProxy);
  this->Internals->PipelineRepresentation = qobject_cast<pqPipelineRepresentation*>(repr);

  //
  // Intensity
  //
  prop = reprProxy->GetProperty("IntensityScalars");
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
  // Colour scalars display control
  //
//  this->Internals->ColorBy->setPropertyArrayName("ColorArrayName");
//  this->Internals->ColorBy->setPropertyArrayComponent("ColorAttributeType");
//  this->Internals->ColorBy->setRepresentation(this->Internals->PipelineRepresentation);
  //
  // 
  //
  this->Internals->VTKConnect->Connect(
      this->Internals->RepresentationProxy->GetProperty("Representation"),
      vtkCommand::ModifiedEvent, this, SLOT(representationTypeChanged()));

  this->Internals->Links.addPropertyLink(
    this->Internals->Brightness, "value", SIGNAL(valueChanged(int)),
    reprProxy, reprProxy->GetProperty("Brightness"));
  
  this->Internals->Links.addPropertyLink(
    this->Internals->Gray, "value", SIGNAL(valueChanged(int)),
    reprProxy, reprProxy->GetProperty("GrayAbsorption"));
  
  this->Internals->Links.addPropertyLink(
    this->Internals->LogIntensity, "checked", SIGNAL(toggled(bool)),
    reprProxy, reprProxy->GetProperty("LogIntensity"));
  
  //
  //
  //
  this->setupGUIConnections();
}
//-----------------------------------------------------------------------------
pqSplotchRaytraceDisplayPanelDecorator::~pqSplotchRaytraceDisplayPanelDecorator()
{
  delete this->Internals;
  this->Internals = 0;
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::setupGUIConnections()
{
  QObject::connect(this->Internals->EditColorMapButton, SIGNAL(clicked()), this,
      SLOT(EditColour()), Qt::QueuedConnection);

  QObject::connect(this->Internals->RepaintButton, SIGNAL(clicked()), this,
      SLOT(RepaintClicked()), Qt::QueuedConnection);

  QObject::connect(this->Internals->TypeSpin, SIGNAL(valueChanged(int)), this,
      SLOT(TypeSpinChanged(int)), Qt::QueuedConnection);

}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::setRepresentation(
    pqPipelineRepresentation* repr)
{
  this->Internals->PipelineRepresentation = repr;
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::representationTypeChanged()
{
  if (this->Internals) {
    const char* reprType = vtkSMPropertyHelper
        ( this->Internals->RepresentationProxy, "Representation" ).GetAsString();
    if ( strcmp(  reprType, "Splotch particles"  ) == 0 ) {
      this->Internals->Frame->setEnabled(true);
      vtkSMPropertyHelper(this->Internals->RepresentationProxy,
        "InterpolateScalarsBeforeMapping").Set(0);
      this->Internals->RepresentationProxy->UpdateVTKObjects();
    }
    else {
      this->Internals->Frame->setEnabled(false);
    }
  }
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::TypeSpinChanged(int v)
{
  this->Internals->TableIndex = v;
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::EditColour()
{
  pqApplicationCore       *core = pqApplicationCore::instance();
  pqLookupTableManager *lut_mgr = core->getLookupTableManager();
  pqScalarsToColors      *pqlut = this->Internals->PipelineRepresentation->getLookupTable();
  vtkSMProxy               *lut = (pqlut)? pqlut->getProxy() : 0;

  pqSplotchColorScaleEditor editor(pqCoreUtilities::mainWidget());
  editor.setActiveColorTable(pqlut);
  editor.setRepresentation(this->Internals->PipelineRepresentation);
  editor.exec();
}
//-----------------------------------------------------------------------------
void pqSplotchRaytraceDisplayPanelDecorator::RepaintClicked()
{
  if (this->Internals->PipelineRepresentation) {
    this->Internals->PipelineRepresentation->renderViewEventually();
  }
}
