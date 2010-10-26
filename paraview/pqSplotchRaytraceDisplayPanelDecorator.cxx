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
//  Copyright (c) CSCS - Swiss National Supercomputing Centre
//                EDF - Electricite de France
//
//  John Biddiscombe, Ugo Varetto (CSCS)
//  Stephane Ploix (EDF)
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

class pqSplotchRaytraceDisplayPanelDecorator::pqInternals: public Ui::pqSplotchRaytraceDisplayPanelDecorator
{
public:
  pqPropertyLinks Links;
  vtkSMPVRepresentationProxy* RepresentationProxy;
  vtkSmartPointer<vtkEventQtSlotConnect> VTKConnect;
  pqPipelineRepresentation* PipelineRepresentation;

  pqWidgetRangeDomain* MaxPixelSizeRangeDomain;
  pqWidgetRangeDomain* OpacityRangeDomain;
  pqWidgetRangeDomain* RadiusRangeDomain;

  pqInternals(QWidget* parent)
  {
/*
    this->RepresentationProxy = 0;
    this->VTKConnect = vtkSmartPointer<vtkEventQtSlotConnect>::New();
    this->TransferFunctionDialog = new pqTransferFunctionDialog(parent);
    this->MaxPixelSizeRangeDomain = NULL;
    this->OpacityRangeDomain = NULL;
    this->RadiusRangeDomain = NULL;
*/
  }
};

//-----------------------------------------------------------------------------
pqSplotchRaytraceDisplayPanelDecorator::pqSplotchRaytraceDisplayPanelDecorator(
    pqDisplayPanel* disp_panel) :
  Superclass(disp_panel)
{
  pqDisplayProxyEditor* panel =
      qobject_cast<pqDisplayProxyEditor*> (disp_panel);
  pqRepresentation* repr = panel->getRepresentation();
  vtkSMProxy* reprProxy = (repr) ? repr->getProxy() : NULL;
  // vtkSMProperty* prop;
  this->Internals = NULL;
/*
  if (!reprProxy || 
    !reprProxy->GetXMLName() ||
    (strcmp(reprProxy->GetXMLName(), "GeometryRepresentation") != 0  &&
     strcmp(reprProxy->GetXMLName(), "UnstructuredGridRepresentation") != 0  &&
     strcmp(reprProxy->GetXMLName(), "UniformGridRepresentation") != 0))
    {
    return;
    }

  if (!pqSMAdaptor::getEnumerationPropertyDomain(
      reprProxy->GetProperty("Representation")).contains("Point Sprite"))
    {
    return;
    }

  // This is not advisable, but we do it nonetheless since that's what the old
  // code was doing. At some point we need to clean this up.
  vtkSMSplotchRaytraceRepresentationProxy::InitializeDefaultValues(reprProxy);
*/
  
  this->Internals = new pqInternals(this);
  QVBoxLayout* vlayout = dynamic_cast<QVBoxLayout*> (panel->layout());
  if (vlayout)
    {
    vlayout->insertWidget(2, this);
    }
  else
    {
    panel->layout()->addWidget(this);
    }
  this->Internals->setupUi(this);
  this->Internals->RepresentationProxy = vtkSMPVRepresentationProxy::SafeDownCast(reprProxy);
/*
  // setup the scaleBy and radiusBy menus
  this->Internals->ScaleBy->setConstantVariableName("Constant Radius");
  this->Internals->ScaleBy->setPropertyArrayName("RadiusArray");
  this->Internals->ScaleBy->setPropertyArrayComponent("RadiusVectorComponent");
  this->Internals->ScaleBy->setToolTip(
    "select method for scaling the point sprites.");

  this->Internals->OpacityBy->setConstantVariableName("Constant Opacity");
  this->Internals->OpacityBy->setPropertyArrayName("OpacityArray");
  this->Internals->OpacityBy->setPropertyArrayComponent(
    "OpacityVectorComponent");
  this->Internals->OpacityBy->setToolTip(
    "select method for setting the opacity of the point sprites.");

  this->Internals->ScaleBy->reloadGUI();
  this->Internals->OpacityBy->reloadGUI();

  this->setupGUIConnections();

  this->setRepresentation(
    static_cast<pqPipelineRepresentation*> (panel->getRepresentation()));
  QObject::connect(&this->Internals->Links, SIGNAL(smPropertyChanged()), panel,
      SLOT(updateAllViews()), Qt::QueuedConnection);

  this->connect(this->Internals->OpacityMapping, SIGNAL(clicked()), this,
      SLOT(showOpacityDialog()));

  this->connect(this->Internals->RadiusMapping, SIGNAL(clicked()), this,
      SLOT(showRadiusDialog()));

  this->Internals->TransferFunctionDialog->setRepresentation(
      static_cast<pqPipelineRepresentation*> (panel->getRepresentation()));

  this->reloadGUI();
*/
}

//-----------------------------------------------------------------------------
pqSplotchRaytraceDisplayPanelDecorator::~pqSplotchRaytraceDisplayPanelDecorator()
{
  delete this->Internals;
  this->Internals = 0;
}
