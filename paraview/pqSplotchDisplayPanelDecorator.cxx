/*=========================================================================

 Program:   Visualization Toolkit
 Module:    pqSplotchDisplayPanelDecorator.cxx

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/

// .NAME pqSplotchDisplayPanelDecorator
// .SECTION Thanks
// <verbatim>
//
//  This file is part of the Splotchs plugin developed and contributed by
//
// </verbatim>

#include "pqSplotchDisplayPanelDecorator.h"
#include "ui_pqSplotchDisplayPanelDecorator.h"

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
#include "vtkPVDataSetAttributesInformation.h"
#include "vtkPVArrayInformation.h"

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
#include "pqComboBoxDomain.h"

// splotch enhanced qt classes
//#include "pqSplotchColorScaleEditor.h"

class pqSplotchDisplayPanelDecorator::pqInternals: public Ui::pqSplotchDisplayPanelDecorator
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
    this->PipelineRepresentation = 0;
  }
};

//-----------------------------------------------------------------------------
pqSplotchDisplayPanelDecorator::pqSplotchDisplayPanelDecorator(
  pqDisplayPanel* _panel) : Superclass(_panel)
{
  pqDisplayProxyEditor *panel = qobject_cast<pqDisplayProxyEditor*> (_panel);
  pqRepresentation     *repr  = panel->getRepresentation();
  vtkSMProxy      *reprProxy  = (repr) ? repr->getProxy() : NULL;
  this->Internals             = NULL;

  QGroupBox* wid = new QGroupBox(panel);
  this->Internals = new pqInternals(this);
  this->Internals->Frame = wid;
  this->Internals->setupUi(wid);
  QVBoxLayout* l = qobject_cast<QVBoxLayout*>(panel->layout());
  l->addWidget(wid);
  //
  this->setRepresentation(
    static_cast<pqPipelineRepresentation*> (panel->getRepresentation()));
  //
  this->setupGUIConnections();
  //
  if (!this->Internals->RepresentationProxy) return;
  this->UpdateParticleTypes();
}
//-----------------------------------------------------------------------------
pqSplotchDisplayPanelDecorator::~pqSplotchDisplayPanelDecorator()
{
  delete this->Internals;
  this->Internals = 0;
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::setupGUIConnections()
{
  QObject::connect(this->Internals->EditColorMapButton, SIGNAL(clicked()), this,
      SLOT(EditColour()), Qt::QueuedConnection);

  QObject::connect(this->Internals->RepaintButton, SIGNAL(clicked()), this,
      SLOT(RepaintClicked()), Qt::QueuedConnection);

  QObject::connect(this->Internals->ActiveParticleType, SIGNAL(valueChanged(int)), this,
      SLOT(ActiveParticleTypeChanged(int)), Qt::QueuedConnection);

  QObject::connect(this->Internals->TypeArray, SIGNAL(currentIndexChanged(int)), this,
      SLOT(UpdateParticleTypes()), Qt::QueuedConnection);
  
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::setRepresentation(pqPipelineRepresentation* repr)
{
  if (this->Internals->PipelineRepresentation == repr) {
    return;
  }

  if (this->Internals->PipelineRepresentation) {
    // break all old links.
    this->Internals->Links.removeAllPropertyLinks();
  }

  this->Internals->PipelineRepresentation = repr;
  if (!repr) {
//    this->Internals->TransferFunctionDialog->hide();
    return;
  }

  vtkSMProperty* prop;
  vtkSMProxy *reprProxy  = repr->getProxy();
  this->Internals->RepresentationProxy = vtkSMPVRepresentationProxy::SafeDownCast(reprProxy);
  reprProxy->GetProperty("Input")->UpdateDependentDomains();

  //
  // Field array controls
  //
  // Intensity
  //
  prop = reprProxy->GetProperty("IntensityScalars");

  // Some representations are not compatible with our mapper, jump out before segfaults occur
  //
  if (prop) {
    // adaptor from combo to property
    pqSignalAdaptorComboBox *adaptor = new pqSignalAdaptorComboBox(this->Internals->IntensityArray);
    // domain to control the combo contents
    new pqComboBoxDomain(this->Internals->IntensityArray, prop, "array_list");
    // link gui changes to property and vice versa
    this->Internals->Links.addPropertyLink(adaptor, "currentText", SIGNAL(currentTextChanged(const QString&)), reprProxy, prop);
    prop->UpdateDependentDomains();

    //
    // Radius
    //
    prop = reprProxy->GetProperty("RadiusScalars");
    // adaptor from combo to property
    adaptor = new pqSignalAdaptorComboBox(this->Internals->RadiusArray);
    // domain to control the combo contents
    new pqComboBoxDomain(this->Internals->RadiusArray, prop, "array_list");
    // link gui changes to property and vice versa
    this->Internals->Links.addPropertyLink(adaptor, "currentText", SIGNAL(currentTextChanged(const QString&)), reprProxy, prop);
    prop->UpdateDependentDomains();

    //
    // Type
    //
    prop = reprProxy->GetProperty("TypeScalars");
    // adaptor from combo to property
    adaptor = new pqSignalAdaptorComboBox(this->Internals->TypeArray);
    // domain to control the combo contents
    new pqComboBoxDomain(this->Internals->TypeArray, prop, "array_list");
    // link gui changes to property and vice versa
    this->Internals->Links.addPropertyLink(adaptor, "currentText", SIGNAL(currentTextChanged(const QString&)), reprProxy, prop);
    prop->UpdateDependentDomains();

    //
    // Active
    //
    prop = reprProxy->GetProperty("ActiveScalars");
    // adaptor from combo to property
    adaptor = new pqSignalAdaptorComboBox(this->Internals->ActiveArray);
    // domain to control the combo contents
    new pqComboBoxDomain(this->Internals->ActiveArray, prop, "array_list");
    // link gui changes to property and vice versa
    this->Internals->Links.addPropertyLink(adaptor, "currentText", SIGNAL(currentTextChanged(const QString&)), reprProxy, prop);
    prop->UpdateDependentDomains();

    //
    // Simple controls
    //
    this->Internals->Links.addPropertyLink(
      this->Internals->ActiveParticleType, "value", SIGNAL(valueChanged(int)),
      reprProxy, reprProxy->GetProperty("ActiveParticleType"));

    this->Internals->Links.addPropertyLink(
      this->Internals->Brightness, "value", SIGNAL(valueChanged(int)),
      reprProxy, reprProxy->GetProperty("Brightness"));
    
    this->Internals->Links.addPropertyLink(
      this->Internals->BrightnessLOD, "value", SIGNAL(valueChanged(int)),
      reprProxy, reprProxy->GetProperty("BrightnessLOD"));
    
    this->Internals->Links.addPropertyLink(
      this->Internals->RadiusMultiplier, "value", SIGNAL(valueChanged(int)),
      reprProxy, reprProxy->GetProperty("RadiusMultiplier"));
    
    this->Internals->Links.addPropertyLink(
      this->Internals->LogIntensity, "checked", SIGNAL(toggled(bool)),
      reprProxy, reprProxy->GetProperty("LogIntensity"));
    
    this->Internals->Links.addPropertyLink(
      this->Internals->TypeActive, "checked", SIGNAL(toggled(bool)),
      reprProxy, reprProxy->GetProperty("TypeActive"));

    this->Internals->Links.addPropertyLink(
      this->Internals->Gray, "value", SIGNAL(valueChanged(int)),
      reprProxy, reprProxy->GetProperty("GrayAbsorption"));
    
    // Connect Property event to GUI
    this->Internals->VTKConnect->Connect(
      reprProxy->GetProperty("Representation"),
      vtkCommand::ModifiedEvent, this, SLOT(representationTypeChanged()));
  }
  this->representationTypeChanged();
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::representationTypeChanged()
{
  if (this->Internals && this->Internals->RepresentationProxy) {
    const char* reprType = vtkSMPropertyHelper
        ( this->Internals->RepresentationProxy, "Representation" ).GetAsString();
    if ( strcmp(  reprType, "Splotch particles"  ) == 0 ) {
      this->Internals->Frame->setEnabled(true);
      vtkSMPropertyHelper(this->Internals->RepresentationProxy,
        "InterpolateScalarsBeforeMapping").Set(0);
      this->Internals->RepresentationProxy->UpdateVTKObjects();
    }
    else {
      this->Internals->Frame->setEnabled(true);
    }
  }
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::UpdateParticleTypes()
{
  vtkPVDataInformation *dataInfo = 
    this->Internals->PipelineRepresentation->getInputDataInformation();
  vtkPVDataInformation* geomInfo = 
    this->Internals->RepresentationProxy->GetRepresentedDataInformation();
  vtkPVDataSetAttributesInformation *pointInfo = 
    dataInfo->GetPointDataInformation();
  vtkPVArrayInformation *arrayInfo = pointInfo->GetArrayInformation(
    this->Internals->TypeArray->currentText().toAscii().data());
  if (!arrayInfo) return;
  //
  QString valstr;
  double *range = arrayInfo->GetComponentRange(0);
  int ntypes = 1+static_cast<int>(range[1]);
  if (ntypes>9) {
    ntypes = 9;
    valstr = "(Error) Clamped to 9";
  }
  else {
    valstr = QString::number(ntypes);
  }
  this->Internals->ActiveParticleType->setMaximum(ntypes-1);
  this->Internals->typeslabel->setText(valstr);
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::ActiveParticleTypeChanged(int v)
{
  this->Internals->TableIndex = v;
  //
  vtkSMProperty* SettingsProperty = this->Internals->RepresentationProxy->GetProperty("ActiveParticleSettings");
  this->Internals->RepresentationProxy->UpdatePropertyInformation(SettingsProperty);
  QList<QVariant> ActiveParticleSettings = pqSMAdaptor::getMultipleElementProperty(SettingsProperty);
  //
  int ptype = ActiveParticleSettings.at(0).toString().toInt();
  if (ptype==this->Internals->ActiveParticleType->value()) {
    double brightness = ActiveParticleSettings.at(1).toString().toDouble();
    this->Internals->Brightness->setValue(log10(brightness)*100);
    //
    bool logI = ActiveParticleSettings.at(2).toString().toInt();
    this->Internals->LogIntensity->setChecked(logI);
    //
    int t = this->Internals->IntensityArray->findText(ActiveParticleSettings.at(3).toString());
    if (t==-1) { t=0; }
    this->Internals->IntensityArray->setCurrentIndex(t);
    //
    t = this->Internals->RadiusArray->findText(ActiveParticleSettings.at(4).toString());
    if (t==-1) { t=0; }
    this->Internals->RadiusArray->setCurrentIndex(t);
    //
    bool active = ActiveParticleSettings.at(5).toString().toInt();
    this->Internals->TypeActive->setChecked(active);
  }  
  for (int i=0; i<ActiveParticleSettings.size(); i++) {
    std::cout << ActiveParticleSettings.at(i).toString().toAscii().data() << std::endl;
  }
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::EditColour()
{
  pqApplicationCore       *core = pqApplicationCore::instance();
  pqLookupTableManager *lut_mgr = core->getLookupTableManager();
  pqScalarsToColors      *pqlut = this->Internals->PipelineRepresentation->getLookupTable();
  vtkSMProxy               *lut = (pqlut)? pqlut->getProxy() : 0;

 // pqSplotchColorScaleEditor editor(pqCoreUtilities::mainWidget());
 // editor.setActiveColorTable(pqlut);
//  editor.setRepresentation(this->Internals->PipelineRepresentation);
//  editor.exec();
}
//-----------------------------------------------------------------------------
void pqSplotchDisplayPanelDecorator::RepaintClicked()
{
  if (this->Internals->PipelineRepresentation) {
    this->Internals->PipelineRepresentation->renderViewEventually();
  }
}
