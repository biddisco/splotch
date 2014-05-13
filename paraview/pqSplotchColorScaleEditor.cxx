/*=========================================================================

   Program: ParaView
   Module:    pqSplotchColorScaleEditor.cxx

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2. 

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

/// \file pqSplotchColorScaleEditor.cxx
/// \date 2/14/2007

#include "pqSplotchColorScaleEditor.h"

#include "pqApplicationCore.h"
#include "pqChartValue.h"
#include "pqColorMapModel.h"
#include "pqColorPresetManager.h"
#include "pqColorPresetModel.h"
#include "pqCoreUtilities.h"
#include "pqDataRepresentation.h"
#include "pqLookupTableManager.h"
#include "pqOutputPort.h"
#include "pqPipelineRepresentation.h"
#include "pqPropertyLinks.h"
#include "pqRenderViewBase.h"
#include "pqRescaleRange.h"
#include "pqScalarBarRepresentation.h"
#include "pqScalarOpacityFunction.h"
#include "pqScalarsToColors.h"
#ifdef FIXME
#include "pqScatterPlotRepresentation.h"
#endif
#include "pqSignalAdaptors.h"
#include "pqSMAdaptor.h"
#include "pqStandardColorLinkAdaptor.h"
#include "vtkColorTransferFunction.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPVTemporalDataInformation.h"
#include "vtkSMProperty.h"
#include "vtkSMProxy.h"
#include "vtkSMPVRepresentationProxy.h"
#include "vtkTransferFunctionViewer.h"
#include "vtkType.h"

#include <QCloseEvent>
#include <QColor>
#include <QColorDialog>
#include <QDoubleValidator>
#include <QGridLayout>
#include <QIntValidator>
#include <QItemSelectionModel>
#include <QList>
#include <QMenu>
#include <QMessageBox>
#include <QPointer>
#include <QSpacerItem>
#include <QString>
#include <QtDebug>
#include <QVariant>


//----------------------------------------------------------------------------
pqSplotchColorScaleEditor::pqSplotchColorScaleEditor(QWidget *widgetParent)
  : pqColorMapEditor(widgetParent)
{
  this->ActiveColorTable = NULL;
}
//----------------------------------------------------------------------------
pqSplotchColorScaleEditor::~pqSplotchColorScaleEditor()
{
}
//----------------------------------------------------------------------------
pqScalarsToColors *pqSplotchColorScaleEditor::GetColorTableToEdit()
{
  return this->ActiveColorTable;
}
//----------------------------------------------------------------------------
void pqSplotchColorScaleEditor::setActiveColorTable(pqScalarsToColors *table)
{
  this->ActiveColorTable = table;
}
//----------------------------------------------------------------------------

