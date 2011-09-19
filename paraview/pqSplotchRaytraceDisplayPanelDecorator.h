/*=========================================================================

  Program:   Visualization Toolkit
  Module:    pqSplotchRaytraceDisplayPanelDecorator.h

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

#ifndef __pqSplotchRaytraceDisplayPanelDecorator_h
#define __pqSplotchRaytraceDisplayPanelDecorator_h

#include <QGroupBox>
class pqDisplayPanel;
class pqPipelineRepresentation;
class pqWidgetRangeDomain;
class vtkSMProperty;

#include "pqVariableType.h"

class pqSplotchRaytraceDisplayPanelDecorator : public QGroupBox
{
  Q_OBJECT
  typedef QGroupBox Superclass;
public:
  pqSplotchRaytraceDisplayPanelDecorator(pqDisplayPanel* panel);
  ~pqSplotchRaytraceDisplayPanelDecorator();

  // called when the representation changed
  void setRepresentation(pqPipelineRepresentation* repr);

  // setup the connections between the GUI and the proxies
  void setupGUIConnections();

  void updateAllViews();

protected slots:
  void representationTypeChanged();
  void editColour0();

protected :

private:
  pqSplotchRaytraceDisplayPanelDecorator(const pqSplotchRaytraceDisplayPanelDecorator&); // Not implemented.
  void operator=(const pqSplotchRaytraceDisplayPanelDecorator&); // Not implemented.

  class pqInternals;
  pqInternals* Internals;
};

#endif


