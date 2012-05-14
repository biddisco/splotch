/*=========================================================================

  Program:   Visualization Toolkit
  Module:    pqSplotchDisplayPanelDecorator.h

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

#ifndef __pqSplotchDisplayPanelDecorator_h
#define __pqSplotchDisplayPanelDecorator_h

#include <QGroupBox>
class pqDisplayPanel;
class pqPipelineRepresentation;
class pqWidgetRangeDomain;
class vtkSMProperty;

#include "pqVariableType.h"

class pqSplotchDisplayPanelDecorator : public QGroupBox
{
  Q_OBJECT
  typedef QGroupBox Superclass;
public:
  pqSplotchDisplayPanelDecorator(pqDisplayPanel* panel);
  ~pqSplotchDisplayPanelDecorator();

  // called when the representation changed
  void setRepresentation(pqPipelineRepresentation* repr);

  // setup the connections between the GUI and the proxies
  void setupGUIConnections();

protected slots:
  void representationTypeChanged();
  void EditColour();
  void RepaintClicked();
  void ActiveParticleTypeChanged(int v);
  void UpdateParticleTypes();

protected :

private:
  pqSplotchDisplayPanelDecorator(const pqSplotchDisplayPanelDecorator&); // Not implemented.
  void operator=(const pqSplotchDisplayPanelDecorator&); // Not implemented.

  class pqInternals;
  pqInternals* Internals;
  int TableIndex;
};

#endif


