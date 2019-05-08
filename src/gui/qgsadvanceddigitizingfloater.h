/***************************************************************************
    qgsadvanceddigitizingfloater.cpp  -  floater for CAD tools
    ----------------------
    begin                : May 2019
    copyright            : (C) Olivier Dalang
    email                : olivier.dalang@gmail.com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QGSADVANCEDDIGITIZINGFLOATER
#define QGSADVANCEDDIGITIZINGFLOATER

#include <QWidget>
#include <QString>

#include "ui_qgsadvanceddigitizingfloaterbase.h"
#include "qgsadvanceddigitizingdockwidget.h"
#include "qgis_gui.h"
#include "qgis_sip.h"

class QgsMapCanvas;
class QgsAdvancedDigitizingDockWidget;

class GUI_EXPORT QgsAdvancedDigitizingFloater : public QWidget, private Ui::QgsAdvancedDigitizingFloaterBase
{
    Q_OBJECT

  public:

    explicit QgsAdvancedDigitizingFloater(QgsMapCanvas *canvas, QgsAdvancedDigitizingDockWidget *cadDockWidget);

    /**
    * getter for mActive
    */
    bool active();

  public slots:
    void setActive(bool active);
    void changeX(QString text);
    void changeY(QString text);
    void changeDistance(QString text);
    void changeAngle(QString text);
    void changeLockX(bool enabled);
    void changeLockY(bool enabled);
    void changeLockDistance(bool enabled);
    void changeLockAngle(bool enabled);

  private:

	  //! pointer to map canvas
	  QgsMapCanvas *mMapCanvas = nullptr;

	  //! pointer to map cad docker widget
	  QgsAdvancedDigitizingDockWidget *mCadDockWidget = nullptr;

	//! Whether the floater is enabled
	bool mActive = false;

	/**
	* event filter to track mouse position
	* \note defined as private in Python bindings
	*/
	bool eventFilter(QObject *obj, QEvent *event) override SIP_SKIP;

	/**
	* move the widget
    * \param pos position of the cursor
	*/
	void updatePos(QPoint pos);

};

#endif // QGSADVANCEDDIGITIZINGFLOATER_H
