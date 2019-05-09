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

/**
* \ingroup gui
* \brief The QgsAdvancedDigitizingFloater class is widget that floats
* next to the mouse pointer, and allow interaction with the AdvancedDigitizing
* feature. It proxies display and actions to QgsMapToolAdvancedDigitizingDockWidget.
*/
class GUI_EXPORT QgsAdvancedDigitizingFloater : public QWidget, private Ui::QgsAdvancedDigitizingFloaterBase
{
    Q_OBJECT

  public:

    explicit QgsAdvancedDigitizingFloater(QgsMapCanvas *canvas, QgsAdvancedDigitizingDockWidget *cadDockWidget);

    /**
    * Whether the floater is active or not.
    * Note that the floater may be active but not visible (e.g. if the mouse is not over the canvas).
    */
    bool active();

  public slots:
	  /**
	  * Set whether the floater should be active or not.
    * Note that the floater may be active but not visible (e.g. if the mouse is not over the canvas).
	  *
	  * \param active
	  */
    void setActive(bool active);

private slots:

	void changeX(QString text);
	void changeY(QString text);
	void changeDistance(QString text);
	void changeAngle(QString text);
	void changeLockX(bool locked);
	void changeLockY(bool locked);
	void changeLockDistance(bool locked);
	void changeLockAngle(bool locked);
	void focusOnX();
	void focusOnY();
	void focusOnAngle();
	void focusOnDistance();
	void enabledChangedX(bool enabled);
	void enabledChangedY(bool enabled);
	void enabledChangedAngle(bool enabled);
	void enabledChangedDistance(bool enabled);

  private:

	  //! pointer to map canvas
	  QgsMapCanvas *mMapCanvas = nullptr;

	  //! pointer to map cad docker widget
	  QgsAdvancedDigitizingDockWidget *mCadDockWidget = nullptr;

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

	//! Whether the floater is enabled.
	bool mActive = false;

};

#endif // QGSADVANCEDDIGITIZINGFLOATER_H
