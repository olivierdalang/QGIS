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

#include <QMouseEvent>
#include <QEnterEvent>
#include <QLocale>

#include "qgsadvanceddigitizingfloater.h"
#include "qgsmessagelog.h"
#include "qgsmapcanvas.h"

QgsAdvancedDigitizingFloater::QgsAdvancedDigitizingFloater(QgsMapCanvas *canvas, QgsAdvancedDigitizingDockWidget *cadDockWidget)
	: QWidget(canvas->viewport()), mMapCanvas(canvas), mCadDockWidget(cadDockWidget)
{
	setupUi(this);
	setWindowFlag(Qt::FramelessWindowHint);
	setAttribute(Qt::WA_TransparentForMouseEvents);
  setVisible(false);

  // This is required to be able to track mouse move events
	mMapCanvas->viewport()->installEventFilter(this);
	mMapCanvas->viewport()->setMouseTracking(true);

  // We reuse cadDockWidget's eventFilter for the CAD specific shortcuts
  mAngleLineEdit->installEventFilter( cadDockWidget );
  mDistanceLineEdit->installEventFilter( cadDockWidget );
  mXLineEdit->installEventFilter( cadDockWidget );
  mYLineEdit->installEventFilter( cadDockWidget );

	// Connect all cadDockWidget's signals to update the widget's display
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::valueXChanged, this, &QgsAdvancedDigitizingFloater::changeX);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::valueYChanged, this, &QgsAdvancedDigitizingFloater::changeY);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::valueAngleChanged, this, &QgsAdvancedDigitizingFloater::changeAngle);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::valueDistanceChanged, this, &QgsAdvancedDigitizingFloater::changeDistance);

  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::lockXChanged, this, &QgsAdvancedDigitizingFloater::changeLockX);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::lockYChanged, this, &QgsAdvancedDigitizingFloater::changeLockY);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::lockAngleChanged, this, &QgsAdvancedDigitizingFloater::changeLockAngle);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::lockDistanceChanged, this, &QgsAdvancedDigitizingFloater::changeLockDistance);

  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::focusOnX, this, &QgsAdvancedDigitizingFloater::focusOnX);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::focusOnY, this, &QgsAdvancedDigitizingFloater::focusOnY);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::focusOnAngle, this, &QgsAdvancedDigitizingFloater::focusOnAngle);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::focusOnDistance, this, &QgsAdvancedDigitizingFloater::focusOnDistance);

  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::enabledChangedX, this, &QgsAdvancedDigitizingFloater::enabledChangedX);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::enabledChangedY, this, &QgsAdvancedDigitizingFloater::enabledChangedY);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::enabledChangedAngle, this, &QgsAdvancedDigitizingFloater::enabledChangedAngle);
  connect(cadDockWidget, &QgsAdvancedDigitizingDockWidget::enabledChangedDistance, this, &QgsAdvancedDigitizingFloater::enabledChangedDistance);

	// Connect our line edits signals to update cadDockWidget's tate
  connect(mXLineEdit, &QLineEdit::returnPressed, cadDockWidget, [=]() { cadDockWidget->setX(mXLineEdit->text()); });
  connect(mYLineEdit, &QLineEdit::returnPressed, cadDockWidget, [=]() { cadDockWidget->setY(mYLineEdit->text()); });
  connect(mAngleLineEdit, &QLineEdit::returnPressed, cadDockWidget, [=]() { cadDockWidget->setAngle(mAngleLineEdit->text()); });
  connect(mDistanceLineEdit, &QLineEdit::returnPressed, cadDockWidget, [=]() { cadDockWidget->setDistance(mDistanceLineEdit->text()); });
}

bool QgsAdvancedDigitizingFloater::eventFilter(QObject *obj, QEvent *event)
{
	if (mCadDockWidget->cadEnabled() && mActive) {
		if (event->type() == QEvent::MouseMove)
		{
      // We update the position when mouse moves
			QMouseEvent *mouseEvent = dynamic_cast<QMouseEvent *>(event);
			updatePos(mouseEvent->pos());
		}
		else if (event->type() == QEvent::Enter)
		{
      // We show the widget when mouse enters
			QEnterEvent  *enterEvent = dynamic_cast<QEnterEvent  *>(event);
			updatePos(enterEvent->pos());
			setVisible(true);
		}
		else if (event->type() == QEvent::Leave)
		{
      // We hide the widget when mouse leaves
			setVisible(false);
		}
	}
	return QWidget::eventFilter(obj, event);
}

bool QgsAdvancedDigitizingFloater::active()
{
	return mActive;
}

void QgsAdvancedDigitizingFloater::setActive(bool active)
{
	mActive = active;
	if (!active) {
		setVisible(false);
	}
}

void QgsAdvancedDigitizingFloater::updatePos(QPoint pos)
{
  // We hardcode a small delta between the mouse position and the widget's position
	move(pos + QPoint(15,5));
}

void QgsAdvancedDigitizingFloater::changeX(QString text)
{
	mXLineEdit->setText(text);
}

void QgsAdvancedDigitizingFloater::changeY(QString text)
{
	mYLineEdit->setText(text);
}

void QgsAdvancedDigitizingFloater::changeDistance(QString text)
{
	mDistanceLineEdit->setText(text);
}

void QgsAdvancedDigitizingFloater::changeAngle(QString text)
{
	mAngleLineEdit->setText(text);
}

void QgsAdvancedDigitizingFloater::changeLockX(bool locked)
{
	if (!locked) {
		mXLineEdit->setStyleSheet("");
	}
	else {
		mXLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::changeLockY(bool locked)
{
	if (!locked) {
		mYLineEdit->setStyleSheet("");
	}
	else {
		mYLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::changeLockDistance(bool locked)
{
	if (!locked) {
		mDistanceLineEdit->setStyleSheet("");
	}
	else {
		mDistanceLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::changeLockAngle(bool locked)
{
	if (!locked) {
		mAngleLineEdit->setStyleSheet("");
	}
	else {
		mAngleLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::focusOnX()
{
	if(mActive){
	mXLineEdit->setFocus();
    mXLineEdit->selectAll();
	}
}

void QgsAdvancedDigitizingFloater::focusOnY()
{
	if (mActive) {
		mYLineEdit->setFocus();
    mYLineEdit->selectAll();
}
}

void QgsAdvancedDigitizingFloater::focusOnDistance()
{
	if (mActive) {
		mDistanceLineEdit->setFocus();
    mDistanceLineEdit->selectAll();
}
}

void QgsAdvancedDigitizingFloater::focusOnAngle()
{
	if (mActive) {
		mAngleLineEdit->setFocus();
    mAngleLineEdit->selectAll();
  }
}


void QgsAdvancedDigitizingFloater::enabledChangedX(bool enabled)
{
	mXLineEdit->setVisible(enabled);
	mXLabel->setVisible(enabled);
}

void QgsAdvancedDigitizingFloater::enabledChangedY(bool enabled)
{
	mYLineEdit->setVisible(enabled);
	mYLabel->setVisible(enabled);
}

void QgsAdvancedDigitizingFloater::enabledChangedDistance(bool enabled)
{
	mDistanceLineEdit->setVisible(enabled);
	mDistanceLabel->setVisible(enabled);
}

void QgsAdvancedDigitizingFloater::enabledChangedAngle(bool enabled)
{
	mAngleLineEdit->setVisible(enabled);
	mAngleLabel->setVisible(enabled);
}
