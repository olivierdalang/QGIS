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
	mMapCanvas->viewport()->installEventFilter(this);
	mMapCanvas->viewport()->setMouseTracking(true);

  // We use the same event filter so that shortcuts are still active
  mAngleLineEdit->installEventFilter( cadDockWidget );
  mDistanceLineEdit->installEventFilter( cadDockWidget );
  mXLineEdit->installEventFilter( cadDockWidget );
  mYLineEdit->installEventFilter( cadDockWidget );

	//connect(mEnableAction, &QAction::triggered, this, &QgsAdvancedDigitizingDockWidget::activateCad);
}

bool QgsAdvancedDigitizingFloater::eventFilter(QObject *obj, QEvent *event)
{
	if (mCadDockWidget->cadEnabled() && mActive) {
		if (event->type() == QEvent::MouseMove)
		{
			QMouseEvent *mouseEvent = dynamic_cast<QMouseEvent *>(event);
			updatePos(mouseEvent->pos());
		}
		else if (event->type() == QEvent::Enter)
		{
			QEnterEvent  *enterEvent = dynamic_cast<QEnterEvent  *>(event);
			updatePos(enterEvent->pos());
			setVisible(true);
		}
		else if (event->type() == QEvent::Leave)
		{
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
	//QPoint newPos = mMapCanvas->mapToGlobal(pos) + QPoint(10, -10);
	//move(newPos);
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

void QgsAdvancedDigitizingFloater::changeLockX(bool enabled)
{
	if (!enabled) {
		mXLineEdit->setStyleSheet("");
	}
	else {
		mXLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::changeLockY(bool enabled)
{
	if (!enabled) {
		mYLineEdit->setStyleSheet("");
	}
	else {
		mYLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::changeLockDistance(bool enabled)
{
	if (!enabled) {
		mDistanceLineEdit->setStyleSheet("");
	}
	else {
		mDistanceLineEdit->setStyleSheet("font-weight: bold");
	}
}

void QgsAdvancedDigitizingFloater::changeLockAngle(bool enabled)
{
	if (!enabled) {
		mAngleLineEdit->setStyleSheet("");
	}
	else {
		mAngleLineEdit->setStyleSheet("font-weight: bold");
	}
}
