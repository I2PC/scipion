/****************************************************************************
** $Id$
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#ifndef TEXTEDIT_H
#define TEXTEDIT_H

#ifdef QT3_SUPPORT
#include <q3mainwindow.h>
#include <q3accel.h>
#else
#include <qmainwindow.h>
#include <qaccel.h>
#endif

#include <qmap.h>

#include <vector>

class QAction;
class QComboBox;
class QTabWidget;

#ifdef QT3_SUPPORT
// MOC_SKIP_BEGIN
class Q3TextEdit;

class TextEdit : public Q3MainWindow

#else
// MOC_SKIP_END
class QTextEdit;

class TextEdit : public QMainWindow
#endif

{
    Q_OBJECT

public:
    TextEdit(QWidget *parent = 0, const char *name = 0);

public slots:
#ifdef QT3_SUPPORT
    Q3TextEdit *currentEditor() const;
#else
    QTextEdit* currentEditor() const;
#endif
    void load(const QString &f);

    void fileNew();
    void fileOpen();
    void fileSave();
    void fileSaveAs();
    void filePrint();
    void fileClose();
    void fileExit();

    void editUndo();
    void editRedo();
    void editCut();
    void editCopy();
    void editPaste();

    void addAccel(int accel, const QString &function);

    void setRemoveFlag()
    {
        remove = true;
    }
private:
    void setupFileActions();
    void setupEditActions();
    void setupTextActions();
#ifdef QT3_SUPPORT
    void doConnections(Q3TextEdit *e);
#else
    void doConnections(QTextEdit* e);
#endif

private slots:
    void editorChanged(QWidget *);
    void execAccel();

private:
    QTabWidget *tabWidget;
#ifdef QT3_SUPPORT
    QMap<Q3TextEdit*, QString> filenames;
    QMap<Q3Accel*, QString> accels;
#else
    QMap< QTextEdit*, QString > filenames;
    QMap< QAccel*, QString > accels;
#endif
    bool remove;
    std::vector< QString > files_open;
};
#endif
