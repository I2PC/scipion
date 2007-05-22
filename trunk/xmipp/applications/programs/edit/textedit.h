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

#include <qmainwindow.h>
#include <qmap.h>
#include <qaccel.h>

#include <vector>

class QAction;
class QComboBox;
class QTabWidget;
class QTextEdit;

class TextEdit : public QMainWindow
{
    Q_OBJECT

public:
    TextEdit(QWidget *parent = 0, const char *name = 0);

public slots:
    QTextEdit *currentEditor() const;
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
    void doConnections(QTextEdit *e);

private slots:
    void editorChanged(QWidget *);
    void execAccel();

private:
    QTabWidget *tabWidget;
    QMap<QTextEdit*, QString> filenames;
    QMap<QAccel*, QString> accels;
    bool remove;
    std::vector< QString > files_open;
};
#endif
