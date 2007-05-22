/****************************************************************************
** $Id$
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#include "textedit.h"

#include <iostream>

#include <qtextedit.h>
#include <qaction.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qtoolbar.h>
#include <qtabwidget.h>
#include <qapplication.h>
#include <qfontdatabase.h>
#include <qcombobox.h>
#include <qlineedit.h>
#include <qfileinfo.h>
#include <qfile.h>
#include <qfiledialog.h>
#include <qprinter.h>
#include <qpaintdevicemetrics.h>
#include <qsimplerichtext.h>
#include <qcolordialog.h>
#include <qpainter.h>
#include <qinputdialog.h>

using namespace std;

TextEdit::TextEdit(QWidget *parent, const char *name)
        : QMainWindow(parent, name)
{
    setupFileActions();
    setupEditActions();
    setupTextActions();

    QPopupMenu *menu = new QPopupMenu(this);

    tabWidget = new QTabWidget(this);
    connect(tabWidget, SIGNAL(currentChanged(QWidget *)),
            this, SLOT(editorChanged(QWidget *)));
    setCentralWidget(tabWidget);

    remove = false;
}

void TextEdit::setupFileActions()
{
    QPopupMenu *menu = new QPopupMenu(this);
    menuBar()->insertItem(tr("&File"), menu);

    QAction *a;
    a = new QAction(tr("New"), tr("&New..."), CTRL + Key_N, this, "actionFileNew");
    connect(a, SIGNAL(activated()), this, SLOT(fileNew()));
    a->addTo(menu);
    a = new QAction(tr("Open"), tr("&Open..."), CTRL + Key_O, this, "actionFileOpen");
    connect(a, SIGNAL(activated()), this, SLOT(fileOpen()));
    a->addTo(menu);
    menu->insertSeparator();
    a = new QAction(tr("Save"), tr("&Save..."), CTRL + Key_S, this, "actionFileSave");
    connect(a, SIGNAL(activated()), this, SLOT(fileSave()));
    a->addTo(menu);
    a = new QAction(tr("Save As"), QPixmap(), tr("Save &As..."), 0, this, "actionFileSaveAs");
    connect(a, SIGNAL(activated()), this, SLOT(fileSaveAs()));
    a->addTo(menu);
    menu->insertSeparator();
    a = new QAction(tr("Print"), tr("&Print..."), CTRL + Key_P, this, "actionFilePrint");
    connect(a, SIGNAL(activated()), this, SLOT(filePrint()));
    a->addTo(menu);
    a = new QAction(tr("Close"), QPixmap(), tr("&Close"), CTRL + Key_W, this, "actionFileClose");
    connect(a, SIGNAL(activated()), this, SLOT(fileClose()));
    a->addTo(menu);
    a = new QAction(tr("Exit"), QPixmap(), tr("E&xit"), CTRL + Key_Q, this, "actionFileExit");
    connect(a, SIGNAL(activated()), this, SLOT(fileExit()));
    a->addTo(menu);
}

void TextEdit::setupEditActions()
{
    QPopupMenu *menu = new QPopupMenu(this);
    menuBar()->insertItem(tr("&Edit"), menu);

    QAction *a;
    a = new QAction(tr("Undo"), tr("&Undo"), CTRL + Key_Z, this, "actionEditUndo");
    connect(a, SIGNAL(activated()), this, SLOT(editUndo()));
    a->addTo(menu);
    a = new QAction(tr("Redo"), tr("&Redo"), CTRL + Key_Y, this, "actionEditRedo");
    connect(a, SIGNAL(activated()), this, SLOT(editRedo()));
    a->addTo(menu);
    menu->insertSeparator();
    a = new QAction(tr("Copy"), tr("&Copy"), CTRL + Key_C, this, "actionEditCopy");
    connect(a, SIGNAL(activated()), this, SLOT(editCopy()));
    a->addTo(menu);
    a = new QAction(tr("Cut"), tr("Cu&t"), CTRL + Key_X, this, "actionEditCut");
    connect(a, SIGNAL(activated()), this, SLOT(editCut()));
    a->addTo(menu);
    a = new QAction(tr("Paste"), tr("&Paste"), CTRL + Key_V, this, "actionEditPaste");
    connect(a, SIGNAL(activated()), this, SLOT(editPaste()));
    a->addTo(menu);
}

void TextEdit::setupTextActions()
{}

void TextEdit::load(const QString &f)
{
    if (!QFile::exists(f))
        return;
    QTextEdit *edit = new QTextEdit(tabWidget);
    edit->setTextFormat(PlainText);
    doConnections(edit);
    tabWidget->addTab(edit, QFileInfo(f).fileName());
    QFile file(f);
    if (!file.open(IO_ReadOnly))
        return;
    QTextStream ts(&file);
    edit->setText(ts.read());
    tabWidget->showPage(edit);
    edit->viewport()->setFocus();
    filenames.replace(edit, f);

    files_open.push_back(f);
}

QTextEdit *TextEdit::currentEditor() const
{
    if (tabWidget->currentPage() &&
        tabWidget->currentPage()->inherits("QTextEdit"))
        return (QTextEdit*)tabWidget->currentPage();
    return 0;
}

void TextEdit::doConnections(QTextEdit *e)
{}

void TextEdit::fileNew()
{
    QTextEdit *edit = new QTextEdit(tabWidget);
    edit->setTextFormat(PlainText);
    doConnections(edit);
    tabWidget->addTab(edit, tr("noname"));
    tabWidget->showPage(edit);
    edit->viewport()->setFocus();
}

void TextEdit::fileOpen()
{
    QString fn = QFileDialog::getOpenFileName(QString::null, tr("All Files (*)"), this);
    if (!fn.isEmpty())
        load(fn);
}

void TextEdit::fileSave()
{
    if (!currentEditor())
        return;
    QString fn;
    if (filenames.find(currentEditor()) == filenames.end())
    {
        fileSaveAs();
    }
    else
    {
        QFile file(*filenames.find(currentEditor()));
        if (!file.open(IO_WriteOnly))
            return;
        QTextStream ts(&file);
        ts << currentEditor()->text();
    }
}

void TextEdit::fileSaveAs()
{
    if (!currentEditor())
        return;
    QString fn = QFileDialog::getSaveFileName(QString::null, tr("HTML-Files (*.htm *.html);;All Files (*)"), this);
    if (!fn.isEmpty())
    {
        filenames.replace(currentEditor(), fn);
        fileSave();
        tabWidget->setTabLabel(currentEditor(), QFileInfo(fn).fileName());
    }
}

void TextEdit::filePrint()
{
    if (!currentEditor())
        return;
#ifndef QT_NO_PRINTER
    QPrinter printer;
    printer.setFullPage(TRUE);
    if (printer.setup(this))
    {
        QPainter p(&printer);
        // Check that there is a valid device to print to.
        if (!p.device()) return;
        QPaintDeviceMetrics metrics(p.device());
        int dpix = metrics.logicalDpiX();
        int dpiy = metrics.logicalDpiY();
        const int margin = 72; // pt
        QRect body(margin * dpix / 72, margin * dpiy / 72,
                   metrics.width() - margin * dpix / 72 * 2,
                   metrics.height() - margin * dpiy / 72 * 2);
        QFont font("Courier", 10);
        QSimpleRichText richText(currentEditor()->text(), font, currentEditor()->context(), currentEditor()->styleSheet(),
                                 currentEditor()->mimeSourceFactory(), body.height());
        richText.setWidth(&p, body.width());
        QRect view(body);
        int page = 1;
        do
        {
            richText.draw(&p, body.left(), body.top(), view, colorGroup());
            view.moveBy(0, body.height());
            p.translate(0 , -body.height());
            p.setFont(font);
            p.drawText(view.right() - p.fontMetrics().width(QString::number(page)),
                       view.bottom() + p.fontMetrics().ascent() + 5, QString::number(page));
            if (view.top()  >= richText.height())
                break;
            printer.newPage();
            page++;
        }
        while (TRUE);
    }
#endif
}

void TextEdit::fileClose()
{
    delete currentEditor();
    if (currentEditor())
        currentEditor()->viewport()->setFocus();
}

void TextEdit::fileExit()
{
    if (remove)
        for (int i = 0; i < files_open.size(); i++)
            system(((QString)"rm " + files_open[i]).ascii());
    qApp->quit();
}

void TextEdit::editUndo()
{
    if (!currentEditor())
        return;
    currentEditor()->undo();
}

void TextEdit::editRedo()
{
    if (!currentEditor())
        return;
    currentEditor()->redo();
}

void TextEdit::editCut()
{
    if (!currentEditor())
        return;
    currentEditor()->cut();
}

void TextEdit::editCopy()
{
    if (!currentEditor())
        return;
    currentEditor()->copy();
}

void TextEdit::editPaste()
{
    if (!currentEditor())
        return;
    currentEditor()->paste();
}

void TextEdit::editorChanged(QWidget *)
{
    if (!currentEditor())
        return;
}

void TextEdit::addAccel(int accel, const QString &function)
{
    QAccel *a = new QAccel(this);
    a->insertItem(accel);
    connect(a, SIGNAL(activated(int)), this, SLOT(execAccel()));
    accels.insert(a, function);
}

void TextEdit::execAccel()
{
    QString func = *accels.find((QAccel*)sender());
}
