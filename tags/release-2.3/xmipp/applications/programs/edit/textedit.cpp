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


#include <qaction.h>
#include <qmenubar.h>
#include <qtabwidget.h>
#include <qapplication.h>
#include <qfontdatabase.h>
#include <qcombobox.h>
#include <qlineedit.h>
#include <qfileinfo.h>
#include <qfile.h>
#include <qprinter.h>
#include <qcolordialog.h>
#include <qpainter.h>
#include <qinputdialog.h>

#ifdef QT3_SUPPORT
#include <q3textedit.h>
#include <q3popupmenu.h>
#include <q3toolbar.h>
#include <q3filedialog.h>
#include <q3paintdevicemetrics.h>
#include <q3simplerichtext.h>
//Added by qt3to4:
#include <QPixmap>
#else
#include <qtextedit.h>
#include <qpopupmenu.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qpaintdevicemetrics.h>
#include <qsimplerichtext.h>
#endif

TextEdit::TextEdit(QWidget *parent, const char *name)
#ifdef QT3_SUPPORT
        : Q3MainWindow(parent, name)
#else
        : QMainWindow(parent, name)
#endif
{
    setupFileActions();
    setupEditActions();
    setupTextActions();

#ifdef QT3_SUPPORT
    Q3PopupMenu *menu = new Q3PopupMenu(this);
#else
    QPopupMenu* menu = new QPopupMenu(this);
#endif

    tabWidget = new QTabWidget(this);
    connect(tabWidget, SIGNAL(currentChanged(QWidget *)),
            this, SLOT(editorChanged(QWidget *)));
    setCentralWidget(tabWidget);

    remove = false;
}

void TextEdit::setupFileActions()
{
#ifdef QT3_SUPPORT
    Q3PopupMenu *menu = new Q3PopupMenu(this);
#else
    QPopupMenu* menu = new QPopupMenu(this);
#endif

    menuBar()->insertItem(tr("&File"), menu);

    QAction *a;
    a = new QAction(tr("&New..."), Qt::CTRL + Qt::Key_N, this, "actionFileNew");
    connect(a, SIGNAL(activated()), this, SLOT(fileNew()));
    a->addTo(menu);
    a = new QAction(tr("&Open..."), Qt::CTRL + Qt::Key_O, this, "actionFileOpen");
    connect(a, SIGNAL(activated()), this, SLOT(fileOpen()));
    a->addTo(menu);
    menu->insertSeparator();
    a = new QAction(tr("&Save..."), Qt::CTRL + Qt::Key_S, this, "actionFileSave");
    connect(a, SIGNAL(activated()), this, SLOT(fileSave()));
    a->addTo(menu);
    a = new QAction(QPixmap(), tr("Save &As..."), 0, this, "actionFileSaveAs");
    connect(a, SIGNAL(activated()), this, SLOT(fileSaveAs()));
    a->addTo(menu);
    menu->insertSeparator();
    a = new QAction(tr("&Print..."), Qt::CTRL + Qt::Key_P, this, "actionFilePrint");
    connect(a, SIGNAL(activated()), this, SLOT(filePrint()));
    a->addTo(menu);
    a = new QAction(QPixmap(), tr("&Close"), Qt::CTRL + Qt::Key_W, this, "actionFileClose");
    connect(a, SIGNAL(activated()), this, SLOT(fileClose()));
    a->addTo(menu);
    a = new QAction(QPixmap(), tr("E&xit"), Qt::CTRL + Qt::Key_Q, this, "actionFileExit");
    connect(a, SIGNAL(activated()), this, SLOT(fileExit()));
    a->addTo(menu);
}

void TextEdit::setupEditActions()
{
#ifdef QT3_SUPPORT
    Q3PopupMenu *menu = new Q3PopupMenu(this);
#else
    QPopupMenu* menu = new QPopupMenu(this);
#endif

    menuBar()->insertItem(tr("&Edit"), menu);

    QAction *a;
    a = new QAction(tr("&Undo"), Qt::CTRL + Qt::Key_Z, this, "actionEditUndo");
    connect(a, SIGNAL(activated()), this, SLOT(editUndo()));
    a->addTo(menu);
    a = new QAction(tr("&Redo"), Qt::CTRL + Qt::Key_Y, this, "actionEditRedo");
    connect(a, SIGNAL(activated()), this, SLOT(editRedo()));
    a->addTo(menu);
    menu->insertSeparator();
    a = new QAction(tr("&Copy"), Qt::CTRL + Qt::Key_C, this, "actionEditCopy");
    connect(a, SIGNAL(activated()), this, SLOT(editCopy()));
    a->addTo(menu);
    a = new QAction(tr("Cu&t"), Qt::CTRL + Qt::Key_X, this, "actionEditCut");
    connect(a, SIGNAL(activated()), this, SLOT(editCut()));
    a->addTo(menu);
    a = new QAction(tr("&Paste"), Qt::CTRL + Qt::Key_V, this, "actionEditPaste");
    connect(a, SIGNAL(activated()), this, SLOT(editPaste()));
    a->addTo(menu);
}

void TextEdit::setupTextActions()
{}

void TextEdit::load(const QString &f)
{
    if (!QFile::exists(f))
        return;
#ifdef QT3_SUPPORT
    Q3TextEdit *edit = new Q3TextEdit(tabWidget);
#else
    QTextEdit* edit = new QTextEdit(tabWidget);
#endif
    edit->setTextFormat(Qt::PlainText);
    doConnections(edit);
    tabWidget->addTab(edit, QFileInfo(f).fileName());
    QFile file(f);
#ifdef QT3_SUPPORT
    if (!file.open(QIODevice::ReadOnly))
#else
    if (!file.open(IO_ReadOnly))
#endif
        return;
    QTextStream ts(&file);
    edit->setText(ts.read());
    tabWidget->showPage(edit);
    edit->viewport()->setFocus();
    filenames.replace(edit, f);

    files_open.push_back(f);
}

#ifdef QT3_SUPPORT
Q3TextEdit *TextEdit::currentEditor() const
#else
QTextEdit* TextEdit::currentEditor() const
#endif
{
    if (tabWidget->currentPage() &&
        tabWidget->currentPage()->inherits("QTextEdit"))
#ifdef QT3_SUPPORT
        return (Q3TextEdit*)tabWidget->currentPage();
#else
        return (QTextEdit*) tabWidget->currentPage();
#endif
    return 0;
}

#ifdef QT3_SUPPORT
void TextEdit::doConnections(Q3TextEdit *e)
#else
void TextEdit::doConnections(QTextEdit* e)
#endif
{}

void TextEdit::fileNew()
{
#ifdef QT3_SUPPORT
    Q3TextEdit *edit = new Q3TextEdit(tabWidget);
#else
    QTextEdit* edit = new QTextEdit(tabWidget);
#endif
    edit->setTextFormat(Qt::PlainText);
    doConnections(edit);
    tabWidget->addTab(edit, tr("noname"));
    tabWidget->showPage(edit);
    edit->viewport()->setFocus();
}

void TextEdit::fileOpen()
{
#ifdef QT3_SUPPORT
    QString fn = Q3FileDialog::getOpenFileName(QString::null, tr("All Files (*)"), this);
#else
    QString fn = QFileDialog::getOpenFileName(QString::null, tr("All Files (*)"), this);
#endif
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
#ifdef QT3_SUPPORT
        if (!file.open(QIODevice::WriteOnly))
#else
        if (!file.open(IO_WriteOnly))
#endif
            return;
        QTextStream ts(&file);
        ts << currentEditor()->text();
    }
}

void TextEdit::fileSaveAs()
{
    if (!currentEditor())
        return;
#ifdef QT3_SUPPORT
    QString fn = Q3FileDialog::getSaveFileName(QString::null, tr("HTML-Files (*.htm *.html);;All Files (*)"), this);
#else
    QString fn = QFileDialog::getSaveFileName(QString::null, tr("HTML-Files (*.htm *.html);;All Files (*)"), this);
#endif
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
#ifdef QT3_SUPPORT
        Q3PaintDeviceMetrics metrics(p.device());
#else
        QPaintDeviceMetrics metrics(p.device());
#endif
        int dpix = metrics.logicalDpiX();
        int dpiy = metrics.logicalDpiY();
        const int margin = 72; // pt
        QRect body(margin * dpix / 72, margin * dpiy / 72,
                   metrics.width() - margin * dpix / 72 * 2,
                   metrics.height() - margin * dpiy / 72 * 2);
        QFont font("Courier", 10);
#ifdef QT3_SUPPORT
        Q3SimpleRichText richText(currentEditor()->text(), font, currentEditor()->context(), currentEditor()->styleSheet(),
                                 currentEditor()->mimeSourceFactory(), body.height());
#else
        QSimpleRichText richText(currentEditor()->text(), font, currentEditor()->context(), currentEditor()->styleSheet(),
                                 currentEditor()->mimeSourceFactory(), body.height());
#endif
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
#include <cstdlib>
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
#ifdef QT3_SUPPORT
    Q3Accel *a = new Q3Accel(this);
#else
    QAccel* a = new QAccel(this);
#endif
    a->insertItem(accel);
    connect(a, SIGNAL(activated(int)), this, SLOT(execAccel()));
    accels.insert(a, function);
}

void TextEdit::execAccel()
{
#ifdef QT3_SUPPORT
    QString func = *accels.find((Q3Accel*)sender());
#else
    QString func = *accels.find((QAccel*) sender());
#endif
}
