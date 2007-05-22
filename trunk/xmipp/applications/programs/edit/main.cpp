/****************************************************************************
** $Id$
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#include <qapplication.h>

#include "textedit.h"

#include <data/args.h>

void usage()
{
    cerr << "edit\n"
    << "   -i <filenames>    : Files to open\n"
    << "  [-remove]          : Remove the files after finishing\n"
    ;
}

int main(int argc, char ** argv)
{
    int ifirst;
    bool remove;

    try
    {
        ifirst = position_param(argc, argv, "-i");
        if (ifirst == -1 && argc != 1) REPORT_ERROR(1, "Edit: Cannot find -i");
        remove = check_param(argc, argv, "-remove");
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        usage();
        exit(1);
    }


    QApplication a(argc, argv);
    TextEdit * mw = new TextEdit(0, "TextEdit");
    if (remove) mw->setRemoveFlag();
    mw->setCaption("Xmipp Edit");

    for (int i = ifirst + 1; i < argc; i++)
        if (exists(argv[i])) mw->load(argv[i]);

    mw->resize(500, 300);
    mw->show();
    a.connect(&a, SIGNAL(lastWindowClosed()), mw, SLOT(fileExit()));
    return a.exec();
}
