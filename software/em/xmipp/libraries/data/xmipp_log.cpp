/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <stdio.h>
#include <iostream>
#include "xmipp_log.h"
#include "xmipp_funcs.h"

/** Global log variable */
XmippLog * __xmippLog;

XmippLog::XmippLog(String filename)
{
  fileHandler = fopen(filename.c_str(), "w+");
}

XmippLog::~XmippLog()
{
  fclose(fileHandler);
}

void XmippLog::logMessage(const String &msg)
{
  //printf("%s%s\t%s\n", getCurrentTimeString(), levelSpaces.c_str(), msg.c_str());
  fprintf(fileHandler, "%s%s   %s\n", getCurrentTimeString(), levelSpaces.c_str(), msg.c_str());
  fflush(fileHandler);
}

void XmippLog::enterLevel(const String &msg)
{
  logMessage(formatString(">> %s", msg.c_str()));
  levelSpaces += "   ";
}

void XmippLog::leaveLevel(const String &msg)
{
  levelSpaces.erase(0, 3);
  if (msg.size())
    logMessage(formatString("<< %s", msg.c_str()));
}

XmippLogBlock::XmippLogBlock(XmippLog * log, const String &block)
{
  this->log = log;
  this->block = block;
  log->enterLevel(block);
}

XmippLogBlock::~XmippLogBlock()
{
  log->leaveLevel();
}

