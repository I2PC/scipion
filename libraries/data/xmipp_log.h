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

#ifndef LOG_H_
#define LOG_H_

#include "xmipp_strings.h"

/**
 * This class will serve to log messages
 */
class XmippLog
{
private:
  String levelSpaces; //This variable will be used to set level spaces
  /* File handler */
  FILE * fileHandler;

public:
  /* Constructor */
  XmippLog(String filename);
  /* Destructor */
  ~XmippLog();

  void logMessage(const String &msg);
  /* Increase by 3 the level spaces */
  void enterLevel(const String &msg="");
  /* Decrease level spaces */
  void leaveLevel(const String &msg="");

};//class XmippLog

/**
 * Class to log entering to functions
 * and other scopes automatically
 */
class XmippLogBlock
{
private:
  XmippLog * log;
  String block;

public:
  /* Constructor */
  XmippLogBlock(XmippLog * log, const String &block);
  /* Destructor */
  ~XmippLogBlock();
};//class XmippLogBlock

/** Global log pointer and macros to use it */
extern XmippLog * __xmippLog;

#define LOG_ENABLED 1


#ifdef LOG_ENABLED
#define LOG_FN(root) formatString("%s_nodo%02d_debug.log", root.c_str(), rank)
#define CREATE_LOG(filename) __xmippLog = new XmippLog(filename)
#define LOG(msg) __xmippLog->logMessage(msg)
#define CLOSE_LOG() delete __xmippLog
#define LOG_LEVEL(var) XmippLogBlock __logBlock##var(__xmippLog, #var)
#else
#define CREATE_LOG(filename) ;
#define LOG(msg) ;
#define CLOSE_LOG() ;
#define LOG_LEVEL(msg) ;
#endif


#endif /* LOG_H_ */
