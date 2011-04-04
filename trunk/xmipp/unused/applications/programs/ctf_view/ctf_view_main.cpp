/***************************************************************************
 *
 * Authors:     Javier Rodr�guez Fern�ndez (javrodri@gmail.com)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Universidad San Pablo CEU (Monteprincipe, Madrid)
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

#include <qapplication.h>

#include <data/program.h>
#include "module.h"

class ProgCtfView: public XmippProgram
{
public:
  FileName fn_ctf;

  void defineParams()
  {
    addUsageLine("Display the CTF parameters");
    addUsageLine("+This program allows you visualize how a CTF fitting really fits the experimental Power Spectrum Density");
    addUsageLine("+long a line whose angle can be selected in the graphical interface. The program can draw the theoretical CTF,");
    addUsageLine("+its damping envelope, and the background noise independently. Measurements can also be done in decibels (10*log10)");
    addUsageLine("+to help visualization. The user can control the plotter number of ticks and visualization region. Press the reset");
    addUsageLine("+button (on the plotter itself, not the menu) to reset the plot settings (sometimes this is needed after resizing or zooming).");
    addParamsLine("   [-i <ctf_file=\"\">] : The file is assumed to be of kind .ctfparam");
  }

  void readParams()
  {
    fn_ctf = getParam("-i");
  }

  void createGUI()
  {
    QApplication app(argc, argv);
    CTFViewer * Viewer = new CTFViewer(0, 0, fn_ctf);
    app.setMainWidget(Viewer);
    app.exec();
  }

  void run()
  {
    createGUI();
  }
};//end of class ProgCtfView

int main(int argc, char *argv[])
{
    ProgCtfView program;
    program.read(argc, argv);
    return program.tryRun();
}
