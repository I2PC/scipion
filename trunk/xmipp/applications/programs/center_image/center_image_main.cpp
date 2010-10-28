/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/filters.h>

class ProgCenterImage: public XmippMetadataProgram
{
public:
    int Niter;
    bool limitShift;

    void defineParams()
    {
        XmippMetadataProgram::defineParams();
       addParamsLine("  [-iter <N=10>]           : Number of iterations");
       addParamsLine("  [-limitShift]            : Limit the maximum shift allowed");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        Niter = getIntParam("-iter");
        limitShift = checkParam("-limitShift");
    }

    void show()
    {
        XmippMetadataProgram::show();
        std::cout << "Iterations = " << Niter << std::endl;
        std::cout << "LimitShift = " << limitShift << std::endl;
    }

    void processImage()
    {
      img.read(fnImg);
      img().checkDimensionWithDebug(2,__FILE__,__LINE__);
      centerImage(img(), Niter, limitShift);
    }


};///end of class ProgCenterImage

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  try
  {
      ProgCenterImage program;
      program.read(argc, argv);
      program.run();
  }
  catch (XmippError xe)
  {
      std::cerr << xe;
  }
}

