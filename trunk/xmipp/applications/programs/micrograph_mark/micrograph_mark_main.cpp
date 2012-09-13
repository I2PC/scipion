/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
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

#include <data/args.h>

#include <qapplication.h>

#include "main_widget_mark.h"
#include "widget_micrograph.h"

class ProgMicrographMark: public XmippProgram
{
private:

protected:
    void defineParams()
    {
        addUsageLine("Mark particles in a micrograph");
        addUsageLine("+This utility allows you to mark and cut micrographs ");
        addUsageLine("+as well as to generate a list of coordinates with the positions to cut.");
        addUsageLine("+The program admits single or tilted couples of micrographs (markpair mode).");
        addUsageLine("+ ");
        addUsageLine("+*Manual particle selection:*");
        addUsageLine("+1. Click with the mouse left button on the overview in order to select ");
        addUsageLine("+the part of the micrograph you desire");
        addUsageLine("+1. Mark particles with the mouse left button on the micrograph window");
        addUsageLine("+If you are in tilt pair mode, select particles in the untilted image; the corresponding point will be drawn in the tilted micrograph. If you don't like the position of the automatically generated particle you can move by clicking with the left mouse button");
        addUsageLine("+1. To remove a particle, click again with the left button on the particle.");
        addUsageLine("+1. To move a particle, click with the right button on the particle and drag it with the left button.");
        addUsageLine("+ ");
        addUsageLine("+*Automatic particle selection:*");
        addUsageLine("+The algorithm is designed to learn the particles from the user, as well as from its own errors.");
        addUsageLine("+The algorithm is fully described in [[http://www.ncbi.nlm.nih.gov/pubmed/19555764][this paper]].");
        addUsageLine("+1. Pick one particle in the first micrograph");
        addUsageLine("+1. Adjust the mark size to cover the particle diameter");
        addUsageLine("+1. Pick the rest of the particles for the first micrograph.");
        addUsageLine("+1. Click on Learn, Save and Quit");
        addUsageLine("+1. Open the second micrograph, pick all the particles in this micrograph.");
        addUsageLine("+1. Click on Learn, Save and Quit");
        addUsageLine("+1. Repeat this process until about 150 particles have been picked");
        addUsageLine("+1. From this point, you may teach the algorithm about difficult particles, for doing so");
        addUsageLine("+1. Open the micrograph");
        addUsageLine("+1. Click on AutoSelect and wait for the results");
        addUsageLine("+1. Correct the algorithm by removing the wrongly selected particles and by adding the unselected particles");
        addUsageLine("+1. Click on Learn, Save and Quit");
        addUsageLine("+1. Repeat this process until the algorithm produces satisfactory results (normally this state is achieved after about 300 particles have been manually picked).");
        addUsageLine("+1. From this point, you may autoselect the rest of the micrographs without supervision.");
        addUsageLine("+ ");
        addUsageLine("+*Key shortcuts:*");
        addUsageLine("+* Ctrl-+ and Ctrl-- in the overview window to Zoom-in and Zoom-out respectively");
        addUsageLine("+* Ctrl-G to adjust the contrast");
        addUsageLine("+* Ctrl-S to save the particle coordinates. If you are marking a tilt pair, you must press Ctrl-S on both images.");
        addUsageLine("+* Ctrl-R to change the radius of the mark. NOTE: if you are working in tilted mode left and right windows have independent radius");
        addUsageLine("+ ");
        addUsageLine("+*Notes:*");
        addUsageLine("+* When loading particles in two micrographs (tilting pairs) the conversion matrix (the one which passes coordinates from one micrograph to the other) is not computed. For computing it, click on Save angles");
        addUsageLine("+* Save angles is an option for paired micrographs which computes the tilting angle of the tilted micrograph, and the angle from the Y axis to the tilt axis (clockwise angles are positive) for both images. The result is written in a file called as the untilted image with the extension .ang");
        addUsageLine("+* When cutting the micrograph into particles you can choose a different file from which to cut. This is used if you want to mark on a downsampled image and cut from the original one. You must supply the root name for the particles and the window size. This window size will be applied directly without any scaling to the micrograph to be cut. The popup window that appears after selecting the option generate images will allow you to invert contrast and compute optical density from transmitance data");
        addUsageLine("+* You can apply filters to the image being visualized, so that you can see particles better. Filters are added into a queue, so whenever you change the visualized area, they are applied again. Select the menu filter in the overview and add so many filters as you like. You can clean the filter queue, if you want to return to the original image");
        addUsageLine("+* IMPORTANT: when working in markpair mode you the tilt and untilted image are independent. For example: if you save the position of the particles in the untilted window you are not saving it in the tilted (by the way, be sure to provide two DIFFERENT file names for the tilt and untilted positions)");
        addUsageLine("+ ");
        addUsageLine("+*Conversion from other programs:*");
        addUsageLine("+* [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/BoxerToXmippMark][From Boxer (Eman)]]");
        addUsageLine("+* [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebToXmippMark][From Web (Spider)]]");
        addParamsLine("  -i <input_untilted_micrograph>      : File with untilted image");
        addParamsLine("     alias --untilted;");
        addParamsLine("  [-t <input_tilted_micrograph>]      : File with tilted image");
        addParamsLine("     alias --tilted;");
        addParamsLine("  [--auto <model_rootname>]           : For autoselection\n");
        addParamsLine("  [--autoSelect]                      : Autoselect without user interaction\n");
        addParamsLine("  [--thr <p=1>]                       : Number of threads for automatic picking\n");
        addParamsLine("  [--outputRoot+ <rootname>]          : Output rootname\n");
        addParamsLine("                                      : If not given, the micrograph name is taken\n");
        addExampleLine("Mark a single image:",false);
        addExampleLine("xmipp_micrograph_mark -i micrograph.tif");
        addExampleLine("Mark a tilt pair:",false);
        addExampleLine("xmipp_micrograph_mark -i untilted.tif --tilted tilted.tif");
        addExampleLine("Train the autoselection algorithm:",false);
        addExampleLine("xmipp_micrograph_mark -i untilted.tif --auto model");
        addExampleLine("Automatically select particles in a micrograph (without supervision):",false);
        addExampleLine("xmipp_micrograph_mark -i untilted.tif --auto model --autoSelect");
    }
    FileName fnRaw;
    FileName fnRawTilted;
    FileName fnAutomaticModel;
    FileName outputRoot;
    bool     autoSelect;
    int      numThreads;

    void readParams()
    {
        autoSelect = false;
        // Get input parameters .................................................
        fnRaw         = getParam( "-i");
        if(checkParam("--tilted"))
            fnRawTilted   = getParam("--tilted");
        else
            fnRawTilted   ="";
        if(checkParam("--auto"))
            fnAutomaticModel = getParam("--auto");
        else
           fnAutomaticModel ="";
        autoSelect = checkParam("--autoSelect");
        if (fnRawTilted!="" && autoSelect)
            REPORT_ERROR(ERR_VALUE_INCORRECT,"Automatic particle picking cannot be performed on tilt pairs");
        numThreads = getIntParam("--thr");
        if (checkParam("--outputRoot"))
            outputRoot = getParam("--outputRoot");
        else
            outputRoot = fnRaw;
    }
public:
    void run()
    {
        Micrograph m, mTilted;
        FileName fn8bits="", fn8bitsTilted="";

        m.open_micrograph(fnRaw);
        std::vector<FileName> filesToDelete;
        fn8bits=fnRaw+"_8bits.raw";
        m.write(fn8bits+"%uint8", CW_ADJUST);
        m.close_micrograph();

        m.open_micrograph(fn8bits);
        m.set_micrograph_name(fnRaw);
        filesToDelete.push_back(fn8bits+"*");

        if (fnRawTilted != "")
        {
            mTilted.open_micrograph(fnRawTilted);
            fn8bitsTilted=fnRawTilted+"_8bits.raw";
            mTilted.write(fn8bitsTilted+"%uint8");
            mTilted.close_micrograph();
            mTilted.open_micrograph(fn8bitsTilted);
            mTilted.set_micrograph_name(fnRawTilted);
            filesToDelete.push_back(fn8bitsTilted+"*");
        }

        // Configure application .............................................
        AutoParticlePickingQt *autoPicking=NULL;
        if (fnAutomaticModel!="")
            autoPicking=new AutoParticlePickingQt(&m);

        QApplication *app=NULL;
        QtMainWidgetMark *mainWidget=NULL;
        if (!autoSelect)
        {
            app=new QApplication(argc, argv);
            if (fnRawTilted == "")
            {
                mainWidget = new QtMainWidgetMark(&m);
                mainWidget->untilted_widget()->
                setAutoParticlePicking(autoPicking);
            }
            else
                mainWidget = new QtMainWidgetMark(&m, &mTilted);
            for (int i=0; i<filesToDelete.size(); i++)
                mainWidget->__filesToDelete.push_back(filesToDelete[i]);
        }

        // Check if a model has been provided ................................
        if (autoPicking!=NULL)
            autoPicking->setNumThreads(numThreads);
        if (fnAutomaticModel!="")
            autoPicking->loadModels(fnAutomaticModel);

        // Run application ...................................................
        if (!autoSelect)
        {
            mainWidget->setOutputRoot(outputRoot);
            if (autoPicking!=NULL)
                autoPicking->setOutputRoot(outputRoot);
            app->setMainWidget(mainWidget);
            mainWidget->openAllWindows();
            app->exec();
        }
        else
        {
            autoPicking->setOutputRoot(outputRoot);
            autoPicking->automaticallySelectParticles();
            autoPicking->saveAutoParticles();
        }
    }
};

int main(int argc, char **argv)
{
    ProgMicrographMark program;
    program.read(argc, argv);
    return program.tryRun();
}
