/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
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


#include "volume_validate_pca.h"
#include <numeric>
#include "data/geometry.h"

// Define params
void ProgVolumeValidationPCA::defineParams()
{
    //usage
    addUsageLine("Validate an obtained volume from a set of class averages");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input classes");
    addParamsLine("  [ --vol <file=\"\">]         : Input volume");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--numVols <N=5>]            : Number of intermediate volumes to generate");
    addParamsLine("  [--numClasses <N=8>]         : Number of classes to generate the intermediate volumes");
}


// Read arguments ==========================================================
void ProgVolumeValidationPCA::readParams()
{
    fnClasses = getParam("-i");
    fnSym = getParam("--sym");
    NVols = getIntParam("--numVols");
    NClasses = getIntParam("--numClasses");
}

// Show ====================================================================
void ProgVolumeValidationPCA::show()
{
    if (verbose > 0)
    {
        std::cout << "Input classes metadata           : "  << fnClasses  << std::endl;
        std::cout << "Number of intermediate volumes to generate : "  << NVols      << std::endl;
        std::cout << "Number of classes to be used               : "  << NClasses   << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
    }
}

void ProgVolumeValidationPCA::produceSideinfo()
{
    mdClasses.read(fnClasses);
    mdClasses.removeDisabled();
    getImageSize(mdClasses,xdim, ydim, zdim, ndim);
}

void ProgVolumeValidationPCA::modifyAngles()
{
    //String args=formatString("-i %s -o %s --operate random_subset %d --mode overwrite",fnClasses.c_str(),fnAngles.c_str(),NClasses);
    String args=formatString("-i %s -o %s --operate bootstrap --mode overwrite",fnClasses.c_str(),fnAngles.c_str());
    String cmd=(String)"xmipp_metadata_utilities "+args;
    system(cmd.c_str());

    //Metadata with the well sampled projection and random projections assigned
    MetaData md;
    MDRow row;
    double rot, tilt, psi;
    size_t id;

    md.read(fnAngles);
    randomize_random_generator();

    int anglePertur = 30;
    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        id = __iter.objId;

        md.getRow(row, id);

        row.getValue(MDL_ANGLE_ROT, rot);
        row.getValue(MDL_ANGLE_TILT,tilt);
        row.getValue(MDL_ANGLE_PSI,psi);

        row.setValue(MDL_ANGLE_ROT, (rnd_unif(-anglePertur,anglePertur))/2+rot);
        row.setValue(MDL_ANGLE_TILT,(rnd_unif(-anglePertur,anglePertur))/2+tilt);
        row.setValue(MDL_ANGLE_PSI, (rnd_unif(-anglePertur,anglePertur))/2+psi);

        md.setRow(row, id);
    }

    md.write(fnAngles);
}

void ProgVolumeValidationPCA::reconstruct()
{
    //number of threads in the reconstruction
    int Nthr = 4;
    float AngularSampling = 30;
    int NProj = 10;

    String args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVol.c_str(),fnSym.c_str(),Nthr);
    String cmd=(String)"xmipp_reconstruct_fourier "+args;
    system(cmd.c_str());

    args=formatString("-i %s -o %s --mask circular %d -v 0",fnVol.c_str(),fnVol.c_str(),-xdim/2);
    cmd=(String)"xmipp_transform_mask "+args;
    system(cmd.c_str());


    if (doProjMatch)
    {

        for (int ind=0; ind<NProj; ind++)
        {

            args=formatString("-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s",fnVol.c_str(),fnGallery.c_str(),AngularSampling,fnSym.c_str(),fnAngles.c_str());
            cmd=(String)"xmipp_angular_project_library "+args;
            system(cmd.c_str());

            args=formatString("-i %s -o %s --ref %s --Ri 0 --Ro %f --max_shift %f --append",fnAngles.c_str(),fnAngles.c_str(),fnGallery.c_str(),xdim/2,xdim/20);
            cmd=(String)"xmipp_angular_projection_matching "+args;
            system(cmd.c_str());

            args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVol.c_str(),fnSym.c_str(),Nthr);
            cmd=(String)"xmipp_reconstruct_fourier "+args;
            system(cmd.c_str());

            args=formatString("-i %s -o %s --mask circular %d -v 0",fnVol.c_str(),fnVol.c_str(),-xdim/2);
            cmd=(String)"xmipp_transform_mask "+args;
            system(cmd.c_str());

        }
    }

    args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVol.c_str(),fnSym.c_str(),Nthr);
    cmd=(String)"xmipp_reconstruct_fourier "+args;
    system(cmd.c_str());

    args=formatString("-i %s -o %s --mask circular %d -v 0",fnVol.c_str(),fnVol.c_str(),-xdim/2);
    cmd=(String)"xmipp_transform_mask "+args;
    system(cmd.c_str());

    args=formatString("-i %s -o %s --select below 0 --substitute value 0 -v 0",fnVol.c_str(),fnVol.c_str());
    cmd=(String)"xmipp_transform_threshold "+args;
    system(cmd.c_str());
}

void ProgVolumeValidationPCA::evaluate()
{
    MetaData md, md2;
    MDRow row;
    FileName imag;
    double rot, tilt, psi, rot2, tilt2, psi2;
    std::vector<double> error;
    size_t id, id2;

    md.read(fnClasses);
    std::stringstream ss;

    MDMultiQuery query;
    MetaData auxMetadata;
    double distance = 0;

    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        id = __iter.objId;

        md.getRow(row, id);
        row.getValue(MDL_IMAGE, imag);
        row.getValue(MDL_ANGLE_ROT, rot);
        row.getValue(MDL_ANGLE_TILT,tilt);
        row.getValue(MDL_ANGLE_PSI,psi);

        for (int index=0; index<NVols;index++)
        {
            ss << index;
            fnAngles=fnClasses.removeAllExtensions()+ss.str();
            fnAngles+=".xmd";
            ss.str(std::string());
            md2.read(fnAngles);
            MDValueEQ eq(MDL_IMAGE, imag);
            query.addAndQuery(eq);
            auxMetadata.importObjects(md2, eq);

            for(MDIterator __iter2(auxMetadata); __iter2.hasNext(); __iter2.moveNext())
            {
                id2 = __iter2.objId;

                auxMetadata.getRow(row, id2);
                row.getValue(MDL_ANGLE_ROT, rot2);
                row.getValue(MDL_ANGLE_TILT,tilt2);
                row.getValue(MDL_ANGLE_PSI,psi2);
                distance = Euler_distanceBetweenAngleSets(rot,tilt, psi,rot2, tilt2, psi2, true);
                error.push_back(distance);
            }
        }
    }

    double sum = std::accumulate(error.begin(), error.end(), 0.0);
    double mean = sum / error.size();
    double sq_sum = std::inner_product(error.begin(), error.end(), error.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / error.size() - mean * mean);

    std::cout << std::endl;
    std::cout << "stdev of errors : " << stdev << std::endl;
    std::cout << "mean of errors : " << mean << std::endl;

}

void ProgVolumeValidationPCA::run()
{
    show();
    produceSideinfo();
    int index=0;

    std::stringstream ss;

    //Here we read the "mean" volume
    //fnAngles=fnClasses;
    //fnVol=fnClasses.removeAllExtensions();
    //fnVol+=".vol";
    //reconstruct();

    doProjMatch = true;

    for( int index = 0; index< NVols; index++)
    {

        ss << index;
        fnAngles=fnClasses.removeAllExtensions()+ss.str();
        fnAngles+=".xmd";
        ss.str(std::string());

        modifyAngles();

        ss << index;
        fnVol=fnClasses.removeAllExtensions()+ss.str();
        fnVol+=".vol";
        ss.str(std::string());

        ss << index;
        fnGallery=fnClasses.removeAllExtensions()+ss.str();
        fnGallery+=".stk";
        ss.str(std::string());

        //Reconstruct and ProjMatching of current
        reconstruct();

        //Some clean up
        FileName tempSampling = fnVol.removeAllExtensions();
        tempSampling+="_sampling.xmd";
        FileName tempDoc = fnAngles.removeAllExtensions();
        tempDoc+=".doc";

        //String args=formatString(" %s %s %s %s",fnVol.c_str(),fnGallery.c_str(),tempSampling.c_str(),tempDoc.c_str());
        //String cmd=(String)"rm "+args;
        //system(cmd.c_str());
    }

    evaluate();
}
