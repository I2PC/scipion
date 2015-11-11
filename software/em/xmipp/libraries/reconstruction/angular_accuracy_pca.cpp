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

#include "angular_accuracy_pca.h"

ProgAngularAccuracyPCA::ProgAngularAccuracyPCA()
{
	rank=0;
	Nprocessors=1;
}

void ProgAngularAccuracyPCA::readParams()
{
	fnParticles = getParam("--i1");
	fnNeighbours = getParam("--i2");
}

void ProgAngularAccuracyPCA::defineParams()
{
    addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
    addParamsLine("  [--i1 <md_file=\"\">]      	: Metadata file with input projections");
    addParamsLine("  [--i2 <md_file=\"\">]    		: Metadata file with neighbour projections");
    addParamsLine("  [--o <outputDir=\".\">] 		: Output directory");
}

void ProgAngularAccuracyPCA::run()
{
	MetaData md;
    StringVector blocks;
    std::cout << "Blocks in " << fnParticles << ": " << std::endl;
    getBlocksInMetaDataFile(fnNeighbours, blocks);
    for (size_t i = 0; i < blocks.size(); ++i)
    	md=read(blocks[i]+fnParticles);

    std::cout << md << std::endl;

    	//std::cout << blocks[i] << std::endl;

}
