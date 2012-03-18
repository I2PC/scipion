/***************************************************************************
 * Authors:     Juanjo Vega
 * 				J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

package xmipp.viewer.models;

import java.awt.Dimension;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.ArrayList;

import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.viewer.ImageDimension;
import xmipp.viewer.RotSpectraItem;

public class RotSpectraGallery extends MetadataGallery {

	final static int DEFAULT_DIM = 128;
	final static String blockKerDenSOM = "KerDenSOM_Layout"
			+ Filename.SEPARATOR;
	final static String blockVectorContent = "vectorContent"
			+ Filename.SEPARATOR;
	int vectorsSize;
	ArrayList<RotSpectraItem> vectorItems;

	String fnClasses, fnVectors, fnVectorsData;
	boolean showLabels;

	public RotSpectraGallery(GalleryData data) throws Exception {
		super(data);
		DEBUG.printMessage("creating rotspectra.....");
	}

	// Load initial dimensions
	@Override
	protected ImageDimension loadDimension() throws Exception {
		fnClasses = data.filename;
		fnVectors = fnClasses.replace("classes", "vectors");
		fnVectorsData = fnVectors.replace(".xmd", ".vec");
		load(fnClasses, fnVectors, fnVectorsData);
		return new ImageDimension(DEFAULT_DIM, DEFAULT_DIM, data.ids.length);
	}

	void load(String filenameClasses, String filenameVectors,
			String filenameData) throws Exception {
		MetaData mdClasses = new MetaData(blockKerDenSOM + filenameClasses);
		MetaData mdVectors = new MetaData(filenameVectors);

		long id = mdClasses.firstObject();
		cols = mdClasses.getValueInt(MDLabel.MDL_XSIZE, id);
		rows = mdClasses.getValueInt(MDLabel.MDL_YSIZE, id);
		int nvectors = cols * rows;

		id = mdVectors.firstObject();
		vectorsSize = mdVectors.getValueInt(
				MDLabel.MDL_CLASSIFICATION_DATA_SIZE, id);

		double vectors[][] = loadVectors(filenameData, nvectors, vectorsSize);

		vectorItems = new ArrayList<RotSpectraItem>(nvectors);
		String key;
		RotSpectraItem item;
		
		for (int i = 0; i < vectors.length; i++) {
			key = String.format("class %d", i);
			item = new RotSpectraItem(key, key, vectors[i]);
			item.cellDim = new Dimension();
			vectorItems.add(item);
		}
	}

	static double[][] loadVectors(String filename, int nvectors, int size)
			throws Exception {
		double vectors[][] = new double[nvectors][size];

		// Load vectors.
		DataInputStream in = new DataInputStream(new FileInputStream(filename));

		for (int i = 0; i < nvectors; i++) {
			for (int j = 0; j < size; j++) {
				vectors[i][j] = in.readFloat();
			}
		}
		in.close();

		return vectors;
	}

	public Object getValueAt(int row, int col) {
		int index = getIndex(row, col);
		if (index < n) {
			RotSpectraItem item = vectorItems.get(index);
			item.showLabel = data.showLabel;
			item.cellDim.setSize(thumb_width, thumb_height);
			return item;
		}
		return null;
	}
	
	// Return True if columns were changed
	public boolean adjustColumn(int width) {
		return false;
	}
}// class RotSpectraGallery
