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

import java.awt.Image;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.viewer.ImageDimension;

public class RotSpectraGalleryTableModel extends MetadataGalleryTableModel {

	final static int DEFAULT_DIM = 128;
	final static String blockKerDenSOM = "KerDenSOM_Layout"
			+ Filename.SEPARATOR;
	final static String blockVectorContent = "vectorContent"
			+ Filename.SEPARATOR;
	int vectorsSize;
	ArrayList<RotSpectraItem> vectorItems;

	String fnClasses, fnVectors, fnVectorsData;
	boolean showLabels;

	public RotSpectraGalleryTableModel(GalleryData data) throws Exception {
		super(data);
		DEBUG.printMessage("creating rotspectra.....");
	}

	// Load initial dimensions
	@Override
	protected ImageDimension loadDimension() throws Exception {
		fnClasses = data.getFileName();
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
		cols = (int)mdClasses.getValueLong(MDLabel.MDL_XSIZE, id);
		rows = (int)mdClasses.getValueLong(MDLabel.MDL_YSIZE, id);
		int nvectors = cols * rows;

		id = mdVectors.firstObject();
		vectorsSize = (int)mdVectors.getValueLong(
				MDLabel.MDL_CLASSIFICATION_DATA_SIZE, id);

		double vectors[][] = loadVectors(filenameData, nvectors, vectorsSize);

		vectorItems = new ArrayList<RotSpectraItem>(nvectors);
		
		for (int i = 0; i < vectors.length; i++) 
			vectorItems.add(new RotSpectraItem(i, vectors[i]));
		
		mdClasses.destroy();
		mdVectors.destroy();
	}

	static double[][] loadVectors(String filename, int nvectors, int size)
			throws Exception {
		double vectors[][] = new double[nvectors][size];

		// Load vectors.
		DataInputStream in = new DataInputStream(new FileInputStream(filename));
		byte[] buffer = new byte[4]; // This is for reading little endian 
		for (int i = 0; i < nvectors; i++) {
			for (int j = 0; j < size; j++) {
				in.read(buffer);
				int a = (buffer[0] & 0xFF) | (buffer[1] & 0xFF) << 8 | (buffer[2] & 0xFF) << 16 | (buffer[3] & 0xFF) << 24;
				vectors[i][j] = Float.intBitsToFloat(a);
			}
		}
		in.close();

		return vectors;
	}

	public Object getValueAt(int row, int col) {
		int index = getIndex(row, col);
		if (index < n) {
			RotSpectraItem item = vectorItems.get(index);
			//item.cellDim.setSize(thumb_width, thumb_height);
			return item;
		}
		return null;
	}
	
	// Return True if columns were changed
	public boolean adjustColumn(int width) {
		return false;
	}
	
	@Override
	public boolean handleDoubleClick(int row, int col){
		int index = getIndex(row, col);
		if (index < n) {
			RotSpectraItem item = vectorItems.get(index);
			//item.cellDim.setSize(thumb_width, thumb_height);
			JFrame chart = item.getChart();
			chart.setVisible(true);
			return true;
		}
		return false;
	}//function handleDoubleClick
	
	public class RotSpectraItem extends ImageItem{

	    JFreeChart chart;
	    double vector[];

	    public RotSpectraItem(int index, double vector[]) {
	        super(index);
	        this.vector = vector;
	        chart = createChart();
	    }

	    /** Create chart based on vector data */
	    private JFreeChart createChart() {
	        XYSeries series = new XYSeries(getLabel());

	        for (int i = 0; i < vector.length; i++) 
	            series.add((double) i, vector[i]);

	        XYSeriesCollection dataset = new XYSeriesCollection(series);

	        JFreeChart chart = ChartFactory.createXYLineChart(
	                "", "", "",
	                dataset, PlotOrientation.VERTICAL,
	                true, true, false);

	        chart.removeLegend();

	        return chart;
	    }//function createChart

	    @Override
	    public Image getImage() {
	        return chart.createBufferedImage(cellDim.width, cellDim.height - font_height);
	    }

	    public JFrame getChart() {
	        ChartPanel panel = new ChartPanel(chart);
	        JFrame frame = new JFrame(getLabel());
	        frame.getContentPane().add(panel);
	        frame.pack();
	        return frame;
	    }
	}//class RotSpectraItem
}// class RotSpectraGallery
