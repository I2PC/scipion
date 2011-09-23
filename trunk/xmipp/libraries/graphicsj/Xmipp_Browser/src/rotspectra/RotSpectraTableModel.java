/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rotspectra;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.LinkedList;
import javax.swing.table.AbstractTableModel;
import xmipp.Filename;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraTableModel extends AbstractTableModel {

    final static String blockKerDenSOM = "KerDenSOM_Layout" + Filename.SEPARATOR;
    final static String blockVectorContent = "vectorContent" + Filename.SEPARATOR;
    int w, h;
    int vectorsSize;
    LinkedList<RotSpectraVector> data = new LinkedList<RotSpectraVector>();

    public RotSpectraTableModel(String filenameClasses, String filenameVectors, String filenameData) throws Exception {
        load(filenameClasses, filenameVectors, filenameData);

        //int size = data.size();
        //System.out.println(size);
//        for (int i = 0; i < data.size(); i++) {
//            RotSpectraVector rsv = data.get(i);
//            Image plot = rsv.getPlotImage(128, 128);
//            ImagePlus imp = new ImagePlus(filenameData, plot);
//            imp.show();
//        }
    }

    public void load(String filenameClasses, String filenameVectors, String filenameData) throws Exception {
        MetaData mdClasses = new MetaData(blockKerDenSOM + filenameClasses);
        MetaData mdVectors = new MetaData(filenameVectors);

        long id = mdClasses.firstObject();
        w = mdClasses.getValueInt(MDLabel.MDL_XSIZE, id);
        h = mdClasses.getValueInt(MDLabel.MDL_YSIZE, id);
        int nvectors = w * h;

        String blocks[] = getBlocks(blockVectorContent + filenameVectors);

        id = mdVectors.firstObject();
        vectorsSize = mdVectors.getValueInt(MDLabel.MDL_CLASSIFICATION_DATA_SIZE, id);

        double vectors[][] = loadVectors(filenameData, nvectors, vectorsSize);

        for (int i = 0; i < vectors.length; i++) {
            double vector[] = vectors[i];
            data.add(new RotSpectraVector(
                    blocks[i], filenameData, vector));
        }
    }

    static String[] getBlocks(String filename) throws Exception {
        MetaData mdVectors = new MetaData(filename);

        long ids[] = mdVectors.findObjects();
        String blocks[] = new String[ids.length];
        for (int i = 0; i < ids.length; i++) {
            blocks[i] = mdVectors.getValueString(MDLabel.MDL_IMAGE, ids[i]);
        }

        return blocks;
    }

    static double[][] loadVectors(String filename, int nvectors, int size) throws Exception {
        double vectors[][] = new double[nvectors][size];

        // Load vectors.
        DataInputStream in = new DataInputStream(
                new FileInputStream(filename));

        for (int i = 0; i < nvectors; i++) {
            for (int j = 0; j < size; j++) {
                vectors[i][j] = in.readFloat();
            }
        }

        in.close();

        return vectors;
    }

    @Override
    public String getColumnName(int i) {
        return String.valueOf(i + 1);
    }

    public int getSize() {
        return data.size();
    }

    public int getRowCount() {
        return w;
    }

    public int getColumnCount() {
        return h;
    }

    @Override
    public Class getColumnClass(int i) {
        return getSize() > 0 ? getValueAt(0, 0).getClass() : Object.class;
    }

    public Object getValueAt(int row, int column) {
        int index = row * getColumnCount() + column;

        return data.get(index);
    }
}
