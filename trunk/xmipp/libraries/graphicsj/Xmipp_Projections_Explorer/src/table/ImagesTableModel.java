/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.util.Vector;
import javax.swing.table.AbstractTableModel;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableModel extends AbstractTableModel {

    protected Vector<ScoreItem> data = new Vector<ScoreItem>();
    protected ImagePlus focusedItem;

    public ImagesTableModel() {
        super();
    }

    public void clear() {
        data.clear();
    }

    public void addItem(ScoreItem scoreItem) {
        data.add(scoreItem);

        fireTableStructureChanged();
    }

    @Override
    public String getColumnName(int column) {
        return "" + column;
    }

    public int getRowCount() {
        return 1;
    }

    public int getColumnCount() {
        return data.size();
    }

    public Object getValueAt(int rowIndex, int columnIndex) {
        return columnIndex < data.size() ? data.get(columnIndex) : null;
    }

    @Override
    public Class getColumnClass(int c) {
        Object object = getValueAt(0, c);

        return object != null ? object.getClass() : null;
    }

    public Vector<ScoreItem> getData() {
        return data;
    }

    public static ImagePlus mean(String images[]) {
        ImagePlus current = IJ.openImage(images[0]);
        int w = current.getWidth();
        int h = current.getHeight();
        float sum[][] = new float[w][h];

        // For all images...
        for (int k = 1; k < images.length; k++) {
            // Adds current image to sum.
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    sum[i][j] += current.getProcessor().getPixelValue(i, j);
                }
            }

            current.close();
            current = IJ.openImage(images[k]);
        }

        float mean[][] = new float[w][h];

        // Calculates mean...
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                mean[i][j] = sum[i][j] / images.length; // ...by dividing with the #images
            }
        }

        FloatProcessor processor = new FloatProcessor(mean);
        ImagePlus meanIp = new ImagePlus();
        meanIp.setProcessor("", processor);

        return meanIp;
    }
}
