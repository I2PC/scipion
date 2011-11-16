/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.table.AbstractTableModel;
import xmipp.ImageGeneric;
import xmippij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableModel extends AbstractTableModel {

    protected ArrayList<ScoreItem> data = new ArrayList<ScoreItem>();
    protected ImagePlus focusedItem;
    int good, bad;

    public ImagesTableModel() {
        super();
    }

    public void clear() {
        data.clear();
    }

    public void addItem(ScoreItem scoreItem) {
        data.add(scoreItem);

        // Counts "good" and "bad" items.
        if (scoreItem.good) {
            good++;
        } else {
            bad++;
        }

        fireTableStructureChanged();
    }

    @Override
    public String getColumnName(int column) {
        return "" + column;
    }

    public int getGoodCount() {
        return good;
    }

    public int getBadCount() {
        return bad;
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

    public ArrayList<ScoreItem> getData() {
        return data;
    }

    public void sortItems() {
        Collections.sort(data);
    }

    public static ImagePlus mean(ArrayList<String> images) {
        try {
            ImageGeneric firstImage = new ImageGeneric();

            firstImage.readHeader(images.get(0));
            int w = firstImage.xSize;
            int h = firstImage.ySize;
            float mean[][] = new float[w][h];

            String currentFileName;

            // For all images...
            for (int k = 0; k < images.size(); k++) {
                ImageGeneric currentImage = new ImageGeneric();
                currentFileName = images.get(k);

                currentImage.readData(currentFileName);
                ImagePlus imp = XmippImageConverter.convertToImagej(currentImage);
                imp.setTitle(currentFileName);

                // TODO: Sum efficiently.
                // Adds current image to sum.
                for (int j = 0; j < h; j++) {
                    for (int i = 0; i < w; i++) {
                        mean[i][j] += imp.getProcessor().getf(i, j);
                    }
                }

                imp.close();
            }

            // Calculates mean...
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    mean[i][j] /= images.size(); // ...by dividing with the #images
                }
            }

            FloatProcessor processor = new FloatProcessor(mean);
            ImagePlus meanIp = new ImagePlus();
            meanIp.setProcessor("", processor);

            return meanIp;
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}
