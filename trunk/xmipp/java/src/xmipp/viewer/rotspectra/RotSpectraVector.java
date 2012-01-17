/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.rotspectra;

import java.awt.Image;
import java.util.ArrayList;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraVector {

    JFreeChart chart;
    String filename;
    String block;
    double vector[];
    ArrayList<String> images;
    protected boolean selected;

    public RotSpectraVector(String block, String filename, double vector[]) {
        this.block = block;
        this.filename = filename;
        this.vector = vector;

        chart = createChart(block, vector);
        images = loadFilenames(block + Filename.SEPARATOR + filename);
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
    }

    public boolean isSelected() {
        return selected;
    }

    static JFreeChart createChart(String label, double vector[]) {
        XYSeries series = new XYSeries(label);

//        double min = Double.MAX_VALUE;
//        double max = Double.MIN_VALUE;
        for (int i = 0; i < vector.length; i++) {
            series.add((double) i, vector[i]);
            /*            if (vector[i] > max) {
            max = vector[i];
            }
            if (vector[i] < min) {
            min = vector[i];
            }
             */
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "", "", "",
                dataset, PlotOrientation.VERTICAL,
                true, true, false);

//        chart.getXYPlot().getRangeAxis().setRange(min, max);
        chart.removeLegend();

        return chart;
    }

    static ArrayList<String> loadFilenames(String filename) {
        ArrayList<String> list = new ArrayList<String>();

        try {
            MetaData md = new MetaData(filename);

            long ids[] = md.findObjects();
            for (int i = 0; i < ids.length; i++) {
                list.add(md.getValueString(MDLabel.MDL_IMAGE, ids[i], true));
            }
        } catch (Exception ex) {
        }

        return list;
    }

    public Image getPreview(int w, int h) {
        return chart.createBufferedImage(w, h);
    }

    public JFrame getChart() {
        ChartPanel panel = new ChartPanel(chart);
        JFrame frame = new JFrame();
        frame.setTitle(block + Filename.SEPARATOR + filename);
        frame.getContentPane().add(panel);
        frame.pack();

        return frame;
    }

    public int getNImages() {
        return images.size();
    }

    public ArrayList<String> getImagesFilenames() {
        return images;
    }

    public String getTooltipText() {
        return block;
    }
}
