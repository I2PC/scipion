/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rotspectra;

import java.awt.Image;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraVector {

    JFreeChart chart;
    String filename;
    String block;
    double vector[];

    public RotSpectraVector(String block, String filename, double vector[]) {
        this.block = block;
        this.filename = filename;
        this.vector = vector;

        chart = createChart(block, vector);
    }

    static JFreeChart createChart(String label, double vector[]) {
        XYSeries series = new XYSeries(label);

        for (int i = 0; i < vector.length; i++) {
            series.add(i, vector[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "", "", "",
                dataset, PlotOrientation.VERTICAL,
                true, true, false);
        chart.removeLegend();

        return chart;
    }

    public Image getPlotImage(int w, int h) {
        return chart.createBufferedImage(w, h);
    }
}
