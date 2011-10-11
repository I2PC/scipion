package micrographs;

import browser.DEBUG;
import browser.LABELS;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import xmipp.MDLabel;
import xmipp.MetaData;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class FSCWindow extends JFrame {

    public FSCWindow(String filename) {
        super(LABELS.TITLE_FSC + filename);

        try {
            setLayout(new BorderLayout());

            MetaData mdout = new MetaData();
            mdout.computeFourierStatistics(filename);

            double xValues[] = mdout.getColumnValues(MDLabel.MDL_RESOLUTION_FREQ);
            double y1s[] = mdout.getColumnValues(MDLabel.MDL_RESOLUTION_FRC);
            double y2s[] = mdout.getColumnValues(MDLabel.MDL_RESOLUTION_DPR);

            //String labelX = MetaData.label2Str(MDLabel.MDL_RESOLUTION_FREQ);
            String labelY1 = MetaData.label2Str(MDLabel.MDL_RESOLUTION_FRC);
            String labelY2 = MetaData.label2Str(MDLabel.MDL_RESOLUTION_DPR);

            XYSeriesCollection dataset1 = new XYSeriesCollection(
                    createSeries(labelY1, xValues, y1s));
            XYSeriesCollection dataset2 = new XYSeriesCollection(
                    createSeries(labelY2, xValues, y2s));

            JFreeChart chart = ChartFactory.createXYLineChart(
                    "", LABELS.LABEL_DIGITAL_FREQUENCY, "",
                    dataset1, PlotOrientation.VERTICAL,
                    true, true, false);

            // Sets second axis.
            XYPlot plot = chart.getXYPlot();
            ValueAxis axisDPR = plot.getRangeAxis();
            axisDPR.setLabel(labelY1);

            NumberAxis axisFSC = new NumberAxis();
            axisFSC.setLabel(labelY2);
            axisFSC.setAutoRangeIncludesZero(false);

            plot.setRangeAxis(1, axisFSC);
            plot.setDataset(1, dataset2);
            plot.mapDatasetToRangeAxis(1, 1);

            // Sets a new renderer, so series color will be different.
            plot.setRenderer(1, new StandardXYItemRenderer());

            add(createChartPanel(chart), BorderLayout.CENTER);

            Button bOk = new Button(LABELS.BUTTON_OK);
            bOk.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent ae) {
                    dispose();
                }
            });

            Panel pButtons = new Panel();
            pButtons.add(bOk);
            add(pButtons, BorderLayout.SOUTH);

            pack();
            setLocationRelativeTo(null);
        } catch (Exception e) {
            DEBUG.printException(e);
        }
    }

    private static ChartPanel createChartPanel(JFreeChart chart) {
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setDomainPannable(true);
        plot.setRangePannable(true);

        java.util.List list = Arrays.asList(new Integer[]{
                    new Integer(0), new Integer(1)
                });
        plot.mapDatasetToDomainAxes(0, list);
        plot.mapDatasetToRangeAxes(0, list);
        ChartUtilities.applyCurrentTheme(chart);

        return new ChartPanel(chart);
    }

    private static XYSeries createSeries(String name, double[] xs, double[] values) {
        XYSeries series = new XYSeries(name);

        for (int i = 0; i < xs.length; i++) {
            series.add(xs[i], values[i]);
        }

        return series;
    }
}
