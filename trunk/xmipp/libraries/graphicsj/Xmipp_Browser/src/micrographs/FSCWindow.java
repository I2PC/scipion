package micrographs;

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
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
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
//
//    public static void main(String args[]) {
//        String filename = "/home/juanjo/temp/angles_save.sel";
//
//        try {
//            MetaData md = new MetaData(filename);
//            FSCWindow frame = new FSCWindow(md, filename);
//
//            frame.setVisible(true);
//        } catch (Exception ex) {
//            ex.printStackTrace();
//        }
//    }

    public FSCWindow(MetaData md, String filename) {
        super("FSC: " + md.getFilename());

        try {
            setLayout(new BorderLayout());

            md.computeFourierStatistics(filename);

            double xValues[] = md.getColumnValues(MDLabel.MDL_RESOLUTION_FREQ);
            double y1s[] = md.getColumnValues(MDLabel.MDL_RESOLUTION_FRC);
            double y2s[] = md.getColumnValues(MDLabel.MDL_RESOLUTION_DPR);

            String labelX = MetaData.label2Str(MDLabel.MDL_RESOLUTION_FREQ);
            String labelY1 = MetaData.label2Str(MDLabel.MDL_RESOLUTION_FRC);
            String labelY2 = MetaData.label2Str(MDLabel.MDL_RESOLUTION_DPR);

            XYDataset dataset = createSeriesCollection(xValues, y1s, y2s, labelY1, labelY2);

            JFreeChart chart = ChartFactory.createXYLineChart(
                    "", "?", labelX,
                    dataset, PlotOrientation.VERTICAL,
                    true, true, false);

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
            e.printStackTrace();
        }
    }
//
//    public void actionPerformed(ActionEvent ae) {
//        if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
//            File file = fc.getSelectedFile();
//            if (file != null) {
//                JFreeChart chart = null;
//                if (ae.getSource() == bExportProfile) {
//                    chart = chartPanelProfile.getChart();
//                } else if (ae.getSource() == bExportAVG) {
//                    chart = chartPanelAVG.getChart();
//                }
//
//                if (export(file.getAbsolutePath(), chart)) {
//                    IJ.showMessage("Save", "File saved sucesfully.");
//                }
//            }
//        }
//    }

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
//    private void customizeSeriesRenderes(XYPlot plot, boolean show_bgnoise, boolean show_envelope,
//            boolean show_psd, boolean show_ctf) {
//
//        XYItemRenderer renderer = plot.getRenderer();
//
//        if (show_ctf) {
//            renderer.setSeriesPaint(0, COLOR_CTF);
//            plot.getRenderer().setSeriesStroke(0, plotsStroke);
//        } else {
//            // 0: Profile.
//            renderer.setSeriesPaint(0, COLOR_PROFILE);
//            plot.getRenderer().setSeriesStroke(0, plotsStroke);
//
//            int index = 1;
//            // 1: BGNoise.
//            if (show_bgnoise) {
//                renderer.setSeriesPaint(index, COLOR_BACKGROUND_NOISE);
//                plot.getRenderer().setSeriesStroke(index, plotsStroke);
//                index++;
//            }
//
//            // 2: BGNoise.
//            if (show_envelope) {
//                renderer.setSeriesPaint(index, COLOR_ENVELOPE);
//                plot.getRenderer().setSeriesStroke(index, plotsStroke);
//                index++;
//            }
//
//            // 3: PSD.
//            if (show_psd) {
//                renderer.setSeriesPaint(index, COLOR_PSD);
//                plot.getRenderer().setSeriesStroke(index, plotsStroke);
//                index++;
//            }
//        }
//    }

    private static XYSeriesCollection createSeriesCollection(
            double[] xs, double ys1[], double ys2[],
            String y1Label, String y2Label) {
        XYSeriesCollection collection = new XYSeriesCollection();

        collection.addSeries(createSeries(y1Label, xs, ys1));
        collection.addSeries(createSeries(y2Label, xs, ys2));

        return collection;
    }

    private static XYSeries createSeries(String name, double[] xs, double[] values) {
        XYSeries series = new XYSeries(name);

        for (int i = 0; i < xs.length; i++) {
            series.add(xs[i], values[i]);
        }

        return series;
    }
}
