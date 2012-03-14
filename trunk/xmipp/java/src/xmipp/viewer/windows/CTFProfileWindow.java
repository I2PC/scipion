package xmipp.viewer.windows;

import xmipp.utils.WindowUtil;
import xmipp.utils.XmippLabel;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Line;
import ij.gui.ProfilePlot;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.UIManager;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import xmipp.jni.CTFDescription;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.ij.commons.XmippImageConverter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
@SuppressWarnings("serial")
public class CTFProfileWindow extends ImageWindow implements ItemListener, ActionListener {

    private final static BasicStroke plotsStroke = new BasicStroke(2.0f);
    private final static int MAX_VALUE = 360;
    private final static Color COLOR_PROFILE = Color.BLACK;
    private final static Color COLOR_BACKGROUND_NOISE = Color.RED;
    private final static Color COLOR_ENVELOPE = Color.GREEN;
    private final static Color COLOR_PSD = Color.BLUE;
    private final static Color COLOR_CTF = Color.MAGENTA;
    private final static Color COLOR_DISABLED = UIManager.getColor("ComboBox.disabledForeground");
    private ImagePlus psdimage;
    private Scrollbar scrollbar;
    private CTFDescription ctfmodel;
    private Point centerDisplay, centerPSD;
    private int displayLength, samples;
    private double angle;
    private JCheckBox cbBgNoise, cbEnvelope, cbPSD, cbCTF;
    private JPanel panelPlot, panelPlotAVG;
    private XYSeriesCollection datasetProfile, datasetAVG;
    private ChartPanel chartPanelProfile, chartPanelAVG;
    private JFileChooser fc = new JFileChooser();
    private JButton bExportProfile, bExportAVG;
    private double xValues[];
    private double Tm;

    public CTFProfileWindow(ImagePlus imp, String CTFFilename, String PSDFilename) {
        super(imp);

        try {
            ImageCanvas imc = imp.getCanvas();
            imc.setEnabled(false);  // Clicks won't remove selection ;)

            Tm = getSamplingRate(CTFFilename);

            ctfmodel = new CTFDescription(CTFFilename);
            ImageGeneric image = new ImageGeneric(PSDFilename);
            image.setUseLogarithm(false);

            psdimage = XmippImageConverter.readToImagePlus(image);
            psdimage.setTitle(PSDFilename);

            removeAll();
            GridBagConstraints gc = new GridBagConstraints();
            setLayout(new GridBagLayout());

            scrollbar = new Scrollbar(Scrollbar.HORIZONTAL);
            scrollbar.setMaximum(MAX_VALUE + scrollbar.getVisibleAmount());
            scrollbar.addAdjustmentListener(new java.awt.event.AdjustmentListener() {

                public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
                    updatePlots();
                }
            });

            Panel panelImage = new Panel(new BorderLayout());
            panelImage.add(imc, BorderLayout.CENTER);
            panelImage.add(scrollbar, BorderLayout.SOUTH);
            addMouseWheelListener(new MouseWheelListener() {

                public void mouseWheelMoved(MouseWheelEvent e) {
                    int MAXSCROLL = scrollbar.getMaximum() - scrollbar.getVisibleAmount();

                    int newValue = (scrollbar.getValue() + e.getUnitsToScroll())
                            % MAXSCROLL;

                    if (newValue < 0) {
                        newValue = newValue + MAXSCROLL;
                    }

                    scrollbar.setValue(newValue);
                    updatePlots();
                }
            });

            JTabbedPane jtpPlots = new JTabbedPane();
            panelPlot = new JPanel();
            panelPlotAVG = new JPanel();

            jtpPlots.add(XmippLabel.LABEL_TAB_PROFILE, panelPlot);
            jtpPlots.add(XmippLabel.LABEL_TAB_RADIAL_AVERAGE, panelPlotAVG);

            Panel pCenter = new Panel(new BorderLayout());
            pCenter.add(jtpPlots);

            // Check boxes for plots
            cbBgNoise = new JCheckBox(XmippLabel.CB_PLOT_BGNOISE, false);
            cbEnvelope = new JCheckBox(XmippLabel.CB_PLOT_ENVELOPE, false);
            cbPSD = new JCheckBox(XmippLabel.CB_PLOT_PSD, true);
            cbCTF = new JCheckBox(XmippLabel.CB_PLOT_CTF, true);

            cbBgNoise.addItemListener(this);
            cbEnvelope.addItemListener(this);
            cbPSD.addItemListener(this);
            cbCTF.addItemListener(this);

            bExportProfile = WindowUtil.getTextButton(XmippLabel.BUTTON_EXPORT_PROFILE, this);
            bExportAVG = WindowUtil.getTextButton(XmippLabel.BUTTON_EXPORT_RADIAL_AVERAGE, this);

            JPanel panelCBoxes = new JPanel(new GridBagLayout());
            //BoxLayout layout = new BoxLayout(panelCBoxes, BoxLayout.PAGE_AXIS);
            //layout.
            

            gc.anchor = GridBagConstraints.WEST;
            panelCBoxes.add(cbBgNoise, WindowUtil.getConstraints(gc, 0, 0));
            panelCBoxes.add(cbEnvelope, WindowUtil.getConstraints(gc, 0, 1));
            panelCBoxes.add(cbPSD, WindowUtil.getConstraints(gc, 0, 2));
            panelCBoxes.add(cbCTF, WindowUtil.getConstraints(gc, 0, 3));
            gc.weightx = 1.0;
            gc.insets = new Insets(0, 0, 5, 5);
            panelCBoxes.add(bExportProfile, WindowUtil.getConstraints(gc, 0, 4, 2));
            panelCBoxes.add(bExportAVG, WindowUtil.getConstraints(gc, 0, 5, 2));

            // Adds panels.
            //setLayout(new BorderLayout());
            add(panelImage, WindowUtil.getConstraints(gc, 0, 0));
            add(pCenter, WindowUtil.getConstraints(gc, 1, 0, 2, 2));
            gc.anchor = GridBagConstraints.EAST;
            add(panelCBoxes,WindowUtil.getConstraints(gc, 0, 1));
//          add(panelImage);
//          add(pCenter);
//          add(panelCBoxes);           

            // Store initial values.
            centerDisplay = new Point(imp.getWidth() / 2, imp.getHeight() / 2);
            centerPSD = new Point(psdimage.getWidth() / 2, psdimage.getHeight() / 2);

            displayLength = (int) (imp.getWidth() / 2);
            samples = (int) (psdimage.getWidth() / 2);

            xValues = getXValues(samples, Tm);

            setCBEnableStatus();
            updatePlots();
            updateRadialAverage();

            // Trick: setting no resizable allows to pack properly.
            setResizable(false);
            pack();
            //setResizable(true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * When checkbox changes, plots will be redisplayed, but data is the same,
     * so there is no need to recalculate anything.
     */
    public void itemStateChanged(ItemEvent ie) {
        if (ie.getSource() == cbCTF) {
            setCBEnableStatus();
        }

        // Restores general view.
        chartPanelProfile.restoreAutoBounds();
        chartPanelAVG.restoreAutoBounds();

        updatePlots();
        updateRadialAverage();
    }

    public void actionPerformed(ActionEvent ae) {
        if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            if (file != null) {
                JFreeChart chart = null;
                if (ae.getSource() == bExportProfile) {
                    chart = chartPanelProfile.getChart();
                } else if (ae.getSource() == bExportAVG) {
                    chart = chartPanelAVG.getChart();
                }

                if (export(file.getAbsolutePath(), chart)) {
                    IJ.showMessage("Save", "File saved sucesfully.");
                }
            }
        }
    }

    private void setCBEnableStatus() {
        cbBgNoise.setEnabled(!cbCTF.isEnabled());
        cbEnvelope.setEnabled(!cbCTF.isEnabled());
        cbPSD.setEnabled(!cbCTF.isEnabled());

        cbBgNoise.setForeground(cbBgNoise.isEnabled() ? COLOR_BACKGROUND_NOISE : COLOR_DISABLED);
        cbEnvelope.setForeground(cbEnvelope.isEnabled() ? COLOR_ENVELOPE : COLOR_DISABLED);
        cbPSD.setForeground(cbPSD.isEnabled() ? COLOR_PSD : COLOR_DISABLED);
        cbCTF.setForeground(cbCTF.isEnabled() ? COLOR_CTF : COLOR_DISABLED);
    }

    private boolean showBGNoise() {
        return cbBgNoise.isEnabled() && cbBgNoise.isEnabled();
    }

    private boolean showEnvelope() {
        return cbEnvelope.isEnabled() && cbEnvelope.isEnabled();
    }

    private boolean showPSD() {
        return cbPSD.isEnabled() && cbPSD.isEnabled();
    }

    private boolean showCTF() {
        return cbCTF.isEnabled() && cbCTF.isEnabled();
    }

    private void updatePlots() {
        calculatePlots();

        boolean show_ctf = showCTF();

        if (chartPanelProfile != null) {
            XYPlot plot = ((XYPlot) chartPanelProfile.getChart().getPlot());
            plot.setDataset(0, datasetProfile);
            plot.getRangeAxis().setLabel(show_ctf ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD);
        } else {
//            JFreeChart chart = createChart("",
//                    LABELS.LABEL_SAMPLING,
//                    show_ctf ? LABELS.LABEL_CTF : LABELS.LABEL_PSD,
//                    datasetProfile);

            JFreeChart chart = ChartFactory.createXYLineChart(
                    "", XmippLabel.LABEL_SAMPLING,
                    show_ctf ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD,
                    datasetProfile,
                    PlotOrientation.VERTICAL,
                    true, true, false);

            chartPanelProfile = createChartPanel(chart);

//            chartPanelProfile.setPreferredSize(new Dimension(chartPanelProfile.getPreferredSize().width, getImagePlus().getHeight()));

            panelPlot.add(chartPanelProfile);
//            pack();
        }

        customizeSeriesRenderes((XYPlot) chartPanelProfile.getChart().getPlot(),
                showBGNoise(), showEnvelope(), showPSD(), show_ctf);
    }

    private void updateRadialAverage() {
        calculateRadialAverageProfile();

        boolean show_ctf = showCTF();

        if (chartPanelAVG != null) {
            XYPlot plot = ((XYPlot) chartPanelAVG.getChart().getPlot());
            plot.setDataset(0, datasetAVG);
            plot.getRangeAxis().setLabel(show_ctf ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD);
        } else {
//            JFreeChart chart = createChart("",
//                    LABELS.LABEL_SAMPLING,
//                    show_ctf ? LABELS.LABEL_CTF : LABELS.LABEL_PSD,
//                    datasetAVG);
            JFreeChart chart = ChartFactory.createXYLineChart(
                    "",
                    XmippLabel.LABEL_SAMPLING,
                    show_ctf ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD,
                    datasetAVG, PlotOrientation.VERTICAL,
                    true, true, false);

            chartPanelAVG = createChartPanel(chart);

//            chartPanelAVG.setPreferredSize(new Dimension(chartPanelAVG.getPreferredSize().width, getImagePlus().getHeight()));

            panelPlotAVG.add(chartPanelAVG);
//            pack();
        }

        customizeSeriesRenderes((XYPlot) chartPanelAVG.getChart().getPlot(),
                showBGNoise(), showEnvelope(), showPSD(), show_ctf);
    }

    private static ChartPanel createChartPanel(JFreeChart chart) {
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setDomainPannable(true);
        plot.setRangePannable(true);

        List<Integer> list = Arrays.asList(new Integer[]{
                    new Integer(0), new Integer(1)
                });
        plot.mapDatasetToDomainAxes(0, list);
        plot.mapDatasetToRangeAxes(0, list);
        ChartUtilities.applyCurrentTheme(chart);

        return new ChartPanel(chart);
    }

//    private static JFreeChart createChart(String title, String xLabel, String yLabel, XYDataset dataset) {
//        JFreeChart chart = ChartFactory.createXYLineChart(
//                title, xLabel, yLabel, dataset, PlotOrientation.VERTICAL,
//                true, true, false);
//
//        XYPlot plot = (XYPlot) chart.getPlot();
//        plot.setDomainPannable(true);
//        plot.setRangePannable(true);
//
//        java.util.List list = Arrays.asList(new Integer[]{
//                    new Integer(0), new Integer(1)
//                });
//        plot.mapDatasetToDomainAxes(0, list);
//        plot.mapDatasetToRangeAxes(0, list);
//        ChartUtilities.applyCurrentTheme(chart);
//
//        return chart;
//    }
    private void customizeSeriesRenderes(XYPlot plot, boolean show_bgnoise, boolean show_envelope,
            boolean show_psd, boolean show_ctf) {

        XYItemRenderer renderer = plot.getRenderer();

        if (show_ctf) {
            renderer.setSeriesPaint(0, COLOR_CTF);
            plot.getRenderer().setSeriesStroke(0, plotsStroke);
        } else {
            // 0: Profile.
            renderer.setSeriesPaint(0, COLOR_PROFILE);
            plot.getRenderer().setSeriesStroke(0, plotsStroke);

            int index = 1;
            // 1: BGNoise.
            if (show_bgnoise) {
                renderer.setSeriesPaint(index, COLOR_BACKGROUND_NOISE);
                plot.getRenderer().setSeriesStroke(index, plotsStroke);
                index++;
            }

            // 2: BGNoise.
            if (show_envelope) {
                renderer.setSeriesPaint(index, COLOR_ENVELOPE);
                plot.getRenderer().setSeriesStroke(index, plotsStroke);
                index++;
            }

            // 3: PSD.
            if (show_psd) {
                renderer.setSeriesPaint(index, COLOR_PSD);
                plot.getRenderer().setSeriesStroke(index, plotsStroke);
                index++;
            }
        }
    }

    private static XYSeriesCollection createSeriesCollection(double[] xs,
            double profile[], double bgNoise[], double envelope[], double psd[],
            double ctf[], boolean show_bgnoise, boolean show_envelope,
            boolean show_psd, boolean show_ctf) {

        XYSeriesCollection collection = new XYSeriesCollection();

        if (show_ctf) {
            collection.addSeries(createSeries(XmippLabel.CB_PLOT_CTF, xs, ctf));
        } else {
            collection.addSeries(createSeries(XmippLabel.CB_PLOT_PROFILE, xs, profile));

            if (show_bgnoise) {
                collection.addSeries(createSeries(XmippLabel.CB_PLOT_BGNOISE, xs, bgNoise));
            }
            if (show_envelope) {
                collection.addSeries(createSeries(XmippLabel.CB_PLOT_ENVELOPE, xs, envelope));
            }
            if (show_psd) {
                collection.addSeries(createSeries(XmippLabel.CB_PLOT_PSD, xs, psd));
            }
        }

        return collection;
    }

    private static XYSeries createSeries(String name, double[] xs, double[] values) {
        XYSeries series = new XYSeries(name);

        for (int i = 0; i < xs.length; i++) {
            series.add(xs[i], values[i]);
        }

        return series;
    }

    private double getSamplingRate(String CTFFilename) {
        double samplingRate = 1.0;

        try {
            MetaData md = new MetaData(CTFFilename);
            samplingRate = md.getValueDouble(MDLabel.MDL_CTF_SAMPLING_RATE, md.firstObject());
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return samplingRate;
    }

    private void calculatePlots() {
        double value = scrollbar.getValue();
        angle = Math.toRadians(-value);

        // Makes line roi...
        int xDisplay = centerDisplay.x + (int) (displayLength * Math.cos(angle));
        int yDisplay = centerDisplay.y + (int) (displayLength * Math.sin(angle));

        int xPSD = centerPSD.x + (int) (samples * Math.cos(angle));
        int yPSD = centerPSD.y + (int) (samples * Math.sin(angle));

        // Over display...
        //IJ.selectWindow();
        IJ.makeLine(centerDisplay.x, centerDisplay.y, xDisplay, yDisplay);

        // ...and internal PSD.
        Line line = new Line(centerPSD.x, centerPSD.y, xPSD, yPSD);
        psdimage.setRoi(line);

        // Get profile.
        ProfilePlot profilePlot = new ProfilePlot(psdimage);
        double[] plotProfile = profilePlot.getProfile();

        ctfmodel.CTFProfile(angle, samples);

        double[] plotBgNoise = ctfmodel.profiles[CTFDescription.BACKGROUND_NOISE];
        double[] plotEnvelope = ctfmodel.profiles[CTFDescription.ENVELOPE];
        double[] plotPSD = ctfmodel.profiles[CTFDescription.PSD];
        double[] plotCTF = ctfmodel.profiles[CTFDescription.CTF];

        datasetProfile = createSeriesCollection(
                xValues, plotProfile, plotBgNoise, plotEnvelope, plotPSD, plotCTF,
                showBGNoise(), showEnvelope(), showPSD(), showCTF());
    }

    private void calculateRadialAverageProfile() {
        double[] plotAVGprofile = new double[samples];

        for (int i = 0; i < MAX_VALUE; i++) {
            double current_angle = Math.toRadians(-i);

            // Make line roi.
            int xPSD = centerPSD.x + (int) (samples * Math.cos(current_angle));
            int yPSD = centerPSD.y + (int) (samples * Math.sin(current_angle));

            Line line = new Line(centerPSD.x, centerPSD.y, xPSD, yPSD);
            psdimage.setRoi(line);

            // Get profile.
            ProfilePlot profilePlot = new ProfilePlot(psdimage);
            double profile[] = profilePlot.getProfile();

            // Total summatory.
            sum(plotAVGprofile, profile);
        }

        // Get profile average.
        div(plotAVGprofile, MAX_VALUE);

        // Get rest of profiles average.
        ctfmodel.CTFAverageProfiles(samples);

        double[] plotAVGbgnoise = ctfmodel.avgprofiles[CTFDescription.BACKGROUND_NOISE];
        double[] plotAVGenvelope = ctfmodel.avgprofiles[CTFDescription.ENVELOPE];
        double[] plotAVGpsd = ctfmodel.avgprofiles[CTFDescription.PSD];
        double[] plotAVGctf = ctfmodel.avgprofiles[CTFDescription.CTF];

        datasetAVG = createSeriesCollection(
                xValues, plotAVGprofile, plotAVGbgnoise, plotAVGenvelope, plotAVGpsd, plotAVGctf,
                showBGNoise(), showEnvelope(), showPSD(), showCTF());
    }

    private boolean export(String path, JFreeChart chart) {
        XYPlot plot = (XYPlot) chart.getPlot();
        XYDataset dataset = plot.getDataset();
        int series = dataset.getSeriesCount();
        int items = dataset.getItemCount(0);

        try {
            FileWriter fstream = new FileWriter(path);
            BufferedWriter out = new BufferedWriter(fstream);

            // Header.
            out.write("X\t");
            for (int i = 0; i < series; i++) {
                out.write(dataset.getSeriesKey(i) + "\t");
                //System.out.print(dataset.getSeriesKey(i) + "\t");
            }
            out.newLine();

            for (int item = 0; item < items; item++) {
                out.write(dataset.getXValue(0, item) + "\t");
                for (int serie = 0; serie < series; serie++) {
                    out.write(dataset.getYValue(serie, item) + "\t");
                }
                out.newLine();
            }
            out.close();

            return true;
        } catch (Exception ex) {
            //ex.printStackTrace();
            IJ.error(ex.getMessage());
        }

        return false;
    }

    private static double[] getXValues(int length, double Tm) {
        double values[] = new double[length];
        double factor = 1 / (length * 2 * Tm);

        for (int i = 0; i < length; i++) {
            values[i] = (double) i * factor;
        }

        return values;
    }

    private static void sum(double a[], double b[]) {
        for (int i = 0; i < a.length; i++) {
            a[i] += b[i];
        }
    }

    private static void div(double array[], double value) {
        for (int i = 0; i < array.length; i++) {
            array[i] /= value;
        }
    }
}
