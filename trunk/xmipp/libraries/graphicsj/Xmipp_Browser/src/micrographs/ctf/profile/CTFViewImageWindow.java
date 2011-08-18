package micrographs.ctf.profile;

import browser.LABELS;
import browser.imageitems.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Line;
import ij.gui.ProfilePlot;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Checkbox;
import java.awt.CheckboxGroup;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.LinkedList;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.UIManager;
import micrographs.ctf.profile.utils.Plot;
import xmipp.CTFDescription;
import xmipp.ImageDouble;
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
public class CTFViewImageWindow extends ImageWindow implements ItemListener {

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
    private JLabelPlot jlPlot, jlPlotAVG;
    private Point centerDisplay, centerPSD;
    private int displayLength, samples;
    private double angle;
    private Checkbox cbBgNoise, cbEnvelope, cbPSD, cbCTF;
    private CheckboxGroup cbgSampling = new CheckboxGroup();
    private Checkbox cbTdirect, cbTinverse;
    private double xValues[];
    private double plotProfile[], plotBgNoise[], plotEnvelope[], plotPSD[], plotCTF[];
    private double plotAVGprofile[], plotAVGbgnoise[], plotAVGenvelope[], plotAVGpsd[], plotAVGctf[];
    private Plot plot, plotAVG;
    private double Tm;
    private boolean invert = true;

    public CTFViewImageWindow(ImagePlus imp, String CTFFilename, String PSDFilename) {
        super(imp);

        setBackground(Color.WHITE);

        try {
            ImageCanvas imc = imp.getCanvas();
            imc.setEnabled(false);  // Clicks won't remove selection ;)

            Tm = getSamplingRate(CTFFilename);

            ctfmodel = new CTFDescription(CTFFilename);
            ImageDouble image = new ImageDouble(PSDFilename);

            psdimage = ImageConverter.convertToImagej(image, PSDFilename, false);

            removeAll();
            setLayout(new FlowLayout());

            scrollbar = new Scrollbar(Scrollbar.HORIZONTAL);
            scrollbar.setMaximum(MAX_VALUE + scrollbar.getVisibleAmount());
            scrollbar.addAdjustmentListener(new java.awt.event.AdjustmentListener() {

                public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
                    refreshPlots();
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
                    refreshPlots();

                    // Updates tooltips
                    Point p = e.getPoint();
                    jlPlot.setToolTipText(plot.getCoordinates(p.x, p.y));
                    jlPlotAVG.setToolTipText(plotAVG.getCoordinates(p.x, p.y));
                }
            });

            // Buttons to extract graphs.
            Button bExtractGraph = new Button(LABELS.BUTTON_EXTRACT_PROFILE);
            bExtractGraph.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    plot.show();
                }
            });

            Button bExtractAvgGraph = new Button(LABELS.BUTTON_EXTRACT_RADIAL_AVERAGE);
            bExtractAvgGraph.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    plotAVG.show();
                }
            });

            // JLabels for plots.
            jlPlot = new JLabelPlot();
            jlPlotAVG = new JLabelPlot();

            JTabbedPane jtpPlots = new JTabbedPane();
            JPanel panelPlot = new JPanel();
            JPanel panelPlotAVG = new JPanel();

            panelPlot.add(jlPlot);
            panelPlotAVG.add(jlPlotAVG);

            panelPlot.setBackground(Color.WHITE);
            panelPlotAVG.setBackground(Color.WHITE);

            jtpPlots.add(LABELS.LABEL_TAB_PROFILE, panelPlot);
            jtpPlots.add(LABELS.LABEL_TAB_RADIAL_AVERAGE, panelPlotAVG);

            Panel pCenter = new Panel(new BorderLayout());
            pCenter.add(jtpPlots);

            // Check boxes for plots
            cbBgNoise = new Checkbox(LABELS.CB_PLOT_BGNOISE, false);
            cbEnvelope = new Checkbox(LABELS.CB_PLOT_ENVELOPE, false);
            cbPSD = new Checkbox(LABELS.CB_PLOT_PSD, true);
            cbCTF = new Checkbox(LABELS.CB_PLOT_CTF, true);

            cbTdirect = new Checkbox(LABELS.SAMPLING_DIRECT, false, cbgSampling);
            cbTinverse = new Checkbox(LABELS.SAMPLING_INVERSE, true, cbgSampling);

            cbBgNoise.addItemListener(this);
            cbEnvelope.addItemListener(this);
            cbPSD.addItemListener(this);
            cbCTF.addItemListener(this);
            cbTdirect.addItemListener(this);
            cbTinverse.addItemListener(this);

            Panel panelCBoxes = new Panel();
            panelCBoxes.setLayout(new BoxLayout(panelCBoxes, BoxLayout.PAGE_AXIS));
            panelCBoxes.add(cbBgNoise);
            panelCBoxes.add(cbEnvelope);
            panelCBoxes.add(cbPSD);
            panelCBoxes.add(cbCTF);

            Panel panelSampling = new Panel();
            panelSampling.add(new Label(LABELS.LABEL_SAMPLING));
            panelSampling.add(cbTdirect);
            panelSampling.add(cbTinverse);

            panelCBoxes.add(panelSampling);
            panelCBoxes.add(bExtractGraph, BorderLayout.SOUTH);
            panelCBoxes.add(bExtractAvgGraph, BorderLayout.SOUTH);

            // Adds panels: image + Plot + Controls
            setLayout(new BorderLayout());
            add(panelImage, BorderLayout.WEST);
            add(pCenter, BorderLayout.CENTER);
            add(panelCBoxes, BorderLayout.EAST);

            // Store initial values.
            centerDisplay = new Point(imp.getWidth() / 2, imp.getHeight() / 2);
            centerPSD = new Point(psdimage.getWidth() / 2, psdimage.getHeight() / 2);

            displayLength = (int) (imp.getWidth() / 2);
            samples = (int) (psdimage.getWidth() / 2);

            xValues = getXValues(samples, Tm);

            setCBEnableStatus();
            refreshPlots();
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
        } else if (ie.getSource() == cbTinverse) {
            invert = true;
        } else if (ie.getSource() == cbTdirect) {
            invert = false;
        }

        updatePlots();
        updateRadialAverage();
    }

    private void setCBEnableStatus() {
        cbBgNoise.setEnabled(!cbCTF.getState());
        cbEnvelope.setEnabled(!cbCTF.getState());
        cbPSD.setEnabled(!cbCTF.getState());

        cbBgNoise.setForeground(cbBgNoise.isEnabled() ? COLOR_BACKGROUND_NOISE : COLOR_DISABLED);
        cbEnvelope.setForeground(cbEnvelope.isEnabled() ? COLOR_ENVELOPE : COLOR_DISABLED);
        cbPSD.setForeground(cbPSD.isEnabled() ? COLOR_PSD : COLOR_DISABLED);
        cbCTF.setForeground(cbCTF.isEnabled() ? COLOR_CTF : COLOR_DISABLED);
    }

    /*
     * Completely refreshes the main profile plot: calculates new data + refreshes display.
     */
    private void refreshPlots() {
        calculatePlots();
        updatePlots();
    }

    private boolean showBGNoise() {
        return cbBgNoise.isEnabled() && cbBgNoise.getState();
    }

    private boolean showEnvelope() {
        return cbEnvelope.isEnabled() && cbEnvelope.getState();
    }

    private boolean showPSD() {
        return cbPSD.isEnabled() && cbPSD.getState();
    }

    private boolean showCTF() {
        return cbCTF.isEnabled() && cbCTF.getState();
    }

    private void updatePlots() {
        boolean show_ctf = showCTF();

        // Store it for "extract graph" operation
        plot = buildPlot(LABELS.LABEL_PROFILES,
                invert ? LABELS.LABEL_RADIAL_FREQ_INVERSE : LABELS.LABEL_RADIAL_FREQ_DIRECT,
                show_ctf ? LABELS.LABEL_CTF : LABELS.LABEL_PSD,
                xValues,
                plotProfile, plotBgNoise, plotEnvelope, plotPSD, plotCTF,
                showBGNoise(), showEnvelope(), showPSD(),
                show_ctf, !invert);

        jlPlot.setImage(plot.getImagePlus());
        jlPlot.addMouseMotionListener(new MouseMotionListener() {

            public void mouseDragged(MouseEvent me) {
            }

            public void mouseMoved(MouseEvent me) {
                Point p = me.getPoint();
                jlPlot.setToolTipText(plot.getCoordinates(p.x, p.y));
            }
        });
    }

    private void updateRadialAverage() {
        calculateRadialProfile();

        boolean show_ctf = showCTF();
        plotAVG = buildPlot(LABELS.LABEL_RADIAL_AVERAGE,
                invert ? LABELS.LABEL_RADIAL_FREQ_INVERSE : LABELS.LABEL_RADIAL_FREQ_DIRECT,
                show_ctf ? LABELS.LABEL_CTF : LABELS.LABEL_PSD,
                xValues,
                plotAVGprofile, plotAVGbgnoise, plotAVGenvelope, plotAVGpsd, plotAVGctf,
                showBGNoise(), showEnvelope(), showPSD(),
                show_ctf, !invert);

        jlPlotAVG.setImage(plotAVG.getImagePlus());
        jlPlotAVG.addMouseMotionListener(new MouseMotionListener() {

            public void mouseDragged(MouseEvent me) {
            }

            public void mouseMoved(MouseEvent me) {
                Point p = me.getPoint();
                jlPlotAVG.setToolTipText(plotAVG.getCoordinates(p.x, p.y));
            }
        });
    }

//    private void setToolTipValues(Point p) {
//        System.out.println(System.currentTimeMillis() + " >> " + p + " > " + plot.getCoordinates(p.x, p.y));
//        jlPlot.setToolTipText(plot.getCoordinates(p.x, p.y));
//        jlPlotAVG.setToolTipText(plotAVG.getCoordinates(p.x, p.y));
//    }
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
//        psdimage.show();

        // Get profile.
        ProfilePlot profilePlot = new ProfilePlot(psdimage);
        plotProfile = profilePlot.getProfile();

        ctfmodel.CTFProfile(angle, samples);

        plotBgNoise = ctfmodel.profiles[CTFDescription.BACKGROUND_NOISE];
        plotEnvelope = ctfmodel.profiles[CTFDescription.ENVELOPE];
        plotPSD = ctfmodel.profiles[CTFDescription.PSD];
        plotCTF = ctfmodel.profiles[CTFDescription.CTF];

//        storePlotsMinAndMax(plotProfile, plotBgNoise, plotEnvelope, plotPSD, plotCTF);
//        CTFLowerLimit = min(plotCTF);
//        CTFUpperLimit = max(plotCTF);
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

    private void calculateRadialProfile() {
        plotAVGprofile = new double[samples];

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

        plotAVGbgnoise = ctfmodel.avgprofiles[CTFDescription.BACKGROUND_NOISE];
        plotAVGenvelope = ctfmodel.avgprofiles[CTFDescription.ENVELOPE];
        plotAVGpsd = ctfmodel.avgprofiles[CTFDescription.PSD];
        plotAVGctf = ctfmodel.avgprofiles[CTFDescription.CTF];

//        storePlotsMinAndMax(plotAVGprofile, plotAVGbgnoise, plotAVGenvelope, plotAVGpsd, plotAVGctf);
//        CTFLowerLimit = min(plotAVGctf);
//        CTFUpperLimit = max(plotAVGctf);
    }

    private static Plot buildPlot(String label, String xLabel, String yLabel,
            double xValues[],
            //double lowerLimit, double upperLimit,
            double profile[], double bgNoise[], double envelope[], double psd[],
            double ctf[], boolean show_bgnoise, boolean show_envelope,
            boolean show_psd, boolean show_ctf, boolean invert) {

        LinkedList<double[]> plots = new LinkedList<double[]>();
        LinkedList<Color> colors = new LinkedList<Color>();

        if (show_ctf) {
            plots.add(ctf);
            colors.add(COLOR_CTF);
        } else {
            plots.add(profile);
            colors.add(COLOR_PROFILE);

            if (show_bgnoise) {
                plots.add(bgNoise);
                colors.add(COLOR_BACKGROUND_NOISE);
            }

            if (show_envelope) {
                plots.add(envelope);
                colors.add(COLOR_ENVELOPE);
            }

            if (show_psd) {
                plots.add(psd);
                colors.add(COLOR_PSD);
            }
        }

        Plot plot = new Plot(label, xLabel, yLabel, xValues, plots.getFirst());
        plot.setInvert(invert);

        for (int i = 1; i < plots.size(); i++) {
            plot.setColor(colors.get(i));

            plot.addPoints(xValues, plots.get(i), Plot.LINE);
        }
        plot.setColor(colors.getFirst());

        return plot;
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
//
//    private static double min(double array[]) {
//        double min = array[0];
//
//        for (int i = 1; i < array.length; i++) {
//            double value = array[i];
//            if (value < min) {
//                min = value;
//            }
//        }
//
//        return min;
//    }
//
//    private static double max(double array[]) {
//        double max = array[0];
//
//        for (int i = 1; i < array.length; i++) {
//            double value = array[i];
//            if (value > max) {
//                max = value;
//            }
//        }
//
//        return max;
//    }
}
