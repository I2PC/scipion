package browser.table.micrographs.ctf.profile;

import browser.LABELS;
import browser.imageitems.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Line;
import ij.gui.Plot;
import ij.gui.ProfilePlot;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.LinkedList;
import javax.swing.BoxLayout;
import javax.swing.JTabbedPane;
import xmipp.CTFDescription;
import xmipp.ImageDouble;

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
    private ImagePlus psdimage;
    private Scrollbar scrollbar;
    private CTFDescription ctfmodel;
    private JLabelPlot jlPlot, jlPlotAVG;
    private Point centerDisplay, centerPSD;
    private int displayLength, samples;
    private double angle;
    private Checkbox cbBgNoise, cbEnvelope, cbPSD, cbCTF;
    private double xValues[];
    private double plotProfile[], plotBgNoise[], plotEnvelope[], plotPSD[], plotCTF[];
    private double plotAVGprofile[], plotAVGbgnoise[], plotAVGenvelope[], plotAVGpsd[], plotAVGctf[];
//    private double minPlots, minAVGPlots, maxPlots, maxAVGPlots;
    private double CTFUpperLimit, CTFLowerLimit;
    private Plot plot, plotAVG;

    public CTFViewImageWindow(ImagePlus imp, String CTFFilename, String PSDFilename) {
        super(imp);

        setBackground(Color.WHITE);

        try {
            ImageCanvas imc = imp.getCanvas();
            imc.setEnabled(false);  // Clicks won't remove selection ;)

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
                }
            });

            // JLabels for plots.
            jlPlot = new JLabelPlot();
            jlPlotAVG = new JLabelPlot();

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

            JTabbedPane jtpPlots = new JTabbedPane();
            jtpPlots.add(LABELS.LABEL_TAB_PROFILE, jlPlot);
            jtpPlots.add(LABELS.LABEL_TAB_RADIAL_AVERAGE, jlPlotAVG);

            Panel pCenter = new Panel(new BorderLayout());
            pCenter.add(jtpPlots);

            // Check boxes for plots
            cbBgNoise = new Checkbox(LABELS.CB_PLOT_BGNOISE, true);
            cbEnvelope = new Checkbox(LABELS.CB_PLOT_ENVELOPE, true);
            cbPSD = new Checkbox(LABELS.CB_PLOT_PSD, true);
            cbCTF = new Checkbox(LABELS.CB_PLOT_CTF, true);

            cbBgNoise.setForeground(COLOR_BACKGROUND_NOISE);
            cbEnvelope.setForeground(COLOR_ENVELOPE);
            cbPSD.setForeground(COLOR_PSD);
            cbCTF.setForeground(COLOR_CTF);

            cbBgNoise.addItemListener(this);
            cbEnvelope.addItemListener(this);
            cbPSD.addItemListener(this);
            cbCTF.addItemListener(this);

            Panel panelCBoxes = new Panel();
            panelCBoxes.setLayout(new BoxLayout(panelCBoxes, BoxLayout.PAGE_AXIS));
            panelCBoxes.add(cbBgNoise);
            panelCBoxes.add(cbEnvelope);
            panelCBoxes.add(cbPSD);
            panelCBoxes.add(cbCTF);
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

            xValues = getXValues(samples);

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
        }

        updatePlots();
        updateRadialAverage();
    }

    private void setCBEnableStatus() {
        cbBgNoise.setEnabled(!cbCTF.getState());
        cbEnvelope.setEnabled(!cbCTF.getState());
        cbPSD.setEnabled(!cbCTF.getState());
    }

    /*
     * Completely refreshes the main profile plot: calculates new data + refreshes display.
     */
    private void refreshPlots() {
        calculatePlots();
        updatePlots();
    }

    private void updatePlots() {
        // Store it for "extract graph" operation
        plot = buildPlot(LABELS.LABEL_PROFILES, xValues, CTFLowerLimit, CTFUpperLimit,
                plotProfile,
                plotBgNoise,
                plotEnvelope,
                plotPSD,
                plotCTF,
                cbBgNoise.isEnabled() && cbBgNoise.getState(),
                cbEnvelope.isEnabled() && cbEnvelope.getState(),
                cbPSD.isEnabled() && cbPSD.getState(),
                cbCTF.isEnabled() && cbCTF.getState());

        jlPlot.setImage(plot.getImagePlus());
    }

    private void updateRadialAverage() {
        calculateRadialProfile();

        plotAVG = buildPlot(LABELS.LABEL_RADIAL_AVERAGE, xValues, CTFLowerLimit, CTFUpperLimit,
                plotAVGprofile,
                plotAVGbgnoise,
                plotAVGenvelope,
                plotAVGpsd,
                plotAVGctf,
                cbBgNoise.isEnabled() && cbBgNoise.getState(),
                cbEnvelope.isEnabled() && cbEnvelope.getState(),
                cbPSD.isEnabled() && cbPSD.getState(),
                cbCTF.isEnabled() && cbCTF.getState());

        jlPlotAVG.setImage(plotAVG.getImagePlus());
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
        CTFLowerLimit = min(plotCTF);
        CTFUpperLimit = max(plotCTF);
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
        CTFLowerLimit = min(plotAVGctf);
        CTFUpperLimit = max(plotAVGctf);
    }
    /*
    private double[] storePlotsMinAndMax(double profile[], double bg[], double envelope[], double psd[], double ctf[]) {
    double profileMin = min(profile);
    double bgMin = min(bg);
    double envelopeMin = min(envelope);
    double psdMin = min(psd);
    double ctfMin = min(ctf);

    double min = min(new double[]{profileMin, bgMin, envelopeMin, psdMin, ctfMin});

    double profileMax = max(profile);
    double bgMax = max(bg);
    double envelopeMax = max(envelope);
    double psdMax = max(psd);
    double ctfMax = max(ctf);

    double max = min(new double[]{profileMax, bgMax, envelopeMax, psdMax, ctfMax});

    return new double[]{min, max};
    }*/

    private static Plot buildPlot(String label, double xValues[], double lowerLimit, double upperLimit,
            double profile[], double bgNoise[], double envelope[], double psd[], double ctf[],
            boolean show_bgnoise, boolean show_envelope, boolean show_psd, boolean show_ctf) {

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

        String yAxis = show_ctf ? LABELS.LABEL_CTF : LABELS.LABEL_PSD;
        Plot plot = new Plot(label, LABELS.LABEL_RADIAL_XLEGEND, yAxis, xValues, plots.getFirst());

        if (show_ctf) {
            plot.setLimits(0, profile.length,
                    lowerLimit, upperLimit);
        }

        for (int i = 0; i < plots.size(); i++) {
            plot.setColor(colors.get(i));

            plot.addPoints(xValues, plots.get(i), Plot.LINE);
        }
        plot.setColor(colors.getFirst());

        return plot;
    }

    private static double[] getXValues(int size) {
        double values[] = new double[size];

        for (int i = 0; i < values.length; i++) {
            values[i] = i;
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

    private static double min(double array[]) {
        double min = array[0];

        for (int i = 1; i < array.length; i++) {
            double value = array[i];
            if (value < min) {
                min = value;
            }
        }

        return min;
    }

    private static double max(double array[]) {
        double max = array[0];

        for (int i = 1; i < array.length; i++) {
            double value = array[i];
            if (value > max) {
                max = value;
            }
        }

        return max;
    }
}
