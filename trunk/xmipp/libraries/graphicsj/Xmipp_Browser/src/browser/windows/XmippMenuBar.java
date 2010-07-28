/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.windows;

import browser.LABELS;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 *
 * @author juanjo
 */
public class XmippMenuBar extends MenuBar implements ActionListener {

    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    protected ImagePlus imp;
    protected Menu menuFlip = new Menu(LABELS.OPERATION_MENU_FLIP);
    protected MenuItem flipV = new MenuItem(LABELS.OPERATION_FLIP_VERTICAL);
    protected MenuItem flipH = new MenuItem(LABELS.OPERATION_FLIP_HORIZONTAL);
    protected Menu menuRotate = new Menu(LABELS.OPERATION_MENU_ROTATE);
    protected MenuItem rotateL = new MenuItem(LABELS.OPERATION_ROTATE_LEFT);
    protected MenuItem rotateR = new MenuItem(LABELS.OPERATION_ROTATE_RIGHT);
    protected Menu menuInfo = new Menu(LABELS.OPERATION_MENU_INFO);
    protected MenuItem histogram = new MenuItem(LABELS.OPERATION_HISTOGRAM);
    protected MenuItem plot_profile = new MenuItem(LABELS.OPERATION_PLOT_PROFILE);
    protected MenuItem measurements = new MenuItem(LABELS.OPERATION_MEASUREMENTS);
    protected Menu menuProcess = new Menu(LABELS.OPERATION_MENU_PROCESS);
    protected MenuItem brigth_and_contrast = new MenuItem(LABELS.OPERATION_BC);
    protected MenuItem threshold = new MenuItem(LABELS.OPERATION_THRESHOLD);
    protected MenuItem subtract_background = new MenuItem(LABELS.OPERATION_SUBTRACTBG);
    protected MenuItem gaussian_blur = new MenuItem(LABELS.OPERATION_GAUSSIAN_BLUR);
    protected MenuItem fft = new MenuItem(LABELS.OPERATION_FFT);
    protected MenuItem fft_band_pass_filter = new MenuItem(LABELS.OPERATION_FFT_BAND_PASS_FILTER);
    protected Menu menuReslice = new Menu(LABELS.OPERATION_MENU_RESLICE);
    protected MenuItem resliceTop = new MenuItem(LABELS.OPERATION_RESLICE_TOP);
    protected MenuItem resliceRight = new MenuItem(LABELS.OPERATION_RESLICE_RIGHT);
    protected Menu menuOpenAs = new Menu(LABELS.OPERATION_MENU_OPEN_AS);
    protected MenuItem openAs3D = new MenuItem(LABELS.OPERATION_OPEN_AS_3D);
    protected Menu menuGoTo = new Menu(LABELS.OPERATION_MENU_GO_TO);
    protected MenuItem goToSlice = new MenuItem(LABELS.LABEL_GO_TO_SLICE);

    public XmippMenuBar(ImagePlus imp) {
        super();

        this.imp = imp;

        createMenuBar();

        menuReslice.setEnabled(imp.getStackSize() > 1);
        menuOpenAs.setEnabled(imp.getStackSize() > 1);
        menuGoTo.setEnabled(imp.getStackSize() > 1);
    }

    private void createMenuBar() {
        // Flip.
        menuFlip.add(flipV);
        menuFlip.add(flipH);

        flipV.addActionListener(this);
        flipH.addActionListener(this);

        // Rotate.
        menuRotate.add(rotateL);
        menuRotate.add(rotateR);

        rotateL.addActionListener(this);
        rotateR.addActionListener(this);

        // Info.
        menuInfo.add(histogram);
        menuInfo.add(plot_profile);
        menuInfo.add(measurements);

        histogram.addActionListener(this);
        plot_profile.addActionListener(this);
        measurements.addActionListener(this);

        // Process.
        menuProcess.add(brigth_and_contrast);
        menuProcess.add(threshold);
        menuProcess.add(subtract_background);
        menuProcess.add(gaussian_blur);
        menuProcess.add(fft);
        menuProcess.add(fft_band_pass_filter);

        brigth_and_contrast.addActionListener(this);
        threshold.addActionListener(this);
        subtract_background.addActionListener(this);
        gaussian_blur.addActionListener(this);
        fft.addActionListener(this);
        fft_band_pass_filter.addActionListener(this);

        // Reslice.
        menuReslice.add(resliceTop);
        menuReslice.add(resliceRight);

        resliceTop.addActionListener(this);
        resliceRight.addActionListener(this);

        // Open As.
        menuOpenAs.add(openAs3D);

        openAs3D.addActionListener(this);

        // Go to.
        menuGoTo.add(goToSlice);

        goToSlice.addActionListener(this);

        // Adds menus.
        add(menuFlip);
        add(menuRotate);
        add(menuInfo);
        add(menuProcess);
        add(menuReslice);
        add(menuOpenAs);
        add(menuGoTo);
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == flipV) {
            imp.getProcessor().flipVertical();
        } else if (e.getSource() == flipH) {
            imp.getProcessor().flipHorizontal();
        } else if (e.getSource() == rotateL) {
            IJ.run("Rotate 90 Degrees Left");
        } else if (e.getSource() == rotateR) {
            IJ.run("Rotate 90 Degrees Right");
        } else if (e.getSource() == histogram) {
            IJ.run("Histogram");
        } else if (e.getSource() == plot_profile) {
            // Check if there is a line selected, so , if not, it
            // selects the right tool for the next time.
            Roi roi = imp.getRoi();
            if (roi == null || !roi.isLine()) {
                IJ.setTool(Toolbar.LINE);
            }

            IJ.run("Plot Profile");
        } else if (e.getSource() == measurements) {
            IJ.run("Set Measurements...");
            IJ.run("Measure");
        } else if (e.getSource() == brigth_and_contrast) {
            IJ.run("Brightness/Contrast...");
        } else if (e.getSource() == threshold) {
            IJ.run("Threshold...");
        } else if (e.getSource() == subtract_background) {
            IJ.run("Subtract Background...");
        } else if (e.getSource() == gaussian_blur) {
            IJ.run("Gaussian Blur...");
        } else if (e.getSource() == fft) {
            IJ.run("FFT");
        } else if (e.getSource() == fft_band_pass_filter) {
            IJ.run("Bandpass Filter...");//, "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");
        } else if (e.getSource() == resliceTop) {
            IJ.run("Reslice [/]...", "slice=1.000 start=Top");
        } else if (e.getSource() == resliceRight) {
            IJ.run("Reslice [/]...", "slice=1.000 start=Right");
        } else if (e.getSource() == openAs3D) {
            run3DViewer(imp);
        } else if (e.getSource() == goToSlice) {
            GenericDialog dialog = new GenericDialog(LABELS.TITLE_GO_TO_SLICE);
            dialog.addNumericField("Slice: ", imp.getCurrentSlice(), String.valueOf(imp.getStackSize()).length());
            dialog.showDialog();

            int slice = (int) dialog.getNextNumber();
            if (slice < 1) {
                slice = 1;
            } else if (slice > imp.getStackSize()) {
                slice = imp.getStackSize();
            }

            imp.setSlice(slice);
        }

        imp.updateAndDraw();
    }

    private void run3DViewer(ImagePlus ip) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);

        // Adds the sphere image plus to universe.
        new StackConverter(ip).convertToRGB();
        Content c = universe.addVoltex(ip);
        c.displayAs(Content.VOLUME);

        universe.show();    // Shows...
    }
}
