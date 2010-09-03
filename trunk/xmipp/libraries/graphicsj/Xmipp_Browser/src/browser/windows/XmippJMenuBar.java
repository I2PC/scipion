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
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.util.Vector;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

/**
 *
 * @author juanjo
 */
public class XmippJMenuBar extends JMenuBar implements ActionListener {

    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    protected ImagePlus imp;
    protected Vector<JMenu> menus = new Vector<JMenu>();
    // Extension menu to display menus that cannot be displayed in the menubar
    protected JMenu extension = new JMenu(">>");
    protected JMenu menuFlip = new JMenu(LABELS.OPERATION_MENU_FLIP);
    protected JMenuItem flipV = new JMenuItem(LABELS.OPERATION_FLIP_VERTICAL);
    protected JMenuItem flipH = new JMenuItem(LABELS.OPERATION_FLIP_HORIZONTAL);
    protected JMenu menuRotate = new JMenu(LABELS.OPERATION_MENU_ROTATE);
    protected JMenuItem rotateL = new JMenuItem(LABELS.OPERATION_ROTATE_LEFT);
    protected JMenuItem rotateR = new JMenuItem(LABELS.OPERATION_ROTATE_RIGHT);
    protected JMenu menuInfo = new JMenu(LABELS.OPERATION_MENU_INFO);
    protected JMenuItem histogram = new JMenuItem(LABELS.OPERATION_HISTOGRAM);
    protected JMenuItem plot_profile = new JMenuItem(LABELS.OPERATION_PLOT_PROFILE);
    protected JMenuItem measurements = new JMenuItem(LABELS.OPERATION_MEASUREMENTS);
    protected JMenu menuProcess = new JMenu(LABELS.OPERATION_MENU_PROCESS);
    protected JMenuItem brigth_and_contrast = new JMenuItem(LABELS.OPERATION_BC);
    protected JMenuItem threshold = new JMenuItem(LABELS.OPERATION_THRESHOLD);
    protected JMenuItem subtract_background = new JMenuItem(LABELS.OPERATION_SUBTRACTBG);
    protected JMenuItem gaussian_blur = new JMenuItem(LABELS.OPERATION_GAUSSIAN_BLUR);
    protected JMenuItem fft = new JMenuItem(LABELS.OPERATION_FFT);
    protected JMenuItem fft_band_pass_filter = new JMenuItem(LABELS.OPERATION_FFT_BAND_PASS_FILTER);
    protected JMenu menuReslice = new JMenu(LABELS.OPERATION_MENU_RESLICE);
    protected JMenuItem resliceTop = new JMenuItem(LABELS.OPERATION_RESLICE_TOP);
    protected JMenuItem resliceRight = new JMenuItem(LABELS.OPERATION_RESLICE_RIGHT);
    protected JMenu menuOpenAs = new JMenu(LABELS.OPERATION_MENU_OPEN_AS);
    protected JMenuItem openAs3D = new JMenuItem(LABELS.OPERATION_OPEN_AS_3D);
    protected JMenu menuGoTo = new JMenu(LABELS.OPERATION_MENU_GO_TO);
    protected JMenuItem goToSlice = new JMenuItem(LABELS.LABEL_GO_TO_SLICE);

    public XmippJMenuBar(ImagePlus imp) {
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

        addComponentListener(new ResizeListener());

        // Adds menus to list.
        menus.add(menuFlip);
        menus.add(menuRotate);
        menus.add(menuInfo);
        menus.add(menuProcess);
        menus.add(menuReslice);
        menus.add(menuOpenAs);
        menus.add(menuGoTo);
    }

    // Every time the menu is resized, check what menus can be displayed at their preferred size and move the rest to
    class ResizeListener extends ComponentAdapter {

        @Override
        public void componentResized(ComponentEvent evt) {
            XmippJMenuBar menuBar = XmippJMenuBar.this;

            // Clears menus.
            menuBar.removeAll();
            extension.removeAll();

            // Adds as much items as possible to the normal menu.
            int index = 0;
            while (index < menus.size()
                    && menuBar.getPreferredSize().width + extension.getPreferredSize().width
                    < menuBar.getWidth()) {

                menuBar.add(menus.get(index++));
            }

            if (index != menus.size()) {
                if (index > 0) {
                    menuBar.remove(--index);
                }
                menuBar.add(extension);
            }

            while (index < menus.size()) {
                extension.add(menus.get(index++));
            }

            menuBar.updateUI();
        }
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

    public static void main(String args[]) {
        JFrame frame = new JFrame("-- TEST --");
        ImagePlus ip = IJ.openImage("/home/juanjo/Descargas/NATURE-IsArutas_1024x768.jpg");
        XmippJMenuBar menuBar = new XmippJMenuBar(ip);
        frame.setLayout(new BorderLayout());
        frame.add(new JLabel(new ImageIcon(ip.getImage())), BorderLayout.CENTER);

        frame.setJMenuBar(menuBar);
        frame.setSize(400, 400);

        frame.setVisible(true);
    }
}
