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
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

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
    protected Menu menuGoTo = new Menu(LABELS.OPERATION_MENU_GO_TO);
    protected MenuItem goToSlice = new MenuItem(LABELS.LABEL_GO_TO_SLICE);
    protected Menu menuOpenAs = new Menu(LABELS.OPERATION_MENU_OPEN_AS);
    protected MenuItem openAs3D = new MenuItem(LABELS.OPERATION_OPEN_AS_3D);
    protected Menu menuPlugins = new Menu(LABELS.OPERATION_MENU_PLUGINS);
    protected Menu menuAnalysis = new Menu(LABELS.OPERATION_MENU_ANALYSIS);
    protected MenuItem itemAutocorrelation = new MenuItem("Autocorrelation");
    protected MenuItem itemLineAnalyzer = new MenuItem("Line Analyzer");
    protected MenuItem itemOvalProfilePlot = new MenuItem("Oval Profile Plot");
    protected MenuItem itemRadialProfilePlotAngle = new MenuItem("Radial Profile Plot Angle");
    protected MenuItem itemRadialProfilePlotHeight = new MenuItem("Radial Profile Plot Height");
    protected MenuItem itemSyncMeasure3D = new MenuItem("Sync Measure 3D");
    protected MenuItem itemTemplateMatching = new MenuItem("Template Matching");
    protected MenuItem itemConcentricCircles = new MenuItem("Concentric Circles (non-destructive overlay)");
    protected MenuItem itemAnisotropicDiffusion = new MenuItem("Anisotropic Diffusion (edge-preserving noise reduction)");
    protected MenuItem itemMultiOtsuThreshold = new MenuItem("Multi Otsu Threshold");
    protected Menu menuConvert = new Menu(LABELS.OPERATION_MENU_CONVERT);
    protected MenuItem itemFloatMorphology = new MenuItem("Float Morphology");
    protected MenuItem itemMeanShiftFilter = new MenuItem("Mean Shift Filter (edge-preserving smoothing)");
    protected MenuItem itemCLAHE = new MenuItem("CLAHE (Contrast Limited Adaptive Histogram Equalization)");
    protected Menu menuThreshold = new Menu(LABELS.OPERATION_MENU_THRESHOLD);
    protected MenuItem itemMixtureModelingThresholding = new MenuItem("Mixture Modeling Thresholding");
    protected MenuItem itemOtsuThresholding = new MenuItem("Otsu Thresholding");
    protected MenuItem itemWatershed = new MenuItem("Watershed");
    protected MenuItem itemMaximumEntropyThresholding = new MenuItem("Maximum Entropy Thresholding");
    protected MenuItem itemMultiThresholder = new MenuItem("MultiThresholder");
    protected MenuItem itemMultiOtsuThreshold_ = new MenuItem("Multi Otsu Threshold");
    protected MenuItem itemSIOX = new MenuItem("SIOX (Simple Interactive Object Extraction)");
    protected MenuItem itemRATS = new MenuItem("RATS (Robust Automatic Threshold Selection)");
    protected Menu menuOther = new Menu(LABELS.OPERATION_MENU_OTHER);
    protected MenuItem itemContourPlotter = new MenuItem("Contour Plotter");
    protected MenuItem itemDottedandDashedLines = new MenuItem("Dotted and Dashed Lines");
    protected MenuItem itemRadialGrid = new MenuItem("Radial Grid");
    protected MenuItem itemPolarTransformer = new MenuItem("Polar Transformer (convert to polar coordinates)");
    protected MenuItem itemVoxelCounter = new MenuItem("Voxel Counter");
    protected MenuItem itemSubstackMaker = new MenuItem("Substack Maker");
    protected MenuItem itemVolumeViewer3DSlicer = new MenuItem("Volume Viewer/3D-Slicer");
    protected MenuItem item3DViewer = new MenuItem("3D Viewer");
    protected MenuItem itemRadialReslice = new MenuItem("Radial Reslice (orthogonal reconstructions)");
    protected MenuItem itemStackSlicer = new MenuItem("StackSlicer (dynamic orthoganal views)");
    protected MenuItem itemUntiltStack = new MenuItem("Untilt Stack (removes tilt from a stack of images)");
    protected MenuItem itemStraightenCurvedObjects = new MenuItem("Straighten Curved Objects");
    protected MenuItem itemSyncWindows = new MenuItem("Sync Windows");
    protected MenuItem itemVolumeJ = new MenuItem("VolumeJ");
    protected MenuItem itemFlowJ = new MenuItem("FlowJ");
    protected MenuItem itemFlow3J = new MenuItem("Flow3J");
    protected MenuItem itemSurfaceJ = new MenuItem("SurfaceJ");
    protected MenuItem itemPCA = new MenuItem("PCA (Principal Component Analysis)");
    protected MenuItem itemTransformJ = new MenuItem("TransformJ");
    protected MenuItem itemFeatureJ = new MenuItem("FeatureJ");
    protected MenuItem itemObjectJ = new MenuItem("ObjectJ");

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

        // Plugins
        menuPlugins.add(menuAnalysis);
        menuAnalysis.add(itemAnisotropicDiffusion);
        menuAnalysis.add(itemAutocorrelation);
        menuAnalysis.add(itemConcentricCircles);
        menuAnalysis.add(itemLineAnalyzer);
        menuAnalysis.add(itemMultiOtsuThreshold);
        menuAnalysis.add(itemOvalProfilePlot);
        menuAnalysis.add(itemRadialProfilePlotAngle);
        menuAnalysis.add(itemRadialProfilePlotHeight);
        menuAnalysis.add(itemSyncMeasure3D);
        menuAnalysis.add(itemTemplateMatching);

        menuPlugins.add(menuConvert);
        menuConvert.add(itemCLAHE);
        menuConvert.add(itemFloatMorphology);
        menuConvert.add(itemMeanShiftFilter);

        menuPlugins.add(menuThreshold);
        menuThreshold.add(itemMaximumEntropyThresholding);
        menuThreshold.add(itemMixtureModelingThresholding);
        menuThreshold.add(itemMultiThresholder);
        menuThreshold.add(itemMultiOtsuThreshold_);
        menuThreshold.add(itemOtsuThresholding);
        menuThreshold.add(itemRATS);
        menuThreshold.add(itemSIOX);
        menuThreshold.add(itemWatershed);

        menuPlugins.add(menuOther);
        menuOther.add(itemContourPlotter);
        menuOther.add(itemDottedandDashedLines);
        menuOther.add(itemFeatureJ);
        menuOther.add(itemFlowJ);
        menuOther.add(itemFlow3J);
        menuOther.add(itemObjectJ);
        menuOther.add(itemPCA);
        menuOther.add(itemPolarTransformer);
        menuOther.add(itemRadialGrid);
        menuOther.add(itemRadialReslice);
        menuOther.add(itemStackSlicer);
        menuOther.add(itemStraightenCurvedObjects);
        menuOther.add(itemSubstackMaker);
        menuOther.add(itemSurfaceJ);
        menuOther.add(itemSyncWindows);
        menuOther.add(itemTransformJ);
        menuOther.add(itemUntiltStack);
        menuOther.add(itemVolumeJ);
        menuOther.add(itemVolumeViewer3DSlicer);
        menuOther.add(itemVoxelCounter);
        menuOther.add(item3DViewer);

        // ACtion listeners
        itemAutocorrelation.addActionListener(this);
        itemLineAnalyzer.addActionListener(this);
        itemOvalProfilePlot.addActionListener(this);
        itemRadialProfilePlotAngle.addActionListener(this);
        itemRadialProfilePlotHeight.addActionListener(this);
        itemSyncMeasure3D.addActionListener(this);
        itemTemplateMatching.addActionListener(this);
        itemConcentricCircles.addActionListener(this);
        itemAnisotropicDiffusion.addActionListener(this);
        itemMultiOtsuThreshold.addActionListener(this);

        itemFloatMorphology.addActionListener(this);
        itemMeanShiftFilter.addActionListener(this);
        itemCLAHE.addActionListener(this);

        itemMixtureModelingThresholding.addActionListener(this);
        itemOtsuThresholding.addActionListener(this);
        itemWatershed.addActionListener(this);
        itemMaximumEntropyThresholding.addActionListener(this);
        itemMultiThresholder.addActionListener(this);
        itemMultiOtsuThreshold_.addActionListener(this);
        itemSIOX.addActionListener(this);
        itemRATS.addActionListener(this);

        itemContourPlotter.addActionListener(this);
        itemDottedandDashedLines.addActionListener(this);
        itemRadialGrid.addActionListener(this);
        itemPolarTransformer.addActionListener(this);
        itemVoxelCounter.addActionListener(this);
        itemSubstackMaker.addActionListener(this);
        itemVolumeViewer3DSlicer.addActionListener(this);
        item3DViewer.addActionListener(this);
        itemRadialReslice.addActionListener(this);
        itemStackSlicer.addActionListener(this);
        itemUntiltStack.addActionListener(this);
        itemStraightenCurvedObjects.addActionListener(this);
        itemSyncWindows.addActionListener(this);
        itemVolumeJ.addActionListener(this);
        itemFlowJ.addActionListener(this);
        itemFlow3J.addActionListener(this);
        itemSurfaceJ.addActionListener(this);
        itemPCA.addActionListener(this);
        itemTransformJ.addActionListener(this);
        itemFeatureJ.addActionListener(this);
        itemObjectJ.addActionListener(this);

        // Adds menus.
        add(menuFlip);
        add(menuRotate);
        add(menuInfo);
        add(menuProcess);
        add(menuReslice);
        add(menuOpenAs);
        add(menuGoTo);
        add(menuPlugins);
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
        } else if (e.getSource() == itemAnisotropicDiffusion) {
            IJ.run("Anisotropic Diffusion...");
        } else if (e.getSource() == itemAutocorrelation) {
            IJ.run("Auto Corr");
        } else if (e.getSource() == itemConcentricCircles) {
            IJ.run("Concentric Circles");
        } else if (e.getSource() == itemLineAnalyzer) {
            IJ.run("Line Analyzer");    // ?
        } else if (e.getSource() == itemMultiOtsuThreshold) {
            IJ.run("Otsu Thresholding");
        } else if (e.getSource() == itemOvalProfilePlot) {
            IJ.run("Oval Profile");
        } else if (e.getSource() == itemRadialProfilePlotAngle) {
            IJ.run("Radial Profile Angle");
        } else if (e.getSource() == itemRadialProfilePlotHeight) {
            IJ.run("Radial Profile Height");
        } else if (e.getSource() == itemSyncMeasure3D) {
            IJ.run("Sync Measure 3D");
        } else if (e.getSource() == itemTemplateMatching) {
            IJ.run("Create Template");
        } else if (e.getSource() == itemCLAHE) {
            IJ.run("CLAHE ");
        } else if (e.getSource() == itemFloatMorphology) {
            IJ.run("Float Morphology");
        } else if (e.getSource() == itemMeanShiftFilter) {
            IJ.run("Mean Shift");
        } else if (e.getSource() == itemMaximumEntropyThresholding) {
            IJ.run("Entropy Threshold");
        } else if (e.getSource() == itemMixtureModelingThresholding) {
            IJ.run("Mixture Modeling");
        } else if (e.getSource() == itemMultiThresholder) {
            IJ.run("MultiThresholder");
        } else if (e.getSource() == itemMultiOtsuThreshold_) {
            IJ.run("Multi OtsuThreshold");
        } else if (e.getSource() == itemOtsuThresholding) {
            IJ.run("Otsu Thresholding");
        } else if (e.getSource() == itemRATS) {
            IJ.run("RATS ");
        } else if (e.getSource() == itemSIOX) {
            IJ.run("SIOX Segmentation");
        } else if (e.getSource() == itemWatershed) {
            IJ.run("Watershed");
        } else if (e.getSource() == itemContourPlotter) {
            IJ.run("ContourPlotter ");
        } else if (e.getSource() == itemDottedandDashedLines) {
            IJ.run("Dotted Line");
        } else if (e.getSource() == itemFeatureJ) {
            IJ.run("FeatureJ...");
        } else if (e.getSource() == itemFlowJ) {
            IJ.run("FlowJ");
        } else if (e.getSource() == itemFlow3J) {
            IJ.run("Flow3J");
        } else if (e.getSource() == itemObjectJ) {
            IJ.run("ObjectJ");
        } else if (e.getSource() == itemPCA) {
            IJ.run("PCA ");
        } else if (e.getSource() == itemPolarTransformer) {
            IJ.run("Polar Transformer");
        } else if (e.getSource() == itemRadialGrid) {
            IJ.run("Radial Grid");
        } else if (e.getSource() == itemRadialReslice) {
            IJ.run("Radial Reslice");
        } else if (e.getSource() == itemStackSlicer) {
            IJ.run("StackSlicer ");
        } else if (e.getSource() == itemStraightenCurvedObjects) {
            IJ.run("Straighten...");
        } else if (e.getSource() == itemSubstackMaker) {
            IJ.run("Substack Maker");
        } else if (e.getSource() == itemSurfaceJ) {
            IJ.run("SurfaceJ ");
        } else if (e.getSource() == itemSyncWindows) {
            IJ.run("Sync Windows");
        } else if (e.getSource() == itemTransformJ) {
            IJ.run("TransformJ...");
        } else if (e.getSource() == itemUntiltStack) {
            IJ.run("Untilt Stack");
        } else if (e.getSource() == itemVolumeJ) {
            IJ.run("VolumeJ ");
        } else if (e.getSource() == itemVolumeViewer3DSlicer) {
            IJ.run("Volume Viewer");
        } else if (e.getSource() == itemVoxelCounter) {
            IJ.run("Voxel Counter");
        } else if (e.getSource() == item3DViewer) {
            IJ.run("3D Viewer");
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
        ImagePlus ip = IJ.openImage("/home/juanjo/.gnome2/gnome-art/thumbnails/X-Trans-Th.jpg");
        XmippMenuBar menuBar = new XmippMenuBar(ip);
        frame.setLayout(new BorderLayout());
        frame.add(new JLabel(new ImageIcon(ip.getImage())), BorderLayout.CENTER);

        frame.setMenuBar(menuBar);
        frame.setSize(400, 400);

        frame.setVisible(true);
    }
}
