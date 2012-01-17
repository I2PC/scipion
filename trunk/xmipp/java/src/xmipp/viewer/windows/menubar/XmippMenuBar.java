/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package xmipp.viewer.windows.menubar;

import xmipp.utils.LABELS;
import xmipp.viewer.windows.ImageWindowOperations;
import xmipp.viewer.windows.ImagesWindowFactory;
import xmipp.viewer.windows.StackWindowOperations;
import xmipp.viewer.windows.iPollImageWindow;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.io.FileInfo;
import ij.plugin.filter.Duplicater;
import java.awt.CheckboxMenuItem;
import java.awt.Container;
import java.awt.Menu;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class XmippMenuBar extends DynamicMenuBar implements ActionListener {

    private ImagePlus imp;
    boolean converted = false;
    boolean binarized = false;
    // *** Transform ***
    private Menu menuTransform = new Menu(LABELS.OPERATION_MENU_TRANSFORM);
    private MenuItem itemToTable = new MenuItem(LABELS.OPERATION_TO_GALLERY);
    private Menu menuFlip = new Menu(LABELS.OPERATION_MENU_FLIP);
    private MenuItem itemFlipV = new MenuItem(LABELS.OPERATION_FLIP_VERTICAL);
    private MenuItem itemFlipH = new MenuItem(LABELS.OPERATION_FLIP_HORIZONTAL);
    private Menu menuRotate = new Menu(LABELS.OPERATION_MENU_ROTATE);
    private MenuItem itemRotateL = new MenuItem(LABELS.OPERATION_ROTATE_LEFT);
    private MenuItem itemRotateR = new MenuItem(LABELS.OPERATION_ROTATE_RIGHT);
    private MenuItem itemCrop = new MenuItem(LABELS.OPERATION_CROP);
    private Menu menuReslice = new Menu(LABELS.LABEL_RESLICE);
    private MenuItem itemResliceTop = new MenuItem(LABELS.OPERATION_RESLICE_TOP);
    private MenuItem itemResliceRight = new MenuItem(LABELS.OPERATION_RESLICE_RIGHT);
    private MenuItem itemStackSlicer = new MenuItem(LABELS.OPERATION_STACK_SLICER);
    private MenuItem itemRadialReslice = new MenuItem(LABELS.OPERATION_RADIAL_RESLICE);
    private MenuItem itemStraightenCurvedObjects = new MenuItem(LABELS.OPERATION_STRAIGHTEN_CURVED_OBJECTS);
    private MenuItem itemUntiltStack = new MenuItem(LABELS.OPERATION_UNTILT_STACK);
    private Menu menuTransformJ = new Menu(LABELS.OPERATION_TRANSFORMJ);
    private MenuItem itemTJAbout = new MenuItem(LABELS.OPERATION_TJABOUT);
    private MenuItem itemTJAffine = new MenuItem(LABELS.OPERATION_TJAFFINE);
    private MenuItem itemTJCrop = new MenuItem(LABELS.OPERATION_TJCROP);
    private MenuItem itemTJEmbed = new MenuItem(LABELS.OPERATION_TJEMBED);
    private MenuItem itemTJMirror = new MenuItem(LABELS.OPERATION_TJMIRROR);
    private MenuItem itemTJOptions = new MenuItem(LABELS.OPERATION_TJOPTIONS);
    private MenuItem itemTJPanel = new MenuItem(LABELS.OPERATION_TJPANEL);
    private MenuItem itemTJRotate = new MenuItem(LABELS.OPERATION_TJROTATE);
    private MenuItem itemTJScale = new MenuItem(LABELS.OPERATION_TJSCALE);
    private MenuItem itemTJShift = new MenuItem(LABELS.OPERATION_TJSHIFT);
    private MenuItem itemTJTurn = new MenuItem(LABELS.OPERATION_TJTURN);
    private MenuItem itemTJWebsite = new MenuItem(LABELS.OPERATION_TJWEBSITE);
    private MenuItem itemSurfaceJ = new MenuItem(LABELS.OPERATION_SURFACEJ);
    // *** Denoising ***
    private Menu menuDenoising = new Menu(LABELS.OPERATION_MENU_DENOISING);
    private MenuItem itemFFT_band_pass_filter = new MenuItem(LABELS.OPERATION_FFT_BAND_PASS_FILTER);
    private MenuItem itemAnisotropicDiffusion = new MenuItem(LABELS.OPERATION_ANISOTROPIC_DIFFUSION);
    private MenuItem itemMeanShiftFilter = new MenuItem(LABELS.OPERATION_MEAN_SHIFT_FILTER);
    // *** Open with ***
    private Menu menuOpenWith = new Menu(LABELS.OPERATION_MENU_OPEN_WITH);
    private MenuItem item3DViewer = new MenuItem(LABELS.OPERATION_OPEN_AS_3D);
    private MenuItem itemVolumeViewer3DSlicer = new MenuItem(LABELS.OPERATION_VOLUME_VIEWER_3D_SLICER);
    private MenuItem itemVolumeJ = new MenuItem(LABELS.OPERATION_VOLUMEJ);
    // *** Info ***
    private Menu menuInfo = new Menu(LABELS.OPERATION_MENU_INFO);
    private CheckboxMenuItem itemPoll = new CheckboxMenuItem(LABELS.OPERATION_POLL);
    private MenuItem itemDuplicate = new MenuItem(LABELS.OPERATION_DUPLICATE);
    private MenuItem itemGoToSlice = new MenuItem(LABELS.OPERATION_GO_TO_SLICE);
    private MenuItem itemSubstackMaker = new MenuItem(LABELS.OPERATION_SUBSTACK_MAKER);
    private MenuItem itemHistogram = new MenuItem(LABELS.OPERATION_HISTOGRAM);
    private MenuItem itemPlotProfile = new MenuItem(LABELS.OPERATION_PLOT_PROFILE);
    private MenuItem itemMeasurements = new MenuItem(LABELS.OPERATION_MEASUREMENTS);
    private MenuItem itemAutocorrelation = new MenuItem(LABELS.OPERATION_AUTOCORRELATION);
    private MenuItem itemConcentricCircles = new MenuItem(LABELS.OPERATION_CONCENTRIC_CIRCLES);
    private MenuItem itemLineAnalyzer = new MenuItem(LABELS.OPERATION_LINE_ANALYZER);
    private MenuItem itemOvalProfilePlot = new MenuItem(LABELS.OPERATION_OVAL_PROFILER_PLOT);
    private MenuItem itemRadialProfilePlotAngle = new MenuItem(LABELS.OPERATION_RADIAL_PROFILE_PLOT_ANGLE);
    private MenuItem itemRadialProfilePlotHeight = new MenuItem(LABELS.OPERATION_RADIAL_PROFILE_PLOT_HEIGHT);
    private MenuItem itemContourPlotter = new MenuItem(LABELS.OPERATION_CONTOUR_PLOTTER);
    private MenuItem itemFeatureJ = new MenuItem(LABELS.OPERATION_FEATUREJ);
    private MenuItem itemSyncWindows = new MenuItem(LABELS.OPERATION_SYNC_WINDOWS);
    private MenuItem itemSyncMeasure3D = new MenuItem(LABELS.OPERATION_SYNC_MEASURE_3D);
    // *** Thresholding ***
    private Menu menuThresholding = new Menu(LABELS.OPERATION_MENU_THRESHOLDING);
    private MenuItem itemThreshold = new MenuItem(LABELS.OPERATION_THRESHOLD);
    private MenuItem itemOtsuThresholding = new MenuItem(LABELS.OPERATION_OTSU_THRESHOLDING);
    private MenuItem itemMultiOtsuThreshold = new MenuItem(LABELS.OPERATION_MULTI_OTSU_THRESHOLDING);
    private MenuItem itemMixtureModelingThresholding = new MenuItem(LABELS.OPERATION_MIXTURE_MODELING_THRESHOLDING);
    private MenuItem itemMaximumEntropyThresholding = new MenuItem(LABELS.OPERATION_MAXIMUM_ENTROPY_THRESHOLDING);
    private MenuItem itemMultiThresholder = new MenuItem(LABELS.OPERATION_MULTI_THRESHOLDER);
    private MenuItem itemSIOX = new MenuItem(LABELS.OPERATION_SIOX);
    private MenuItem itemRATS = new MenuItem(LABELS.OPERATION_RATS);
    // *** Binary ***
    private Menu menuBinary = new Menu(LABELS.OPERATION_MENU_BINARY);
    private MenuItem itemVoxelCounter = new MenuItem(LABELS.OPERATION_VOXEL_COUNTER);
    private MenuItem itemErode = new MenuItem(LABELS.OPERATION_ERODE);
    private MenuItem itemDilate = new MenuItem(LABELS.OPERATION_DILATE);
    private MenuItem itemOpen = new MenuItem(LABELS.OPERATION_OPEN);
    private MenuItem itemClose = new MenuItem(LABELS.OPERATION_CLOSE);
    private MenuItem itemFloatMorphology = new MenuItem(LABELS.OPERATION_FLOAT_MORPHOLOGY);
    private MenuItem itemOutline = new MenuItem(LABELS.OPERATION_OUTLINE);
    private MenuItem itemFillHoles = new MenuItem(LABELS.OPERATION_FILL_HOLES);
    private MenuItem itemSkeletonize = new MenuItem(LABELS.OPERATION_SKELETONIZE);
    private MenuItem itemDistanceMap = new MenuItem(LABELS.OPERATION_DISTANCE_MAP);
    private MenuItem itemUltimatePoints = new MenuItem(LABELS.OPERATION_ULTIMATE_POINTS);
    private MenuItem itemWatershed = new MenuItem(LABELS.OPERATION_WATERSHED);
    private MenuItem itemVoronoi = new MenuItem(LABELS.OPERATION_VORONOI);
    private MenuItem itemOptions = new MenuItem(LABELS.OPERATION_OPTIONS);
    // *** Process ***
    private Menu menuProcess = new Menu(LABELS.OPERATION_MENU_PROCESS);
    private MenuItem itemBC = new MenuItem(LABELS.OPERATION_BC);
    private MenuItem itemEnhanceContrast = new MenuItem(LABELS.OPERATION_ENHANCE_CONTRAST);
    private MenuItem itemSubtractBG = new MenuItem(LABELS.OPERATION_SUBTRACTBG);
    private MenuItem itemGaussianBlur = new MenuItem(LABELS.OPERATION_GAUSSIAN_BLUR);
    private MenuItem itemConvolve = new MenuItem(LABELS.OPERATION_CONVOLVE);
    private MenuItem itemMedian = new MenuItem(LABELS.OPERATION_MEDIAN);
    private MenuItem itemMean = new MenuItem(LABELS.OPERATION_MEAN);
    private MenuItem itemFFT = new MenuItem(LABELS.OPERATION_FFT);
    private MenuItem itemCLAHE = new MenuItem(LABELS.OPERATION_CLAHE);
    private MenuItem itemPCA = new MenuItem(LABELS.OPERATION_PCA);
    private MenuItem itemPolarTransformer = new MenuItem(LABELS.OPERATION_POLAR_TRANSFORMER);
    private MenuItem itemTemplateMatching = new MenuItem(LABELS.OPERATION_TEMPLATE_MATCHING);
    // *** Draw ***
    private Menu menuDraw = new Menu(LABELS.OPERATION_MENU_DRAW);
    private MenuItem itemDottedandDashedLines = new MenuItem(LABELS.OPERATION_DOTTED_AND_DASHED_LINES);
    private MenuItem itemRadialGrid = new MenuItem(LABELS.OPERATION_RADIAL_GRID);
    // *** Masks ***
    private Menu menuMasks = new Menu(LABELS.OPERATION_MENU_MASKS);
    private MenuItem itemMasksToolBar = new MenuItem(LABELS.OPERATION_MASKS_TOOLBAR);
    boolean allowsPolling;

    public XmippMenuBar(Container parent, ImagePlus imp) {
        super(parent);

        this.imp = imp;

        createMenuBar();

        allowsPolling = isPollAllowed();

        enableMenuItems(imp.getStackSize(), allowsPolling);
    }

    private boolean isPollAllowed() {
        FileInfo ofi = imp.getOriginalFileInfo();

        if (ofi != null) {
            // Avoid empty filename.
            //if (!ofi.directory.trim().isEmpty() && !ofi.fileName.trim().isEmpty()) {
            File f = new File(ofi.directory + File.separator + ofi.fileName);

            return f.isFile() && f.exists();
            //}
        }

        return false;
    }

    public boolean allowsPolling() {
        return allowsPolling;
    }

    private void enableMenuItems(int nslices, boolean canPoll) {
        // Enables / disables depending on context (image or stack)
        menuReslice.setEnabled(nslices > 1);
        itemGoToSlice.setEnabled(nslices > 1);
        itemUntiltStack.setEnabled(nslices > 1);
        itemPCA.setEnabled(nslices > 1);
        itemSurfaceJ.setEnabled(nslices > 1);

        itemPoll.setEnabled(canPoll);   // Some images can't poll as they are not loaded from disk (mean, std_avg, ...).
    }

    private void createMenuBar() {
        // Transform ***
        menuTransform.add(itemToTable);
        itemToTable.addActionListener(this);

        menuTransform.add(menuFlip);
        menuFlip.add(itemFlipV);
        menuFlip.add(itemFlipH);
        itemFlipV.addActionListener(this);
        itemFlipH.addActionListener(this);

        menuTransform.add(menuRotate);
        menuRotate.add(itemRotateL);
        menuRotate.add(itemRotateR);
        itemRotateL.addActionListener(this);
        itemRotateR.addActionListener(this);

        menuTransform.add(itemCrop);
        itemCrop.addActionListener(this);

        menuTransform.add(menuReslice);
        menuReslice.add(itemResliceTop);
        menuReslice.add(itemResliceRight);
        menuReslice.add(itemStackSlicer);
        menuReslice.add(itemRadialReslice);
        itemResliceTop.addActionListener(this);
        itemResliceRight.addActionListener(this);
        itemStackSlicer.addActionListener(this);
        itemRadialReslice.addActionListener(this);

        menuTransform.add(itemStraightenCurvedObjects);
        itemStraightenCurvedObjects.addActionListener(this);

        menuTransform.add(itemUntiltStack);
        itemUntiltStack.addActionListener(this);

        menuTransform.add(menuTransformJ);
        menuTransformJ.add(itemTJAbout);
        menuTransformJ.add(itemTJAffine);
        menuTransformJ.add(itemTJCrop);
        menuTransformJ.add(itemTJEmbed);
        menuTransformJ.add(itemTJMirror);
        menuTransformJ.add(itemTJOptions);
        menuTransformJ.add(itemTJPanel);
        menuTransformJ.add(itemTJRotate);
        menuTransformJ.add(itemTJScale);
        menuTransformJ.add(itemTJShift);
        menuTransformJ.add(itemTJTurn);
        menuTransformJ.add(itemTJWebsite);

        menuTransform.add(itemSurfaceJ);

        itemTJAbout.addActionListener(this);
        itemTJAffine.addActionListener(this);
        itemTJCrop.addActionListener(this);
        itemTJEmbed.addActionListener(this);
        itemTJMirror.addActionListener(this);
        itemTJOptions.addActionListener(this);
        itemTJPanel.addActionListener(this);
        itemTJRotate.addActionListener(this);
        itemTJScale.addActionListener(this);
        itemTJShift.addActionListener(this);
        itemTJTurn.addActionListener(this);
        itemTJWebsite.addActionListener(this);
        itemSurfaceJ.addActionListener(this);

        // Denoising ***
        menuDenoising.add(itemFFT_band_pass_filter);
        menuDenoising.add(itemAnisotropicDiffusion);
        menuDenoising.add(itemMeanShiftFilter);

        itemFFT_band_pass_filter.addActionListener(this);
        itemAnisotropicDiffusion.addActionListener(this);
        itemMeanShiftFilter.addActionListener(this);

        // Open with ***
        menuOpenWith.add(item3DViewer);
        menuOpenWith.add(itemVolumeViewer3DSlicer);
        menuOpenWith.add(itemVolumeJ);
        item3DViewer.addActionListener(this);
        itemVolumeViewer3DSlicer.addActionListener(this);
        itemVolumeJ.addActionListener(this);

        // Info ***
        menuInfo.add(itemPoll);
        menuInfo.addSeparator();
        menuInfo.add(itemDuplicate);
        menuInfo.add(itemGoToSlice);
        menuInfo.add(itemSubstackMaker);
        menuInfo.addSeparator();
        menuInfo.add(itemHistogram);
        menuInfo.add(itemPlotProfile);
        menuInfo.add(itemMeasurements);
        menuInfo.add(itemAutocorrelation);
        menuInfo.add(itemLineAnalyzer);
        menuInfo.add(itemOvalProfilePlot);
        menuInfo.add(itemRadialProfilePlotAngle);
        menuInfo.add(itemRadialProfilePlotHeight);
        menuInfo.add(itemContourPlotter);
        menuInfo.add(itemFeatureJ);
        menuInfo.addSeparator();
        menuInfo.add(itemSyncWindows);
        menuInfo.add(itemSyncMeasure3D);

        // Sets initial state for poll check box
        itemPoll.setState(((iPollImageWindow) imp.getWindow()).isPoll());
        itemPoll.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                ((iPollImageWindow) imp.getWindow()).setPoll(itemPoll.getState());
            }
        });
        itemDuplicate.addActionListener(this);
        itemGoToSlice.addActionListener(this);
        itemSubstackMaker.addActionListener(this);
        itemHistogram.addActionListener(this);
        itemPlotProfile.addActionListener(this);
        itemMeasurements.addActionListener(this);
        itemAutocorrelation.addActionListener(this);
        itemLineAnalyzer.addActionListener(this);
        itemOvalProfilePlot.addActionListener(this);
        itemRadialProfilePlotAngle.addActionListener(this);
        itemRadialProfilePlotHeight.addActionListener(this);
        itemContourPlotter.addActionListener(this);
        itemFeatureJ.addActionListener(this);
        itemSyncWindows.addActionListener(this);
        itemSyncMeasure3D.addActionListener(this);

        // Thresholding ***
        menuThresholding.add(itemThreshold);
        menuThresholding.add(itemOtsuThresholding);
        menuThresholding.add(itemMultiOtsuThreshold);
        menuThresholding.add(itemMaximumEntropyThresholding);
        menuThresholding.add(itemMixtureModelingThresholding);
        menuThresholding.add(itemMultiThresholder);
        menuThresholding.add(itemRATS);
        menuThresholding.add(itemSIOX);
        menuThresholding.add(itemWatershed);

        itemThreshold.addActionListener(this);
        itemOtsuThresholding.addActionListener(this);
        itemMultiOtsuThreshold.addActionListener(this);
        itemMaximumEntropyThresholding.addActionListener(this);
        itemMixtureModelingThresholding.addActionListener(this);
        itemMultiThresholder.addActionListener(this);
        itemRATS.addActionListener(this);
        itemSIOX.addActionListener(this);

        // Binary ***
        menuBinary.add(itemVoxelCounter);
        menuBinary.add(itemErode);
        menuBinary.add(itemDilate);
        menuBinary.add(itemOpen);
        menuBinary.add(itemClose);
        menuBinary.addSeparator();
        menuBinary.add(itemFloatMorphology);
        menuBinary.addSeparator();
        menuBinary.add(itemOutline);
        menuBinary.add(itemFillHoles);
        menuBinary.add(itemSkeletonize);
        menuBinary.addSeparator();
        menuBinary.add(itemDistanceMap);
        menuBinary.add(itemUltimatePoints);
        menuBinary.add(itemWatershed);
        menuBinary.add(itemVoronoi);
        menuBinary.addSeparator();
        menuBinary.add(itemOptions);

        itemVoxelCounter.addActionListener(this);
        itemErode.addActionListener(this);
        itemDilate.addActionListener(this);
        itemOpen.addActionListener(this);
        itemClose.addActionListener(this);
        itemFloatMorphology.addActionListener(this);
        itemOutline.addActionListener(this);
        itemFillHoles.addActionListener(this);
        itemSkeletonize.addActionListener(this);
        itemDistanceMap.addActionListener(this);
        itemUltimatePoints.addActionListener(this);
        itemWatershed.addActionListener(this);
        itemVoronoi.addActionListener(this);
        itemOptions.addActionListener(this);

        // Process ***
        menuProcess.add(itemBC);
        menuProcess.add(itemEnhanceContrast);
        menuProcess.add(itemSubtractBG);
        menuProcess.add(itemGaussianBlur);
        menuProcess.add(itemConvolve);
        menuProcess.add(itemMedian);
        menuProcess.add(itemMean);
        menuProcess.add(itemFFT);
        menuProcess.add(itemCLAHE);
        menuProcess.add(itemPCA);
        menuProcess.add(itemPolarTransformer);
        menuProcess.add(itemTemplateMatching);

        itemBC.addActionListener(this);
        itemEnhanceContrast.addActionListener(this);
        itemSubtractBG.addActionListener(this);
        itemGaussianBlur.addActionListener(this);
        itemConvolve.addActionListener(this);
        itemMedian.addActionListener(this);
        itemMean.addActionListener(this);
        itemFFT.addActionListener(this);
        itemCLAHE.addActionListener(this);
        itemPCA.addActionListener(this);
        itemPolarTransformer.addActionListener(this);
        itemTemplateMatching.addActionListener(this);

        // Draw ***
        menuDraw.add(itemDottedandDashedLines);
        menuDraw.add(itemRadialGrid);
        menuDraw.add(itemConcentricCircles);

        itemDottedandDashedLines.addActionListener(this);
        itemRadialGrid.addActionListener(this);
        itemConcentricCircles.addActionListener(this);

        // Masks ***
        menuMasks.add(itemMasksToolBar);

        itemMasksToolBar.addActionListener(this);
        // Others ***
//        menuMisc.add(itemFlowJ);
//        menuMisc.add(itemFlow3J);
//        menuMisc.add(itemObjectJ);

        // Action listeners
//        itemVolumeJ.addActionListener(this);
//        itemFlowJ.addActionListener(this);
//        itemFlow3J.addActionListener(this);
//        itemObjectJ.addActionListener(this);

        // Adds menus.
        addMenu(menuTransform);
        addMenu(menuDenoising);
        addMenu(menuOpenWith);
        addMenu(menuInfo);
        addMenu(menuThresholding);
        addMenu(menuBinary);
        addMenu(menuProcess);
        addMenu(menuDraw);
        addMenu(menuMasks);
    }

    public void actionPerformed(ActionEvent e) {
        boolean captureFrame = true;    // If operation doesn't generates a new image/stack window, set this to false.
        converted = false;
        binarized = false;

        int lineWidth = Line.getWidth();    // Stores line width to restore it later.

        try {
            if (e.getSource() == itemToTable) {
                ImagesWindowFactory.openImagePlusAsGallery(imp);
            } else if (e.getSource() == itemFlipV) {
                imp.getProcessor().flipVertical();
            } else if (e.getSource() == itemFlipH) {
                imp.getProcessor().flipHorizontal();
            } else if (e.getSource() == itemRotateL) {
                IJ.run(imp, "Rotate 90 Degrees Left", "");
            } else if (e.getSource() == itemRotateR) {
                IJ.run(imp, "Rotate 90 Degrees Right", "");
            } else if (e.getSource() == itemCrop) {
                // Check if there is a line selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();
                if (roi == null || !roi.isArea()) {
                    IJ.setTool(Toolbar.RECTANGLE);
                }

                IJ.run(imp, "Crop", "");
            } else if (e.getSource() == itemHistogram) {
                IJ.run(imp, "Histogram", "");
            } else if (e.getSource() == itemPlotProfile) {
                // Check if there is a line selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();
                if (roi == null || !roi.isLine()) {
                    IJ.setTool(Toolbar.LINE);
                }

                IJ.run(imp, "Plot Profile", "");
            } else if (e.getSource() == itemMeasurements) {
                IJ.run(imp, "Set Measurements...", "");
                IJ.run(imp, "Measure", "");
            } else if (e.getSource() == itemBC) {
                IJ.run(imp, "Brightness/Contrast...", "");
            } else if (e.getSource() == itemEnhanceContrast) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Enhance Contrast", "");
            } else if (e.getSource() == itemThreshold) {
                IJ.run(imp, "Threshold...", "");
            } else if (e.getSource() == itemSubtractBG) {
                IJ.run(imp, "Subtract Background...", "");
            } else if (e.getSource() == itemGaussianBlur) {
                IJ.run(imp, "Gaussian Blur...", "");
            } else if (e.getSource() == itemConvolve) {
                IJ.run(imp, "Convolve...", "");
            } else if (e.getSource() == itemMedian) {
                IJ.run(imp, "Median...", "");
            } else if (e.getSource() == itemMean) {
                IJ.run(imp, "Mean...", "");
            } else if (e.getSource() == itemFFT) {
                IJ.run(imp, "FFT", "");
            } else if (e.getSource() == itemFFT_band_pass_filter) {
                IJ.run(imp, "Bandpass Filter...", "");
            } else if (e.getSource() == itemResliceTop) {
                IJ.run(imp, "Reslice [/]...", "slice=1.000 start=Top");
            } else if (e.getSource() == itemResliceRight) {
                IJ.run(imp, "Reslice [/]...", "slice=1.000 start=Right");
            } else if (e.getSource() == itemDuplicate) {
                ImagePlus imp2 = (new Duplicater()).duplicateStack(imp, imp.getTitle() + "_duplicated");
                imp2.show();
            } else if (e.getSource() == itemGoToSlice) {
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
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Anisotropic Diffusion...", "");
            } else if (e.getSource() == itemAutocorrelation) {
                // Check if there is a line selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();
                if (roi == null || roi.isLine()) {
                    IJ.setTool(Toolbar.RECTANGLE);
                }

                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Auto Corr", "");
            } else if (e.getSource() == itemConcentricCircles) {
                IJ.run(imp, "Concentric Circles", "");
            } else if (e.getSource() == itemLineAnalyzer) {
                // Check if there is a line selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();
                if (roi == null || !roi.isLine()) {
                    IJ.setTool(Toolbar.LINE);
                }

                IJ.run(imp, "Line Analyzer", "");
            } else if (e.getSource() == itemOvalProfilePlot) {
                // Check if there is an oval selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();

                if (roi == null || !(roi instanceof OvalRoi)) {
                    IJ.setTool(Toolbar.OVAL);
                }

                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Oval Profile", "");
            } else if (e.getSource() == itemRadialProfilePlotAngle) {
                // Check if there is an oval selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();

                if (roi == null) {
                    IJ.setTool(Toolbar.RECTANGLE);
                }

                IJ.run(imp, "Radial Profile Angle", "");
            } else if (e.getSource() == itemRadialProfilePlotHeight) {
                IJ.run(imp, "Radial Profile Height", "");
            } else if (e.getSource() == itemSyncMeasure3D) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Sync Measure 3D", "");
            } else if (e.getSource() == itemTemplateMatching) {
                IJ.run(imp, "Create Template", "");
            } else if (e.getSource() == itemCLAHE) {
                IJ.run(imp, "CLAHE ", "");
            } else if (e.getSource() == itemFloatMorphology) {
                IJ.run(imp, "Float Morphology", "");
            } else if (e.getSource() == itemMeanShiftFilter) {
                IJ.run(imp, "Mean Shift", "");
            } else if (e.getSource() == itemMaximumEntropyThresholding) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Entropy Threshold", "");
            } else if (e.getSource() == itemMixtureModelingThresholding) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Mixture Modeling", "");
            } else if (e.getSource() == itemMultiThresholder) {
                IJ.run(imp, "MultiThresholder", "");
            } else if (e.getSource() == itemMultiOtsuThreshold) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Multi OtsuThreshold", "");
            } else if (e.getSource() == itemOtsuThresholding) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Otsu Thresholding", "");
            } else if (e.getSource() == itemRATS) {
                IJ.run(imp, "RATS ", "");
            } else if (e.getSource() == itemSIOX) {
                converted = convert(imp, 24, "RGB Color");
                IJ.run(imp, "SIOX Segmentation", "");
            } else if (e.getSource() == itemContourPlotter) {
                IJ.run(imp, "ContourPlotter ", "");
            } else if (e.getSource() == itemDottedandDashedLines) {
                // Check if there is a line selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();
                if (roi == null) {
                    IJ.setTool(Toolbar.RECTANGLE);
                }

                IJ.run(imp, "Dotted Line", "");
            } else if (e.getSource() == itemFeatureJ) {
                IJ.run(imp, "FeatureJ...", "");
//            } else if (e.getSource() == itemFlowJ) {
//                IJ.run(imp, "FlowJ", "");
//            } else if (e.getSource() == itemFlow3J) {
//                IJ.run(imp, "Flow3J", "");
//            } else if (e.getSource() == itemObjectJ) {
//                IJ.run(imp, "ObjectJ", "");
            } else if (e.getSource() == itemPCA) {
                IJ.run(imp, "PCA ", "");
            } else if (e.getSource() == itemPolarTransformer) {
                IJ.run(imp, "Polar Transformer", "");
            } else if (e.getSource() == itemRadialGrid) {
                IJ.run(imp, "Radial Grid", "");
            } else if (e.getSource() == itemRadialReslice) {
                IJ.run(imp, "Radial Reslice", "");
            } else if (e.getSource() == itemStackSlicer) {
                IJ.run(imp, "StackSlicer ", "");
            } else if (e.getSource() == itemStraightenCurvedObjects) {
                // Check if there is a line selected. If not, it
                // selects the right tool for the next time.
                Roi roi = imp.getRoi();
                if (roi == null || !roi.isLine()) {
                    IJ.setTool(Toolbar.LINE);
                }

                IJ.run(imp, "Straighten...", "");
            } else if (e.getSource() == itemSubstackMaker) {
                IJ.run(imp, "Substack Maker", "");
            } else if (e.getSource() == itemSurfaceJ) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "SurfaceJ ", "");
            } else if (e.getSource() == itemSyncWindows) {
                IJ.run(imp, "Sync Windows", "");
            } else if (e.getSource() == itemTJAbout) {
                IJ.run(imp, "TJ About", "");
            } else if (e.getSource() == itemTJAffine) {
                IJ.run(imp, "TJ Affine", "");
            } else if (e.getSource() == itemTJCrop) {
                IJ.run(imp, "TJ Crop", "");
            } else if (e.getSource() == itemTJEmbed) {
                IJ.run(imp, "TJ Embed", "");
            } else if (e.getSource() == itemTJMirror) {
                IJ.run(imp, "TJ Mirror", "");
            } else if (e.getSource() == itemTJOptions) {
                IJ.run(imp, "TJ Options", "");
            } else if (e.getSource() == itemTJPanel) {
                IJ.run(imp, "TJ Panel", "");
            } else if (e.getSource() == itemTJRotate) {
                IJ.run(imp, "TJ Rotate", "");
            } else if (e.getSource() == itemTJScale) {
                IJ.run(imp, "TJ Scale", "");
            } else if (e.getSource() == itemTJShift) {
                IJ.run(imp, "TJ Shift", "");
            } else if (e.getSource() == itemTJTurn) {
                IJ.run(imp, "TJ Turn", "");
            } else if (e.getSource() == itemTJWebsite) {
                IJ.run(imp, "TJ Website", "");
            } else if (e.getSource() == itemUntiltStack) {
                IJ.run(imp, "Untilt Stack", "");
            } else if (e.getSource() == itemVolumeJ) {
                IJ.run(imp, "VolumeJ ", "");
            } else if (e.getSource() == itemVolumeViewer3DSlicer) {
                IJ.run(imp, "Volume Viewer", "");
            } else if (e.getSource() == itemVoxelCounter) {
                converted = convert(imp, 8, "8-bit");
                IJ.run(imp, "Voxel Counter", "");
            } else if (e.getSource() == item3DViewer) {
                ImagesWindowFactory.openImagePlusAs3D(imp);
                captureFrame = false;
            } else if (e.getSource() == itemErode) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Erode", "");
            } else if (e.getSource() == itemDilate) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Dilate", "");
            } else if (e.getSource() == itemOpen) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Open", "");
            } else if (e.getSource() == itemClose) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Close-", "");
            } else if (e.getSource() == itemOutline) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Outline", "");
            } else if (e.getSource() == itemFillHoles) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Fill Holes", "");
            } else if (e.getSource() == itemSkeletonize) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Skeletonize", "");
            } else if (e.getSource() == itemDistanceMap) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Distance Map", "");
            } else if (e.getSource() == itemUltimatePoints) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Ultimate Points", "");
            } else if (e.getSource() == itemWatershed) {
                binarized = makeBinary(imp);
                IJ.run(imp, "Watershed", "");
            } else if (e.getSource() == itemVoronoi) {
                binarized = true;
                IJ.run(imp, "Voronoi", "");
            } else if (e.getSource() == itemOptions) {
                IJ.run(imp, "Options...", "");
            } else if (e.getSource() == itemMasksToolBar) {
                IJ.run("Masks Tool Bar");
            }
        } catch (Exception ex) {
            //throw new RuntimeException(ex);
        }

        // Restores line width.
        Line.setWidth(lineWidth);

        // To avoid windows without our cool menu :)
        ImageWindow iw = WindowManager.getCurrentWindow();
        if (captureFrame && !(iw instanceof StackWindowOperations || iw instanceof ImageWindowOperations)) {
            ImagesWindowFactory.captureFrame(iw.getImagePlus());
        }

        // Tells user about possible internal conversions.
        if (converted) {
            IJ.showMessage("Image changed", "Image has been converted to " + imp.getBitDepth() + " bits.");
        }

        if (binarized) {
            IJ.showMessage("Image changed", "Image has been binarized.");
        }
    }

    private boolean convert(ImagePlus imp, int bytes, String strBytes) {
        if (imp.getBitDepth() != bytes) {
            IJ.run(imp, strBytes, "");
            return true;
        }
        return false;
    }

    private boolean makeBinary(ImagePlus imp) {
        if (imp.getBitDepth() != 8) {
            IJ.run(imp, "Make Binary", "");
            return true;
        }
        return false;
    }

    public void setPollStatus(boolean poll) {
        itemPoll.setState(poll);
    }
}
