/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.ctf;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.Toolbar;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextField;

/**
 *
 * @author Juanjo Vega
 */
public class CTFImageWindowValues extends CTFImageWindow {

    private Label lMajor = new Label("Major: ");
    private Label lMinor = new Label("Minor: ");
    private Label lAngle = new Label("Angle: ");
    private Label lTheta = new Label("Theta: ");
    private Label lxCenter = new Label("xCenter: ");
    private Label lyCenter = new Label("yCenter: ");
    private TextField tfMajor = new TextField(10);
    private TextField tfMinor = new TextField(10);
    private TextField tfAngle = new TextField(10);
    private TextField tfTheta = new TextField(10);
    private TextField tfxCenter = new TextField(10);
    private TextField tfyCenter = new TextField(10);

    public CTFImageWindowValues(ImagePlus imp, String CTFFilename) {
        super(imp, CTFFilename, null);

        GridLayout layout = new GridLayout(7, 2);
        Panel panel = new Panel(layout);

        panel.add(lMajor);
        panel.add(tfMajor);
        panel.add(lMinor);
        panel.add(tfMinor);
        panel.add(lAngle);
        panel.add(tfAngle);
        panel.add(lTheta);
        panel.add(tfTheta);
        panel.add(lxCenter);
        panel.add(tfxCenter);
        panel.add(lyCenter);
        panel.add(tfyCenter);

        // Save previous content...
        Panel previousContent = new Panel();
        previousContent.setLayout(new ImageLayout(ic));

        Component components[] = getComponents();
        for (int i = 0; i < components.length; i++) {
            previousContent.add(components[i]);
        }

        // ...to reset panel...
        removeAll();
        setLayout(new BorderLayout());
        add(previousContent, BorderLayout.CENTER);

        // ..and to add the new stuff.
        add(panel, BorderLayout.EAST);

        pack();
    }

    @Override
    public boolean fitEllipse() {
        boolean fitted = super.fitEllipse();

        tfMajor.setText(String.valueOf(ellipseFitter.major / 2));
        tfMinor.setText(String.valueOf(ellipseFitter.minor / 2));
        tfAngle.setText(String.valueOf(ellipseFitter.angle));
        tfTheta.setText(String.valueOf(ellipseFitter.theta));
        tfxCenter.setText(String.valueOf(imp.getWidth() / 2));
        tfyCenter.setText(String.valueOf(imp.getHeight() / 2));

        return fitted;
    }

    public static void main(String args[]) {
        new ImageJ();
        System.setProperty("plugins.dir", "/home/juanjo/Desktop/ImageJ/plugins");

        int MICROGRAPH_INDEX = 4;

        String CTFFilenames[] = {
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X1/down1_DnaB-DnaC_50000X1_Periodogramavg.ctfparam",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X7/down1_DnaB-DnaC_50000X7_Periodogramavg.ctfparam",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X9/down1_DnaB-DnaC_50000X9_Periodogramavg.ctfparam",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X50/down1_DnaB-DnaC_50000X50_Periodogramavg.ctfparam",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X71/down1_DnaB-DnaC_50000X71_Periodogramavg.ctfparam"
        };
        String ImageFilenames[] = {
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X1/down1_DnaB-DnaC_50000X1_Periodogramavg_ctfmodel_halfplane.xmp",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X7/down1_DnaB-DnaC_50000X7_Periodogramavg_ctfmodel_halfplane.xmp",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X9/down1_DnaB-DnaC_50000X9_Periodogramavg_ctfmodel_halfplane.xmp",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X50/down1_DnaB-DnaC_50000X50_Periodogramavg_ctfmodel_halfplane.xmp",
            "/home/juanjo/Desktop/20_05_2010/DnaB-DnaC_50000X71/down1_DnaB-DnaC_50000X71_Periodogramavg_ctfmodel_halfplane.xmp"
        };

        IJ.setTool(Toolbar.FREEROI);

        ImagePlus ip = IJ.openImage(ImageFilenames[MICROGRAPH_INDEX]);
        new CTFImageWindowValues(ip, CTFFilenames[MICROGRAPH_INDEX]);
    }
}
