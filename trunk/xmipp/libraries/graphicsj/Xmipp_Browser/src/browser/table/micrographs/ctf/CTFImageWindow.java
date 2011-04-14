/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.ctf;

import browser.LABELS;
import browser.table.micrographs.iMicrographsGUI;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.ImageWindow;
import ij.gui.Roi;
import ij.process.EllipseFitter;
import ij.process.ImageStatistics;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Component;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 *
 * @author Juanjo Vega
 */
public class CTFImageWindow extends ImageWindow {

    private final static String XMIPP_CTF_ESTIMATE_FROM_PSD = "xmipp_ctf_estimate_from_psd";
    private final static String XMIPP_CTF_SORT_PSDS = "xmipp_ctf_sort_psds";
    private Button button = new Button(LABELS.LABEL_RECALCULATE_CTF);
    protected EllipseFitter ellipseFitter = new EllipseFitter();
    private EllipseCTF ellipseCTF;
    private iMicrographsGUI micrographsGUI;

    public CTFImageWindow(ImagePlus imp, String CTFFilename, iMicrographsGUI micrographsGUI) {
        super(imp);

        this.micrographsGUI = micrographsGUI;
        ellipseCTF = new EllipseCTF(CTFFilename, imp.getWidth());

        button.setEnabled(false);

        imp.getCanvas().addMouseListener(new MouseListener() {

            public void mouseClicked(MouseEvent e) {
            }

            public void mousePressed(MouseEvent e) {
            }

            public void mouseReleased(MouseEvent e) {
                // If there is no ellipse ROI, button is disabled.
                button.setEnabled(fitEllipse());
            }

            public void mouseEntered(MouseEvent e) {
            }

            public void mouseExited(MouseEvent e) {
            }
        });

        Panel previousContent = new Panel();
        previousContent.setLayout(new ImageLayout(ic));

        Component components[] = getComponents();
        for (int i = 0; i < components.length; i++) {
            previousContent.add(components[i]);
        }

        removeAll();
        setLayout(new BorderLayout());

        add(previousContent, BorderLayout.CENTER);

        Panel panel = new Panel();
        panel.add(button);
        add(panel, BorderLayout.SOUTH);

        setMaximumSize(getPreferredSize());

        button.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                recalculateCTF();
            }
        });

        pack();
        imp.updateImage();
    }

    public boolean fitEllipse() {
        boolean fitted = false;
        Roi roi = imp.getRoi();

        if (roi != null && roi.isArea()) {
            IJ.run(imp, "Fit Ellipse", "");
            ImageStatistics is = imp.getStatistics();

            ellipseFitter = new EllipseFitter();
            ellipseFitter.fit(imp.getProcessor(), is);

            // Centers ellipse in image.
            roi = imp.getRoi();
            if (roi != null) {
                Rectangle r = roi.getBounds();

                double xCenter = r.x + r.width / 2;
                double yCenter = r.y + r.height / 2;

                double dx = imp.getWidth() / 2 - xCenter;
                double dy = imp.getHeight() / 2 - yCenter;

                roi.setLocation((int) (r.x + dx), (int) (r.y + dy));
                imp.draw();

                // Store values for later use.
                fitted = true;
            }
        }

        return fitted;
    }

    private void recalculateCTF() {
        ellipseCTF.calculateDefocus(ellipseFitter.minor / 2, ellipseFitter.major / 2);

        // Invokes xmipp to recalculate image and update file.
//        estimate_from_psd(ellipseCTF);
//        ctf_sort_psds(micrographFile);

        if (micrographsGUI != null) {
            micrographsGUI.refresh();
        }
    }

    private static void estimate_from_psd(EllipseCTF ellipseCTF) {
        String command[] = new String[]{
            XMIPP_CTF_ESTIMATE_FROM_PSD,
            "--sampling_rate", String.valueOf(ellipseCTF.getSamplingRate()),
            "--kV", String.valueOf(ellipseCTF.getVoltage()),
            "--Cs", String.valueOf(ellipseCTF.getSphericalAberration()),
            "--defocusU", String.valueOf(ellipseCTF.getDefocusU()),
            "--psd"
        };

        System.out.print(" >>> ");
        for (int i = 0; i < command.length; i++) {
            System.out.print(command[i] + " ");
        }
        System.out.println();

        runCommand(command);
    }

    private static void ctf_sort_psds(String filename) {
        String command[] = new String[]{
            XMIPP_CTF_SORT_PSDS,
            "-i", filename};

        System.out.print(" >>> ");
        for (int i = 0; i < command.length; i++) {
            System.out.print(command[i] + " ");
        }
        System.out.println();

        runCommand(command);
    }

    private static void runCommand(String command[]) {
        try {
            Process p = Runtime.getRuntime().exec(command);

            BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

            // read the output from the command
            String s;
            while ((s = stdInput.readLine()) != null) {
                System.out.println(s);
            }
        } catch (IOException ex) {
            IJ.error("Error running command: " + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }
}
