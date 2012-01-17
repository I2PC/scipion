/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.micrographs.ctf;

import xmipp.utils.LABELS;
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
import xmipp.viewer.micrographs.ctf.tasks.EstimateFromCTFTask;
import xmipp.viewer.micrographs.ctf.tasks.TasksEngine;

/**
 *
 * @author Juanjo Vega
 */
public class CTFRecalculateImageWindow extends ImageWindow {

    private Button button = new Button(LABELS.LABEL_RECALCULATE_CTF);
    protected EllipseFitter ellipseFitter = new EllipseFitter();
    private EllipseCTF ellipseCTF;
    private TasksEngine tasksEngine;
    private String PSDFilename;
    private int row;

    public CTFRecalculateImageWindow(ImagePlus imp, String CTFFilename, String PSDFilename,
            TasksEngine tasksEngine, int row) {
        super(imp);

        this.PSDFilename = PSDFilename;
        this.tasksEngine = tasksEngine;
        this.row = row;

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

        // Add "estimate..." to tasks.
        EstimateFromCTFTask estimateFromCTFTask = new EstimateFromCTFTask(
                ellipseCTF,
                90 - ellipseFitter.angle, PSDFilename, imp.getWidth(), tasksEngine, row);

        tasksEngine.add(estimateFromCTFTask);

        dispose();
    }
}
