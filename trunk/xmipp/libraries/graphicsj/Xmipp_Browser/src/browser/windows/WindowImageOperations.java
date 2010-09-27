/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.windows;

import browser.windows.menubar.XmippMenuBar;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.ImageWindow;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import java.util.Timer;
import java.util.TimerTask;
import javax.swing.JPanel;

/**
 *
 * @author Juanjo Vega
 */
public class WindowImageOperations extends ImageWindow {

    public WindowImageOperations(ImagePlus imp) {
        this(imp, false);
    }

    public WindowImageOperations(ImagePlus imp, boolean poll) {
        super(imp);

        JPanel previousContent = new JPanel();
        previousContent.setLayout(new ImageLayout(ic));

        Component components[] = getComponents();
        for (int i = 0; i < components.length; i++) {
            previousContent.add(components[i]);
        }

        removeAll();
        setLayout(new BorderLayout());

        add(previousContent, BorderLayout.CENTER);
        //add(new PanelImageOperations(imp), BorderLayout.EAST);
        setMenuBar(new XmippMenuBar(this, imp));
//        add(new XmippJMenuBar(imp), BorderLayout.NORTH);

        setMaximumSize(getPreferredSize());

        fixSize();

        // If poll, starts timer to reload image form disk every period.
        if (poll) {
            startTimer();
        }
    }

    private void startTimer() {
        int period = 5000;  // repeat every 5 seconds
        Timer timer = new Timer();

        timer.scheduleAtFixedRate(new TimerTask() {

            public void run() {
                IJ.run(imp, "Revert", "");

                System.out.println(" >>> Reverting: " + imp.getTitle() + ": " + System.currentTimeMillis());
            }
        }, 0, period);
    }

    private void fixSize() {
        pack();

        Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
        Point loc = getLocation();
        Dimension size = getSize();
        if (loc.y + size.height > screen.height) {
            getCanvas().zoomOut(0, 0);
        }
    }
}
