/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package window;

import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.ImageWindow;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import javax.swing.JPanel;

/**
 *
 * @author Juanjo Vega
 */
public class ProjectionWindow extends ImageWindow {

    public ProjectionWindow(ImagePlus imp) {
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

        fixSize();
    }

    public void update(ImagePlus newImageplus, int n_images) {
        getImagePlus().setImage(newImageplus.getImage());
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
