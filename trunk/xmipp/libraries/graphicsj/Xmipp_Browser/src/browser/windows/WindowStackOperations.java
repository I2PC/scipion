package browser.windows;

import browser.windows.menubar.XmippMenuBar;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageLayout;
import ij.gui.StackWindow;
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
public class WindowStackOperations extends StackWindow {

    public WindowStackOperations(ImagePlus imp, ImageCanvas ic) {
        super(imp, ic);

        JPanel previousContent = new JPanel();
        previousContent.setLayout(new ImageLayout(ic));

        Component components[] = getComponents();
        for (int i = 0; i < components.length; i++) {
            previousContent.add(components[i]);
        }

        removeAll();
        setLayout(new BorderLayout());

        add(previousContent, BorderLayout.CENTER);
        //add(new PanelStackOperations(imp), BorderLayout.EAST);
        setMenuBar(new XmippMenuBar(this, imp));
//        add(new XmippJMenuBar(imp), BorderLayout.NORTH);

        fixSize();
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
