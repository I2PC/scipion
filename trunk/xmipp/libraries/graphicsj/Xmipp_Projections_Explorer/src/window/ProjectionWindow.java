/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package window;

import constants.LABELS;
import explorer.ProjectionsExplorer;
import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.ImageWindow;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JPanel;

/**
 *
 * @author Juanjo Vega
 */
public class ProjectionWindow extends ImageWindow {

    private PanelProjectionAnalyze panelProjectionAnalyze;
    private boolean analyze;

    public ProjectionWindow(ImagePlus imp, int n_images, ProjectionsExplorer projectionProcessor, boolean analyze) {
        super(imp);

        this.analyze = analyze;

        if (analyze) {  // If analyze, then adds button and images counter.
            JPanel previousContent = new JPanel();
            previousContent.setLayout(new ImageLayout(ic));

            Component components[] = getComponents();
            for (int i = 0; i < components.length; i++) {
                previousContent.add(components[i]);
            }

            removeAll();
            setLayout(new BorderLayout());

            add(previousContent, BorderLayout.CENTER);

            panelProjectionAnalyze = new PanelProjectionAnalyze(projectionProcessor, n_images);
            add(panelProjectionAnalyze, BorderLayout.SOUTH);
        }

        fixSize();
    }

    public void update(ImagePlus newImageplus, int n_images) {
        ImagePlus ip = getImagePlus();
        ip.getProcessor().setPixels(newImageplus.getProcessor().getPixels());
        ip.updateAndDraw();

        if (analyze) {
            panelProjectionAnalyze.setLabelImages(n_images);

            // Enables button just when there is any projection.
            panelProjectionAnalyze.bAnalyze.setEnabled(n_images > 0);
        }
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

    class PanelProjectionAnalyze extends JPanel implements ActionListener {

        public Label ln_images;
        protected Button bAnalyze;
        private ProjectionsExplorer projectionsExplorer;

        /** Creates new form PanelProjectionAnalyze */
        public PanelProjectionAnalyze(ProjectionsExplorer projectionsExplorer, int n_images) {
            this.projectionsExplorer = projectionsExplorer;

            initComponents();

            setLabelImages(n_images);
            bAnalyze.addActionListener(this);
        }

        public void setLabelImages(int n_images) {
            ln_images.setText(LABELS.LABEL_N_IMAGES(n_images));
        }

        public void actionPerformed(ActionEvent e) {
            if (e.getSource() == bAnalyze) {
                projectionsExplorer.analyzeProjection();
            }
        }

        private void initComponents() {
            setLayout(new GridLayout(2, 1));

            ln_images = new Label();
            bAnalyze = new Button(LABELS.BUTTON_ANALYZE);

            setLabelImages(0);

            add(ln_images);
            add(bAnalyze);
        }
    }
}
