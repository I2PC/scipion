/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.windows;

import browser.LABELS;
import ij.IJ;
import ij.ImagePlus;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;
import java.awt.Choice;
import java.awt.GridBagConstraints;
import java.awt.Label;
import java.awt.event.ItemEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import javax.swing.JSpinner;
import javax.swing.JSpinner.DefaultEditor;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 *
 * @author Juanjo Vega
 */
public class PanelStackOperations extends PanelImageOperations {

    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    private Choice cReslice;
    private Choice cOpenAs;
    private JSpinner jsGoToSlice;

    public PanelStackOperations(final ImagePlus ip) {
        super(ip);

        // Reslice.
        Label labelReslice = new Label(LABELS.LABEL_RESLICE);
        GridBagConstraints gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 4;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        add(labelReslice, gridBagConstraints);

        cReslice = new Choice();
        cReslice.addItem(" -- Select op --");
        cReslice.addItem(LABELS.OPERATION_RESLICE_TOP);
        cReslice.addItem(LABELS.OPERATION_RESLICE_RIGHT);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 4;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        add(cReslice, gridBagConstraints);

        // Go to slice.
        Label labelGoTo = new Label(" -- Select op --");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 5;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        add(labelGoTo, gridBagConstraints);

        SpinnerModel sliceModel = new SpinnerNumberModel(1, 1, ip.getStackSize(), 1);
        jsGoToSlice = new JSpinner(sliceModel);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 5;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        add(jsGoToSlice, gridBagConstraints);

        // Open As
        Label labelOpenAs = new Label(LABELS.LABEL_OPEN_AS);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 6;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        add(labelOpenAs, gridBagConstraints);

        cOpenAs = new Choice();
        cOpenAs.addItem(" -- Select op --");
        cOpenAs.addItem(LABELS.OPERATION_OPEN_AS_3D);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 6;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        add(cOpenAs, gridBagConstraints);

        jsGoToSlice.setFocusable(true);
        ((DefaultEditor) jsGoToSlice.getEditor()).getTextField().addKeyListener(new KeyListener() {

            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ENTER
                        || e.getKeyCode() == KeyEvent.VK_UP
                        || e.getKeyCode() == KeyEvent.VK_DOWN) {
                    goToSelectedSlice();
                }
            }
        });

        jsGoToSlice.addChangeListener(new ChangeListener() {

            public void stateChanged(ChangeEvent e) {
                Integer value = (Integer) jsGoToSlice.getValue();

                if (value < 1) {
                    value = 1;
                } else if (value > ip.getStackSize()) {
                    value = ip.getStackSize();
                }

                jsGoToSlice.setValue(value);

                goToSelectedSlice();
            }
        });

        cReslice.addItemListener(this);
        cOpenAs.addItemListener(this);
    }

    @Override
    public void itemStateChanged(ItemEvent e) {
        try {
            ImagePlus imp;

            if (e.getSource() == cReslice) {
                switch (cReslice.getSelectedIndex()) {
                    case 1: // Reslice TOP
                        IJ.run("Reslice [/]...", "slice=1.000 start=Top");

                        imp = IJ.getImage();
                        imp.getWindow().dispose();

                        ImagesWindowFactory.openImageWindow(imp);
                        break;
                    case 2: // Reslice RIGHT
                        IJ.run("Reslice [/]...", "slice=1.000 start=Right");

                        imp = IJ.getImage();
                        imp.getWindow().dispose();

                        ImagesWindowFactory.openImageWindow(imp);
                }
            } else if (e.getSource() == cOpenAs) {
                switch (cOpenAs.getSelectedIndex()) {
                    case 1: // 3D
                        run3DViewer(IJ.getImage());
                        break;
                }
            }
        } catch (Exception ex) {
        } finally {
            cReslice.select(0);
            cOpenAs.select(0);
        }

        super.itemStateChanged(e);
    }

    private void goToSelectedSlice() {
        int slice = (Integer) jsGoToSlice.getValue();
        IJ.getImage().setSlice(slice);
    }

    private void run3DViewer(ImagePlus ip) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);

        // Adds the sphere image plus to universe.
        new StackConverter(ip).convertToRGB();
        Content c = universe.addVoltex(ip);
        c.displayAs(Content.VOLUME);

        universe.show();    // Shows...
    }
}
