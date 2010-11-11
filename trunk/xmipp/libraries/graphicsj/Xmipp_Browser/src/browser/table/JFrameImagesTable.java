/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameImagesTable.java
 *
 * Created on Jul 28, 2010, 11:40:56 AM
 */
package browser.table;

import browser.table.normalization.JPanelNormalization;
import browser.LABELS;
import browser.imageitems.listitems.AbstractImageItem;
import java.awt.BorderLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JToggleButton;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends JFrameVolumeTable {

    private javax.swing.JButton jbMean;
    private javax.swing.JToggleButton jbNormalize;
    private javax.swing.JButton jbStd_Deviation;
    protected javax.swing.JCheckBox jcbShowLabels;
    private JPanelNormalization panelNormalization = new JPanelNormalization(this);

    /** Creates new form JFrameImagesTable */
    public JFrameImagesTable() {
        super();

        initComponents2();
    }

    private void initComponents2() {
        jbMean = new JButton();
        jbNormalize = new JToggleButton();
        jcbShowLabels = new JCheckBox();
        jbStd_Deviation = new JButton();

        toolBar.addSeparator();

        jbNormalize.setText(LABELS.BUTTON_NORMALIZE);
        jbNormalize.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbNormalizeActionPerformed(evt);
            }
        });
        toolBar.add(jbNormalize);

        jbMean.setText(LABELS.BUTTON_MEAN);
        jbMean.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbMeanActionPerformed(evt);
            }
        });
        toolBar.add(jbMean);

        jcbShowLabels.setText(LABELS.BUTTON_SHOW_LABELS);
        jcbShowLabels.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbShowLabelsActionPerformed(evt);
            }
        });
        jpControls.add(jcbShowLabels);

        jbStd_Deviation.setText(LABELS.BUTTON_STD_DEVIATION);
        jbStd_Deviation.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbStd_DeviationActionPerformed(evt);
            }
        });
        toolBar.add(jbStd_Deviation);

        pack();
    }
    /*
    public void normalize(double min, double max) {
    jPanelTable.normalize(min, max);
    }

    public void normalizeAuto() {
    jPanelTable.normalizeAuto();
    }

    public void disableNormalization() {
    jPanelTable.disableNormalization();
    }*/

    public void addImage(AbstractImageItem itemImage) {
        jPanelTable.addImage(itemImage);

        setTitle(LABELS.TITLE_TABLE_WINDOW_IMAGES(jPanelTable.getImagesCount()));
    }

    private void jbMeanActionPerformed(java.awt.event.ActionEvent evt) {
        jPanelTable.mean();
    }

    private void jbNormalizeActionPerformed(java.awt.event.ActionEvent evt) {
        showNormalizePanel(jbNormalize.isSelected());
    }

    protected void showNormalizePanel(boolean show) {
        if (show) {
            add(panelNormalization, BorderLayout.EAST);
        } else {
            remove(panelNormalization);
        }

        repaint();
    }

    private void jbStd_DeviationActionPerformed(java.awt.event.ActionEvent evt) {
        jPanelTable.std_deviation();
    }

    private void jcbShowLabelsActionPerformed(java.awt.event.ActionEvent evt) {
        jPanelTable.setShowLabels(jcbShowLabels.isSelected());
    }
}
