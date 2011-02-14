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
import browser.imageitems.listitems.XmippImageItem;
import java.awt.BorderLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JToggleButton;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends JFrameVolumeTable {

    private javax.swing.JButton jbAverage;
    private javax.swing.JToggleButton jbNormalize;
    private javax.swing.JButton jbStd_Deviation;
    protected javax.swing.JCheckBox jcbShowLabels;
    private JPanelNormalization panelNormalization = new JPanelNormalization(this);

    /** Creates new form JFrameImagesTable */
    public JFrameImagesTable(int initialRows, int initialColumns) {
        super(initialRows, initialColumns);

        initComponents2();
    }

    private void initComponents2() {
        jbAverage = new JButton();
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

        jbAverage.setText(LABELS.BUTTON_AVERAGE);
        jbAverage.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAverageActionPerformed(evt);
            }
        });
        toolBar.add(jbAverage);

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

    @Override
    public void addImageItem(XmippImageItem itemImage, int n) {
        jPanelTable.addImageItem(itemImage, n);

        setTitle(LABELS.TITLE_TABLE_WINDOW(jPanelTable.getItemsCount()));
    }

    private void jbAverageActionPerformed(java.awt.event.ActionEvent evt) {
        jPanelTable.average();
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
