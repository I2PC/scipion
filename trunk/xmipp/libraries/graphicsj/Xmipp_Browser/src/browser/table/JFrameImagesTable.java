/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameImagesTable.java
 *
 * Created on 22-abr-2010, 12:35:07
 */
package browser.table;

import browser.LABELS;
import browser.imageitems.ImageItem;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends JFrameVolumeTable {

    /** Creates new form JFrameImagesTable */
    public JFrameImagesTable() {
        super();

        initComponents2();
    }

    private void initComponents2() {
        jbMean = new javax.swing.JButton();
        jbNormalize = new javax.swing.JButton();
        jcbShowLabels = new javax.swing.JCheckBox();
        jbStd_Deviation = new javax.swing.JButton();

        toolBar.addSeparator();

        jbMean.setText(LABELS.BUTTON_MEAN);
        jbMean.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbMeanActionPerformed(evt);
            }
        });
        toolBar.add(jbMean);

        jbNormalize.setText(LABELS.BUTTON_NORMALIZE);
        jbNormalize.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbNormalizeActionPerformed(evt);
            }
        });
        toolBar.add(jbNormalize);

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

    public void addImage(ImageItem itemImage) {
        jPanelTable.addImage(itemImage);

        setTitle(LABELS.TITLE_TABLE_WINDOW_IMAGES(jPanelTable.getImagesCount()));

        /*        jsRows.setValue(jPanelTable.getRows());
        jsColumns.setValue(jPanelTable.getColumns());*/
    }

    private void jbMeanActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbMeanActionPerformed
        jPanelTable.mean();
    }//GEN-LAST:event_jbMeanActionPerformed

    private void jbNormalizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbNormalizeActionPerformed
        jPanelTable.normalize();
    }//GEN-LAST:event_jbNormalizeActionPerformed

    private void jbStd_DeviationActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbStd_DeviationActionPerformed
        jPanelTable.std_deviation();
    }//GEN-LAST:event_jbStd_DeviationActionPerformed

    private void jcbShowLabelsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbShowLabelsActionPerformed
        jPanelTable.setShowLabels(jcbShowLabels.isSelected());
}//GEN-LAST:event_jcbShowLabelsActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jbMean;
    private javax.swing.JButton jbNormalize;
    private javax.swing.JButton jbStd_Deviation;
    protected javax.swing.JCheckBox jcbShowLabels;
    // End of variables declaration//GEN-END:variables
}
