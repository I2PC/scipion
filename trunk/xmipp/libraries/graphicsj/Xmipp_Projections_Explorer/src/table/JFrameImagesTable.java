/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JFrameImagesTable.java
 *
 * Created on 22-abr-2010, 12:35:07
 */
package table;

import constants.LABELS;
import ij.IJ;
import ij.ImagePlus;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import xmipp.ImageGeneric;
import xmipp.MetaData;
import xmippij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends javax.swing.JFrame {

    private JTable jtImages;
    private ImagesTableModel tableModel;
    private ScoreImagesTableRenderer tableRenderer;

    /** Creates new form JFrameImagesTable */
    public JFrameImagesTable() {
        super();

        initComponents();

        tableModel = new ImagesTableModel();    // Creates table model

        jtImages = new JTable(tableModel);  // Creates table using previous table model.

        tableRenderer = new ScoreImagesTableRenderer();
        jtImages.setDefaultRenderer(ScoreItem.class, tableRenderer);

        jtImages.setRowSelectionAllowed(false);
        jtImages.setColumnSelectionAllowed(false);
        jtImages.setCellSelectionEnabled(false);

        jtImages.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);

        jtImages.setTableHeader(null);
        jtImages.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

        jtImages.addMouseListener(new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() > 1) {
                    int row = jtImages.rowAtPoint(e.getPoint());
                    int col = jtImages.columnAtPoint(e.getPoint());

                    ScoreItem item = (ScoreItem) tableModel.getValueAt(row, col);

                    try {
                        String filename = item.fileName;
                        ImageGeneric image = new ImageGeneric(filename);
                        image.read(ImageGeneric.FIRST_IMAGE);

                        ImagePlus imp = XmippImageConverter.convertToImageJ(image);
                        imp.setTitle(filename);
                        imp.show();
                    } catch (Exception ex) {
                        throw new RuntimeException(ex);
                    }
                }
            }
        });

        // Adds components to GUI.
        JScrollPane jsPane = new JScrollPane();
        jsPane.setViewportView(jtImages);
        add(jsPane, BorderLayout.CENTER);
    }

    public void loadScoreFile(String scoreFile) {
        clear();

        buildImagesTableModel(scoreFile, tableModel);

        // Shows mean.
        ImagePlus mean = calculateMean(tableModel.getData());
        jlMean.setIcon(new ImageIcon(mean.getImage()));

        // Finally sets remaining stuff.
        setTitle("Analysis Results: " + tableModel.getColumnCount() + " images ("
                + tableModel.getGoodCount() + " good / " + tableModel.getBadCount() + " bad)");

        if (tableModel.getData().size() > 0) {
            adjustTableSize();
        }
    }

    private static ImagePlus calculateMean(ArrayList<ScoreItem> data) {
        ArrayList<String> files = new ArrayList<String>();

        for (int i = 0; i < data.size(); i++) {
            ScoreItem item = data.get(i);

            if (item.good) {    // Only "good" results are used for mean.
                files.add(item.fileName);
            }
        }

        return ImagesTableModel.mean(files);
    }

    private static void buildImagesTableModel(String filename, ImagesTableModel tablemodel) {
        try {
            // Parse md
            MetaData md = new MetaData(filename);

            long ids[] = md.findObjects();
            for (long id : ids) {
                tablemodel.addItem(new ScoreItem(md, id));
            }

            tablemodel.sortItems();
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }

    public void setWidth(int width) {
        Dimension d = getSize();

        d.width = width;

        setSize(d);
    }

    public void setHeight(int height) {
        Dimension d = getSize();

        d.height = height;

        setSize(d);
        //pack();
    }

    private void adjustTableSize() {
        // Tricky way to calculate table height
        ScoreItem item = ((ScoreItem) tableModel.getValueAt(0, 0));
        item.getImagePlus();

        int fontHeight = getFontMetrics((new JLabel(item.getLabel())).getFont()).getHeight();

        jtImages.setRowHeight(item.getHeight() + fontHeight * 2);
        for (int i = 0; i < jtImages.getColumnCount(); i++) {
            jtImages.getColumnModel().getColumn(i).setPreferredWidth(item.getWidth());
        }

        setHeight(item.getHeight() + fontHeight * 8);  // Window height.
    }

    public void clear() {
        tableModel.clear();
    }

    public int getRows() {
        return tableModel.getRowCount();
    }

    public int getColumns() {
        return tableModel.getColumnCount();
    }

    public void addImage(ScoreItem scoreItem) {
        tableModel.addItem(scoreItem);
    }

    public static void main(String args[]) {
        try {
            JFrameImagesTable frame = new JFrameImagesTable();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.loadScoreFile("/home/juanjo/Desktop/score.xmd");
            frame.pack();
            frame.setLocationRelativeTo(null);

            frame.setVisible(true);
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jlMean = new javax.swing.JLabel();

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosed(java.awt.event.WindowEvent evt) {
                formWindowClosed(evt);
            }
        });

        jlMean.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_AVERAGE_IMAGE));
        getContentPane().add(jlMean, java.awt.BorderLayout.WEST);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void formWindowClosed(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosed
        clear();
    }//GEN-LAST:event_formWindowClosed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jlMean;
    // End of variables declaration//GEN-END:variables
}
