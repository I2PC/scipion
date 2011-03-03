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
import java.io.File;
import java.util.ArrayList;
import java.util.Vector;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import xmipp.MDLabel;
import xmipp.MetaData;

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

                    IJ.open(item.fileName);
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

        // Loads data from score files.
        String items[][] = loadScoreFile("", scoreFile);

        // Creates graphical items.
        Vector<ScoreItem> pcaScoreItems = buildScoreItems(items, true);

        // Adds items to table.
        for (int i = 0; i < pcaScoreItems.size(); i++) {
            tableModel.addItem(pcaScoreItems.elementAt(i));
        }

        // Shows mean.
        ImagePlus mean = calculateMean();
        jlMean.setIcon(new ImageIcon(mean.getImage()));

        // Finally sets remaining stuff.
        setTitle("Analysis Results: " + tableModel.getColumnCount() + " images.");

        if (tableModel.getData().size() > 0) {
            adjustTableSize();
        }
    }

    private ImagePlus calculateMean() {
        Vector<ScoreItem> data = tableModel.getData();
        String files[] = new String[data.size()];

        for (int i = 0; i < data.size(); i++) {
            files[i] = tableModel.getData().elementAt(i).fileName;
        }

        return ImagesTableModel.mean(files);
    }

    private Vector<ScoreItem> buildScoreItems(String score_data[][], boolean good) {
        Vector<ScoreItem> items = new Vector<ScoreItem>();

        for (int i = 0; i < score_data.length; i++) {
            items.add(new ScoreItem(score_data[i][0],
                    Double.parseDouble(score_data[i][1]),
                    good));
        }

        return items;
    }

    private static String[][] loadScoreFile(String directory, String fileName) {
        if ((fileName == null) || (fileName.isEmpty())) {
            return null;
        }

        if (!directory.isEmpty() && !directory.endsWith(File.separator)) {
            directory += File.separator;
        }

        // Parse metadata
        MetaData md = new MetaData(directory + fileName);

        long ids[] = md.findObjects();

        ArrayList<String[]> filesList = new ArrayList<String[]>();

        for (long id : ids) {
            String file = md.getValueString(MDLabel.MDL_IMAGE, id);
            double score = md.getValueDouble(MDLabel.MDL_ZSCORE, id);

            filesList.add(new String[]{file, String.valueOf(score)});
        }

        return filesList.toArray(new String[filesList.size()][2]);
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
        item.getImage();

        int fontHeight = getFontMetrics((new JLabel(item.getLabel())).getFont()).getHeight();

        jtImages.setRowHeight(item.h + fontHeight * 2);
        for (int i = 0; i < jtImages.getColumnCount(); i++) {
            jtImages.getColumnModel().getColumn(i).setPreferredWidth(item.w);
        }

        setHeight(item.h + fontHeight * 8);  // Window height.
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

        jlMean.setBorder(javax.swing.BorderFactory.createTitledBorder(LABELS.LABEL_MEAN_IMAGE));
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
