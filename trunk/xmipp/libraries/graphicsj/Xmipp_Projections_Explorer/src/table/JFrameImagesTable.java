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
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Vector;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;

/**
 *
 * @author Juanjo Vega
 */
public class JFrameImagesTable extends javax.swing.JFrame {

    private final static String IMAGE = "image";
    private final static String SCORE = "zscore";
    private final static int INDEX_IMAGE = 0;
    private final static int INDEX_SCORE = 1;
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

    public void loadScoreFiles(String scoreFiles[]) {
        clear();

        // Loads data from score files.
        String score_pca[][] = loadScoreFile("", scoreFiles[0]);
        String score_outliers[][] = loadScoreFile("", scoreFiles[1]);

        // Creates graphical items.
        Vector<ScoreItem> pcaScoreItems = buildScoreItems(score_pca, true);
        Vector<ScoreItem> outliersScoreItems = buildScoreItems(score_outliers, false);

        // Adds GOOD items to table.
        for (int i = 0; i < pcaScoreItems.size(); i++) {
            tableModel.addItem(pcaScoreItems.elementAt(i));
        }

        // Shows mean (ONLY for GOOD items).
        ImagePlus mean = calculateMean();
        jlMean.setIcon(new ImageIcon(mean.getImage()));

        // Adds BAD items to table.
        for (int i = 0; i < outliersScoreItems.size(); i++) {
            tableModel.addItem(outliersScoreItems.elementAt(i));
        }

        // Finally sets remaining stuff.
        setTitle("Analysis Results: " + tableModel.getColumnCount() + " images ("
                + pcaScoreItems.size() + " good + "
                + outliersScoreItems.size() + " bad)");

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
            items.add(new ScoreItem(score_data[i][INDEX_IMAGE],
                    Double.parseDouble(score_data[i][INDEX_SCORE]),
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

        ArrayList<String[]> filesList = new ArrayList<String[]>();

        try {
            RandomAccessFile f = new RandomAccessFile(directory + fileName, "r");
            String line;
            boolean oldFormat;
            int indexes[] = null;

            // Skips header
            line = f.readLine();//; XMIPP_3 * column_format *
            if (line.toUpperCase().contains("XMIPP_3")) {
                oldFormat = false;

                f.readLine();// 2nd header line ";"

                // 3rd line contains parameters order.
                String metadata[] = f.readLine().split("\\s+");
                indexes = new int[metadata.length];
                for (int i = 0; i < indexes.length; i++) {
                    indexes[i] = -1;
                }

                for (int i = 1; i < metadata.length; i++) { // 0 position is ";"
                    int index = getFieldIndex(metadata[i]);
                    if (index >= 0) {
                        indexes[index] = i - 1;
                    }
                }
            } else {  // Old format.
                oldFormat = true;
                f.seek(0);
            }

            while ((line = f.readLine()) != null && !line.trim().isEmpty()) {
                String[] s = line.trim().split("\\s+");

                // Gets file depending on format version.
                String file = (oldFormat ? s[0] : s[indexes[INDEX_IMAGE]]);
                String score = (oldFormat ? s[1] : s[indexes[INDEX_SCORE]]);

                if (!file.trim().isEmpty()) {
                    // Fixes file name by adding its parent.
                    if (!file.startsWith(File.separator)) {
                        file = directory + File.separator + file;
                    }
                    filesList.add(new String[]{file, score});
                }
            }
        } catch (Exception e) {
            IJ.showStatus("");
            IJ.showMessage("Sel_Reader", "Sel_Reader : " + e);
            e.printStackTrace();
            return null;
        }

        return filesList.toArray(new String[filesList.size()][2]);
    }

    private static int getFieldIndex(String field) {
        if (field.toLowerCase().compareTo(IMAGE) == 0) {
            return INDEX_IMAGE;
        } else if (field.toLowerCase().compareTo(SCORE) == 0) {
            return INDEX_SCORE;
        }

        return -1;
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
