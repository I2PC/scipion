/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.models;

import ij.IJ;
import javax.swing.table.DefaultTableModel;
import browser.Cache;
import metadata.images.TableFileItem;
import metadata.images.TableImageItem;
import xmipp.Filename;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class MetaDataTableModel extends DefaultTableModel {// implements TableModelListener {

    private int ENABLED_COLUMN_INDEX;
    private int MD_LABELS[];
    private String TEXT_LABELS[];
    private String filename;
    private String blocks[];
    private int selectedBlock = 0;
    private long ids[];
    protected static Cache cache = new Cache();
    private MetaData md;

    public MetaDataTableModel(String filename) {
        super();

        this.filename = filename;
        loadBlocks(filename);
    }
//
//    public void print() {
//        System.out.println(" -------------------------------- ");
//        Vector data = getDataVector();
//        for (int i = 0; i < data.size(); i++) {
//            System.out.println(i + ": " + isRowEnabled(i));
//        }
//        System.out.println(" -------------------------------- ");
//    }

    public void selectBlock(int selectedBlock) {
        this.selectedBlock = selectedBlock;
    }

    public void reload() {
        loadBlock(selectedBlock);    // ...restore data...
    }

    private void clear() {
        cache.clear();

        while (getRowCount() > 0) {
            removeRow(0);
        }

        setColumnCount(0);
    }

    public String[] getBlocks() {
        return blocks;
    }

    private void loadBlocks(String filename) {
        try {
            blocks = MetaData.getBlocksInMetaDataFile(filename);
        } catch (Exception ex) {
            System.out.println("Exception: " + ex.getMessage());
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    public void loadBlock(int index) {
        load(getFilename(), blocks[index]);
    }

    public void load(String filename, String block) {
        try {
            clear();    // Clear the whole data.

            block += !block.isEmpty() ? Filename.SEPARATOR : "";
            md = new MetaData(block + filename);

            // Contains field enabled ?
            boolean hasEnabledField = true;
            if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
                md.addLabel(MDLabel.MDL_ENABLED);
                hasEnabledField = false;
            }

            // Store ids.
            ids = md.findObjects();

            // Reads all rows.
            MD_LABELS = md.getActiveLabels();
            // Builds columns structure.
            TEXT_LABELS = new String[MD_LABELS.length];
            for (int i = 0; i < MD_LABELS.length; i++) {
                TEXT_LABELS[i] = MetaData.label2Str(MD_LABELS[i]);
                addColumn(TEXT_LABELS[i]);
            }
            Object row[] = new Object[MD_LABELS.length];

            for (long id : ids) {
                // Field enabled does not exist, so adds it (true by default).
                if (!hasEnabledField) {
                    md.setValueInt(MDLabel.MDL_ENABLED, 1, id);
                }

                for (int i = 0; i < MD_LABELS.length; i++) {
                    int label = MD_LABELS[i];

                    Class class_ = MetaData.getLabelType(label);

                    if (class_ == String.class) {
                        String value = md.getValueString(label, id);

                        if (MetaData.isImage(label)) {
                            String path = md.getValueString(label, id, true);
                            row[i] = new TableImageItem(path, value, cache);
                        } else if (MetaData.isMetadata(label)) {
                            String path = md.getValueString(label, id, true);
                            row[i] = new TableFileItem(path, value);
                        } else if (MetaData.isTextFile(label)) {
                            String path = md.getValueString(label, id, true);
                            row[i] = new TableFileItem(path, value);
                        } else {
                            row[i] = value;
                        }
                    } else if (class_ == Double.class) {
                        row[i] = md.getValueDouble(label, id);
                    } else if (class_ == Integer.class) {
                        int value = md.getValueInt(label, id);

                        // Special case.
                        if (label == MDLabel.MDL_ENABLED) {
                            row[i] = value == 1 ? Boolean.TRUE : Boolean.FALSE;
                            ENABLED_COLUMN_INDEX = i;   // Store it for future changes.
                        } else {
                            row[i] = value;
                        }
                    } else if (class_ == Boolean.class) {
                        row[i] = md.getValueBoolean(label, id);
                    } else {
                        row[i] = "Unknown type";
                    }
                }
                // Adds a new row data to table.
                addRow(row);
            }
        } catch (Exception ex) {
            System.out.println("Exception: " + ex.getMessage());
            ex.printStackTrace();
        }
    }

    public String[] getLabels() {
        return TEXT_LABELS;
    }

    public int getEnabledColumnIndex() {
        return ENABLED_COLUMN_INDEX;
    }

    @Override
    public Class getColumnClass(int column) {
        Object item = null;
        try {
            item = getValueAt(0, column);
        } catch (Exception e) {
        }

        return item != null ? item.getClass() : Object.class;
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        //return MD_LABELS[column] == MDLabel.MDL_ENABLED;
        return column == ENABLED_COLUMN_INDEX;
    }

    @Override
    public int getColumnCount() {
        return MD_LABELS != null ? MD_LABELS.length : 0;
    }

    public String getFilename() {
        return filename;//md.getFilename();
    }

    public boolean isRowEnabled(int row) {
        return md.getValueInt(MDLabel.MDL_ENABLED, ids[row]) == 1;
    }

    public void setRowEnabled(int row, boolean enabled) {
        md.setValueInt(MDLabel.MDL_ENABLED, enabled ? 1 : 0, ids[row]);
        setValueAt(enabled, row, ENABLED_COLUMN_INDEX);
    }

    public void enableAllRows(boolean enabled) {
        // Updates table.
        for (int i = 0; i < getRowCount(); i++) {
            setRowEnabled(i, enabled);
        }
    }

    public boolean save(String fileName) {
        boolean saved = true;

        try {
            md.write(fileName);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            saved = false;
        }

        return saved;
    }
}
