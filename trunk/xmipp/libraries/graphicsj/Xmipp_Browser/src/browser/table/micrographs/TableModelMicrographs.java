/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs;

import browser.Cache;
import browser.imageitems.TableImageItem;
import java.io.File;
import java.util.Vector;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class TableModelMicrographs extends DefaultTableModel implements TableModelListener {

    public final static int MD_LABELS[] = {
        MDLabel.MDL_ENABLED,
        MDLabel.MDL_IMAGE,
        MDLabel.MDL_PSD,
        MDLabel.MDL_ASSOCIATED_IMAGE1,
        MDLabel.MDL_ASSOCIATED_IMAGE2,
        MDLabel.MDL_ASSOCIATED_IMAGE3,
        MDLabel.MDL_CTF_CRITERION_DAMPING,
        MDLabel.MDL_CTF_CRITERION_FIRSTZEROAVG,
        MDLabel.MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT,
        MDLabel.MDL_CTF_CRITERION_FIRSTZERORATIO,
        MDLabel.MDL_CTF_CRITERION_FITTINGSCORE,
        MDLabel.MDL_CTF_CRITERION_FITTINGCORR13,
        MDLabel.MDL_CTF_CRITERION_PSDCORRELATION90,
        MDLabel.MDL_CTF_CRITERION_PSDRADIALINTEGRAL,
        MDLabel.MDL_CTF_CRITERION_PSDVARIANCE,
        MDLabel.MDL_CTF_CRITERION_PSDPCARUNSTEST
    };
    private static final String COLUMNS_NAMES[] = new String[]{
        "ID",
        "Enabled",
        "Micrograph",
        "PSD",
        "Xmipp CTF",
        "Xmipp CTF",
        "CTFFind",
        "Damping",
        "1st Zero AVG.",
        "1st Zero Disagreement",
        "1st Zero Ratio",
        "Fitting Score",
        "Fitting Correlation 1st-3rd Zero",
        "PSD Correlation At 90 degree",
        "PSD Radial Integral",
        "PSD Variance",
        "PSD PCA Runs Test"};
    private static final int EXTRA_COLUMNS_LABELS[] = {
        MDLabel.MDL_CTF_DEFOCUSU,
        MDLabel.MDL_CTF_DEFOCUSV};
    private static final String EXTRA_COLUMNS_NAMES[] = {
        "Defocus U",
        "Defocus V"};
    public final static int INDEX_ID = 0;
    public final static int INDEX_ENABLED = 1;
    public final static int DEFOCUS_U_COL = 7;
    public final static int DEFOCUS_V_COL = 8;
    // Data type contained by columns to set renderes properly.
    public final static int filenameColumnIndex[] = new int[]{2};
    public final static int imagesColumnIndex[] = new int[]{3, 4, 5, 6};
    public final static int doubleColumnIndex[] = new int[]{7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    protected static Cache cache = new Cache();
    private String rootDir;
    private MetaData md;

    public TableModelMicrographs(String filename) {
        super();

        md = new MetaData(filename);

        File f = new File(md.getFilename());
        rootDir = f.getParent();

        for (int i = 0; i < COLUMNS_NAMES.length; i++) {
            addColumn(COLUMNS_NAMES[i]);
        }

        buildTable(md);

        addTableModelListener(this);
    }

    private void buildTable(MetaData md) {
        try {
            // Read metadata.
            Object row[] = new Object[MD_LABELS.length + 1];

            // Contains field enabled ?
            boolean hasEnabledField = true;
            if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
                md.addLabel(MDLabel.MDL_ENABLED);
                hasEnabledField = false;
            }

            long objs[] = md.findObjects();

            for (long id : objs) {
                // Stores id.
                row[0] = id;

                // Field enabled does not exist, so adds it (true by default).
                if (!hasEnabledField) {
                    md.setValueInt(MDLabel.MDL_ENABLED, 1, id);
                }

                for (int i = 0, col = 1; i < MD_LABELS.length; i++, col++) {
                    int label = MD_LABELS[i];

                    if (md.containsLabel(label)) {
                        switch (label) {
                            case MDLabel.MDL_ENABLED:
                                Integer enabled = md.getValueInt(label, id);
                                row[col] = enabled > 0;
                                break;
                            case MDLabel.MDL_IMAGE:
                            case MDLabel.MDL_PSD:
                            case MDLabel.MDL_ASSOCIATED_IMAGE1:
                            case MDLabel.MDL_ASSOCIATED_IMAGE2:
                            case MDLabel.MDL_ASSOCIATED_IMAGE3:
                                String field = md.getValueString(label, id);
                                File f = new File(rootDir, field);
                                row[col] = (Object) new TableImageItem(f, cache);
                                break;
                            case MDLabel.MDL_CTF_CRITERION_DAMPING:
                            case MDLabel.MDL_CTF_CRITERION_FIRSTZEROAVG:
                            case MDLabel.MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT:
                            case MDLabel.MDL_CTF_CRITERION_FIRSTZERORATIO:
                            case MDLabel.MDL_CTF_CRITERION_FITTINGSCORE:
                            case MDLabel.MDL_CTF_CRITERION_FITTINGCORR13:
                            case MDLabel.MDL_CTF_CRITERION_PSDCORRELATION90:
                            case MDLabel.MDL_CTF_CRITERION_PSDRADIALINTEGRAL:
                            case MDLabel.MDL_CTF_CRITERION_PSDVARIANCE:
                            case MDLabel.MDL_CTF_CRITERION_PSDPCARUNSTEST:
                                row[col] = md.getValueDouble(label, id);
                                break;
                            default:
                                row[col] = md.getValueString(label, id);
                        }
                    }
                }
                // Adds a new row data to table.
                addRow(row);
            }

            addExtraColumns();
        } catch (Exception ex) {
            System.out.println("Exception: " + ex.getMessage());
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    private void addExtraColumns() {
        addColumn(EXTRA_COLUMNS_NAMES[0]);  // DEFOCUS_U
        addColumn(EXTRA_COLUMNS_NAMES[1]);  // DEFOCUS_V

        double defocusU, defocusV;

        if (hasCtfData()) {
            for (int i = 0; i < getRowCount(); i++) {
                String ctf_file = getCTFfile(i);

                MetaData ctfMetaData = new MetaData(ctf_file);

                long firstID = ctfMetaData.findObjects()[0];

                // DEFOCUS_U
                defocusU = ctfMetaData.getValueDouble(EXTRA_COLUMNS_LABELS[0], firstID);

                // DEFOCUS_V
                defocusV = ctfMetaData.getValueDouble(EXTRA_COLUMNS_LABELS[1], firstID);

                // Sets values.
                setValueAt(defocusU, i, getColumnCount() - 2);
                setValueAt(defocusV, i, getColumnCount() - 1);
            }
        }
    }

    public Vector<Integer> getColumnsToHide() {
        Vector<Integer> toHide = new Vector<Integer>();

        // Hide ID column.
        toHide.add(INDEX_ID);

        // Adds null columns.
        for (int i = 0; i < MD_LABELS.length; i++) {
            if (!md.containsLabel(MD_LABELS[i])) {
                toHide.add(i + 1);
                System.out.println("To hide: " + COLUMNS_NAMES[i]);
            }
        }

        // Extra columns: DEFOCUS
        if (!hasCtfData()) {
            toHide.add(MD_LABELS.length + 1);
            toHide.add(MD_LABELS.length + 2);
        }

        return toHide;
    }

    public boolean hasCtfData() {
        return md.containsLabel(MDLabel.MDL_CTFMODEL);
    }

    /*
     * JTable uses this method to determine the default renderer/
     * editor for each cell.  If we didn't implement this method,
     * then the first column would contain text ("true"/"false"),
     * rather than a check box.
     */
    @Override
    public Class getColumnClass(int column) {
        Object item = getValueAt(0, column);
        return item != null ? item.getClass() : Object.class;
    }

    @Override
    public Object getValueAt(int row, int column) {
        return super.getValueAt(row, column);
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        return column == 0 ? false : true;
    }

    public int getColumnSize() {
        return COLUMNS_NAMES.length;
    }

    public String getCTFfile(int row) {
        long id = (Long) getValueAt(row, 0);
        String file = md.getValueString(MDLabel.MDL_CTFMODEL, id);

        if (file != null) {
            File f = new File(file);

            if (!f.isAbsolute()) {
                file = rootDir + File.separator + file;
            }
        }

        return file;
    }

    // This method is invoked TWICE for each event (I don't know how fix it for now)
    public void tableChanged(TableModelEvent e) {
        int row = e.getFirstRow();

        long id = (Long) getValueAt(row, 0);
        boolean enabled = (Boolean) getValueAt(row, 1);

        md.setValueInt(MDLabel.MDL_ENABLED, enabled ? 1 : 0, id);
    }

    public void save(String fileName) {
        md.write(fileName);
    }
}
