/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs;

import browser.Cache;
import browser.imageitems.TableImageItem;
import java.io.File;
import java.util.Vector;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class TableModelMicrographs extends DefaultTableModel {

    public final static int MD_LABELS[] = {
        MDLabel.MDL_ENABLED,
        MDLabel.MDL_IMAGE,
        MDLabel.MDL_PSD,
        MDLabel.MDL_CTFINPUTPARAMS,
        MDLabel.MDL_CTFMODEL,
        MDLabel.MDL_ASSOCIATED_IMAGE1,
        MDLabel.MDL_ASSOCIATED_IMAGE2,
        MDLabel.MDL_CTFMODEL2,
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
        "Enabled",
        "Micrograph",
        "PSD",
        "-CTFINPUTPARAMS-",
        "-CTFMODEL-",
        "Xmipp CTF",
        "Xmipp CTF",
        "-CTFMODEL2-",
        "CTFFind",
        "Damping",
        "1st Zero AVG.",
        "1st Zero Disagreement",
        "1st Zero Ratio",
        "Fitting Score",
        "Fitting Correlation 1st-3rd Zero",
        "PSD Correlation At 90º",
        "PSD Radial Integral",
        "PSD Variance",
        "PSD PCA Runs Test"};
    private static final int EXTRA_COLUMNS_LABELS[] = {
        MDLabel.MDL_CTF_DEFOCUSU,
        MDLabel.MDL_CTF_DEFOCUSV};
    private static final String EXTRA_COLUMNS_NAMES[] = {
        "Defocus U",
        "Defocus V"};
    public final static int CTF_COL = 4;
    public final static int DEFOCUS_U_COL = 9;
    public final static int DEFOCUS_V_COL = 10;
    // Columns to hide.
    public int columnsToHide[];
    // Columns containing null data (not saved).
    private Vector<Integer> nullColumns = new Vector<Integer>();
    private final static int invisibleColumns[] = new int[]{3, 4, 7};
    // Data type contained by columns to set renderes properly.
    public final static int filenameColumnIndex[] = new int[]{1};
    public final static int imagesColumnIndex[] = new int[]{2, 5, 6, 8};
    public final static int doubleColumnIndex[] = new int[]{9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    // Columns to save (Defocus parameters are displayed but ignored when saving).
    private final static int columns2SaveIndex[] = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    protected static Cache cache = new Cache();
    private String rootDir;
    private boolean hasCTFData = false;

    public TableModelMicrographs(MetaData md) {
        super();

        File f = new File(md.getFilename());
        rootDir = f.getParent();

        for (int i = 0; i < COLUMNS_NAMES.length; i++) {
            addColumn(COLUMNS_NAMES[i]);
        }

        parseMetaData(md);
    }

    private void parseMetaData(MetaData md) {
        try {
            // Read metadata.
            Object row[] = new Object[MD_LABELS.length];

            // Contains field enabled ?
            int startIndex = 0;
            if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
                startIndex = 1;
            }

            long objs[] = md.findObjects();

            for (long id : objs) {
                // Field enabled does not exist.
                if (startIndex > 0) {
                    row[0] = Boolean.TRUE;
                }

                for (int i = startIndex; i < MD_LABELS.length; i++) {
                    if (md.containsLabel(MD_LABELS[i])) {
                        switch (i) {
                            case 0:
                                Integer enabled = md.getValueInt(MD_LABELS[i], id);
                                row[0] = (enabled > 0 ? true : false);
                                break;
                            case 1: // IMAGE
                            case 2: // MDL_PSD
                            case 5: // ASSOCIATED_IMAGE1
                            case 6: // ASSOCIATED_IMAGE2
                            case 8: // ASSOCIATED_IMAGE3
                                String field = md.getValueString(MD_LABELS[i], id);
                                File f = new File(rootDir, field);
                                row[i] = (Object) new TableImageItem(f, cache);
                                break;
                            case 9:
                            case 10:
                            case 11:
                            case 12:
                            case 13:
                            case 14:
                            case 15:
                            case 16:
                            case 17:
                            case 18:
                            case 19:
                                row[i] = md.getValueDouble(MD_LABELS[i], id);
                                break;
                            default:
                                row[i] = md.getValueString(MD_LABELS[i], id);
                        }
                    }
                }
                // Adds a new row data to table.
                addRow(row);
            }

            addExtraColumns();

            // Stores the list of images to hide.
            columnsToHide = getColumnsToHide(md);
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("Exception: " + ex.getLocalizedMessage());
        }
    }

    private int[] getColumnsToHide(MetaData md) {
        Vector<Integer> toHide = new Vector<Integer>();

        // Adds hidden data.
        for (int i = 0; i < invisibleColumns.length; i++) {
            toHide.add(invisibleColumns[i]);
        }

        // Adds null columns.
        // Skip "enabled" field.
        for (int i = 1; i < MD_LABELS.length; i++) {
            if (!md.containsLabel(MD_LABELS[i])) {
                toHide.add(i);
                nullColumns.add(i);
            }
        }

        if (!md.containsLabel(MDLabel.MDL_CTFMODEL)) {
            nullColumns.add(getColumnCount() - 2);
            nullColumns.add(getColumnCount() - 1);

            toHide.add(getColumnCount() - 2);
            toHide.add(getColumnCount() - 1);
        }

        // Converts into array
        int columns[] = new int[toHide.size()];
        for (int i = 0; i < toHide.size(); i++) {
            columns[i] = toHide.elementAt(i).intValue();
        }

        return columns;
    }

    private void addExtraColumns() {
        addColumn(EXTRA_COLUMNS_NAMES[0]);  // DEFOCUS_U
        addColumn(EXTRA_COLUMNS_NAMES[1]);  // DEFOCUS_V

        double defocusU, defocusV;

        for (int i = 0; i < getRowCount(); i++) {
            String ctf_file = getCTFfile(i);
            //System.out.println(" CTF_File: " + ctf_file);

            if (ctf_file != null) {
                MetaData md = new MetaData(ctf_file);

                long firstID = md.findObjects()[0];

                // DEFOCUS_U
                defocusU = md.getValueDouble(EXTRA_COLUMNS_LABELS[0], firstID);

                // DEFOCUS_V
                defocusV = md.getValueDouble(EXTRA_COLUMNS_LABELS[1], firstID);

                // Sets values.
                setValueAt(defocusU, i, getColumnCount() - 2);
                setValueAt(defocusV, i, getColumnCount() - 1);

                hasCTFData = true; // There is CTF data associated :)
            }
        }
    }

    public boolean hasCtfData() {
        return hasCTFData;
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
        return column > 0 ? false : true;
    }

    public int getColumnSize() {
        return COLUMNS_NAMES.length;
    }

    public String getCTFfile(int row) {
        String file = (String) getValueAt(row, CTF_COL);
        return file != null ? rootDir + File.separator + file : null;
    }

    private boolean[] getColumns2Save() {
        boolean toSave[] = new boolean[MD_LABELS.length];

        for (int i = 0; i < MD_LABELS.length; i++) {
            if (nullColumns.contains(i)) {
                toSave[i] = false;
            } else {
                toSave[i] = true;
            }
        }
        /*
        System.out.println(" --- --- --- --- ");
        for (int i = 0; i < toSave.length; i++) {
        System.out.println(COLUMNS_NAMES[i] + ": " + toSave[i]);
        }
        System.out.println(" --- --- --- --- ");*/

        return toSave;
    }

    public boolean save(JTable table, String fileName) {
        boolean toSave[] = getColumns2Save();

        MetaData md = new MetaData();

        // Builds labels
        for (int i = 0; i < MD_LABELS.length; i++) {
            if (toSave[i]) {
//                System.out.println(i + ": " + MD_LABELS[i]);
                md.addLabel(MD_LABELS[i]);
            }
        }

        // Sets data.
        for (int i = 0; i < table.getModel().getRowCount(); i++) {
            md.addObject();

            for (int j = 0; j < columns2SaveIndex.length; j++) {
                if (toSave[j]) {
                    int index = columns2SaveIndex[j];

                    Object item = table.getModel().getValueAt(i, index);

                    if (MD_LABELS[j] == MDLabel.MDL_ENABLED) {
                        item = (Boolean) item == true ? "1" : "-1";
                    }

                    System.err.println("@TODO Check this!");
                    md.setValueString(MD_LABELS[index], item.toString(), index);
//                    System.out.println(j + ": " + MD_LABELS[index] + " / " + item);
                }
            }
//            System.out.println("---------------------------------------------");
        }

        md.write(fileName);

        return true;
    }
}
