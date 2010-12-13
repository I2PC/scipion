/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss;

import browser.Cache;
import java.io.File;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;
import xmipp.FileName;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class TableModelCOSS extends DefaultTableModel {

    public final static MDLabel MD_LABELS[] = {
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
    private static final String columnsNames[] = new String[]{
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
    public final static int CTF_COL = 3;
    public final static int invisibleColumns[] = new int[]{3, 4, 7};
    public final static int imagesColumnIndex[] = new int[]{2, 3, 4, 5};
    protected static Cache cache = new Cache();

    public TableModelCOSS(String fileName) {
        super();

        for (int i = 0; i < columnsNames.length; i++) {
            addColumn(columnsNames[i]);
        }

        parseMetaData(fileName);
    }

    private void parseMetaData(String fileName) {
        try {
            FileName fn = new FileName(fileName);
            MetaData md = new MetaData(fn);

            String rootDir = new File(fileName).getParent();

            // Read metadata.
            Object row[] = new Object[TableModelCOSS.MD_LABELS.length];
            String field[] = new String[1];
            md.iteratorBegin();

            // Contains field enabled ?
            int startIndex = 0;
            if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
                startIndex = 1;
            }

            while (!md.iteratorEnd()) { // Retrives all values from the current file row.
                // Field enabled does not exist.
                if (startIndex > 0) {
                    row[0] = Boolean.TRUE;
                }

                for (int i = startIndex; i < TableModelCOSS.MD_LABELS.length; i++) {
                    md.getStrFromValue(TableModelCOSS.MD_LABELS[i], field);

                    switch (i) {
                        case 0:
                            Integer enabled = Integer.parseInt(field[0]);
                            row[0] = new Boolean(enabled > 0 ? true : false);
                            break;
                        case 1: // IMAGE
                        case 2: // MDL_PSD
                        case 5: // ASSOCIATED_IMAGE1
                        case 6: // ASSOCIATED_IMAGE2
                        case 8: // ASSOCIATED_IMAGE3
                            row[i] = (Object) new TableImageItemCOSS(
                                    rootDir, field[0], cache, this);
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
                            row[i] = new Double(field[0]);
                            break;
                        default:
                            row[i] = field[0];
                    }
                }

                // Adds a new row data to table.
                addRow(row);

                md.iteratorNext();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("Exception: " + ex.getLocalizedMessage());
        }
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
        return columnsNames.length;
    }

    public String getCTFfile(JTable table, int row) {
        return (String) table.getModel().getValueAt(row, CTF_COL);
    }

    public static boolean save(JTable table, String fileName) {
        MetaData md = new MetaData();

        // Builds labels
        md.addLabel(MDLabel.MDL_ENABLED);
        for (int i = 0; i < table.getModel().getColumnCount() - 1; i++) {
            md.addLabel(MD_LABELS[i]);
        }

        // Sets data.
        for (int i = 0; i < table.getModel().getRowCount(); i++) {
            md.addObject();
            //md.setValueFromStr(MDLabel.MDL_ENABLED, table.getModel().getValueAt(i, 0).toString());
            for (int j = 0; j < table.getModel().getColumnCount(); j++) {
//                System.out.println((j + 1) + ": " + MD_LABELS[j] + " / " + table.getModel().getValueAt(i, j + 1));
                Object item = table.getModel().getValueAt(i, j);

                if (MD_LABELS[j] == MDLabel.MDL_ENABLED) {
                    item = (Boolean) item == true ? "1" : "-1";
                    //System.out.println(table.getModel().getValueAt(i, j).toString());
                }

                md.setValueFromStr(MD_LABELS[j], item.toString());
            }
//            System.out.println("---------------------------------------------");
        }

        md.write(new FileName(fileName));

        return true;
    }
}
