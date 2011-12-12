/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.models;

import javax.swing.table.DefaultTableModel;

/**
 *
 * @author Juanjo Vega
 */
public abstract class XmippTableModelRowDisabler extends DefaultTableModel {

   public abstract boolean isRowEnabled(int row);
}
