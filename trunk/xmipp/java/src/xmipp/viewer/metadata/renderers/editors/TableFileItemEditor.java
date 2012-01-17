/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.metadata.renderers.editors;

import java.awt.Component;
import javax.swing.DefaultCellEditor;
import javax.swing.JTable;
import javax.swing.JTextField;
import xmipp.viewer.metadata.images.TableFileItem;

/**
 *
 * @author Juanjo Vega
 */
public class TableFileItemEditor extends DefaultCellEditor {//JTextField implements TableCellEditor {

    private TableFileItem item;

    public TableFileItemEditor(JTextField jtfield) {
        super(jtfield);
    }

    @Override
    public Component getTableCellEditorComponent(JTable jtable, Object o, boolean bln, int i, int i1) {
        item = (TableFileItem) o;

        return super.getTableCellEditorComponent(jtable, item.getOriginalValue(), bln, i, i1);
    }

    @Override
    public Object getCellEditorValue() {
        return item;
    }

    @Override
    public boolean stopCellEditing() {
        item.setOriginalValue(super.getCellEditorValue().toString());

        fireEditingStopped();

        return true;
    }

    @Override
    public void cancelCellEditing() {
        fireEditingCanceled();
    }
}
