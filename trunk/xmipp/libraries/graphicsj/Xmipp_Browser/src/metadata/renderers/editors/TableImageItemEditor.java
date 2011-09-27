/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package metadata.renderers.editors;

import browser.imageitems.tableitems.GalleryImageItem;
import java.awt.Component;
import javax.swing.DefaultCellEditor;
import javax.swing.JTable;
import javax.swing.JTextField;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItemEditor extends DefaultCellEditor {

    private GalleryImageItem item;
    private String previousValue;

    public TableImageItemEditor() {
        super(new JTextField());
    }

    @Override
    public Component getTableCellEditorComponent(JTable jtable, Object o, boolean bln, int i, int i1) {
        item = (GalleryImageItem) o;

        previousValue = item.getOriginalValue();

        return super.getTableCellEditorComponent(jtable, previousValue, bln, i, i1);
    }

    @Override
    public Object getCellEditorValue() {
        String text = ((JTextField) getComponent()).getText();

        if (text.compareTo(previousValue) != 0) {
            System.out.println(" *** Saving values");

            item.setPath(text);
            item.setOriginalValue(text);
        }

        return item;
    }

    @Override
    public boolean stopCellEditing() {
        // @TODO Weird effect if still in edition mode when activates rendering.
        fireEditingStopped();

        return true;
    }

    @Override
    public void cancelCellEditing() {
        System.out.println(" >>> cancelled");

        super.cancelCellEditing();
    }
}
