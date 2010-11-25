/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import javax.swing.AbstractCellEditor;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.table.TableCellEditor;

/**
 *
 * @author Juanjo Vega
 */
public class ImageCellEditor extends AbstractCellEditor implements TableCellEditor, ActionListener {
    //extends AbstractCellEditor implements TableCellEditor {

    private JFileChooser fileChooser = new JFileChooser();
    private String currentFile;
    private JButton button;
    JFrameCOSS frameCOSS;

    public ImageCellEditor(JFrameCOSS frameCOSS) {
        super();

        this.frameCOSS = frameCOSS;

        button = new JButton();
        button.addActionListener(this);
        button.setBorderPainted(false);
    }

    public void actionPerformed(ActionEvent e) {
        //The user has clicked the cell, so bring up the dialog.
        fileChooser.setSelectedFile(new File(currentFile));

        if (fileChooser.showOpenDialog(button) == JFileChooser.APPROVE_OPTION) {
            currentFile = fileChooser.getSelectedFile().getPath();

            // Update rows height after the image is displayed.
            SwingUtilities.invokeLater(new Runnable() {

                public void run() {
                    frameCOSS.packRows();
                }
            });
        }

        fireEditingStopped(); //Make the renderer reappear.
    }

    // This method is called when editing is completed.
    // It must return the new value to be stored in the cell.
    public Object getCellEditorValue() {
        File selected = fileChooser.getSelectedFile();

        return selected != null ? selected.getPath() : null;
    }

    // This method is called when a cell value is edited by the user.
    public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
        currentFile = (String) value;

        return button;
    }
}
