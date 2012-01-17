/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.metadata.renderers;

import java.awt.Component;
import javax.swing.JCheckBox;
import javax.swing.JList;
import javax.swing.ListCellRenderer;
import javax.swing.UIManager;

/**
 *
 * @author Juanjo Vega
 */
public class CheckBoxRenderer extends JCheckBox implements ListCellRenderer {

    public CheckBoxRenderer() {
        setBackground(UIManager.getColor("List.textBackground"));
        setForeground(UIManager.getColor("List.textForeground"));
    }

    public Component getListCellRendererComponent(JList listBox, Object obj, int currentindex,
            boolean isChecked, boolean hasFocus) {
        Object item[] = (Object[]) obj;

        setText((String) item[0]);
        setSelected((Boolean) item[1]);

        return this;
    }
}
