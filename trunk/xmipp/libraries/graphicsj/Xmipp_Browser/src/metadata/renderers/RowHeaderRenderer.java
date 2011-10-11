package metadata.renderers;

import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.ListCellRenderer;
import javax.swing.UIManager;
import javax.swing.border.Border;

public class RowHeaderRenderer extends JLabel implements ListCellRenderer {

    protected Border border;

    public RowHeaderRenderer() {
        super();

        setOpaque(true);

        border = BorderFactory.createCompoundBorder(
                UIManager.getBorder("TableHeader.cellBorder"),
                BorderFactory.createEmptyBorder(0, 0, 0, 2));

        // Border to show the entire label and with the same look and feel as columns.
        setBorder(border);
    }

    public Component getListCellRendererComponent(JList list, Object value, int index,
            boolean selected, boolean hasFocus) {
        setText(String.valueOf((Integer) value + 1));

        if (selected) {
            setBackground(list.getSelectionBackground());
            setForeground(list.getSelectionForeground());
        } else {
            setBackground(list.getBackground());
            setForeground(list.getForeground());
        }

        return this;
    }
}
