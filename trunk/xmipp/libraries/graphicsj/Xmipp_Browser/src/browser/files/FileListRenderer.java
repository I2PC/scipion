/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.files;

import browser.ICONS_MANAGER;
import browser.REGISTERED_FILE_ITEMS;
import browser.imageitems.listitems.FileItem;
import browser.imageitems.listitems.FolderFileItem;
import java.awt.BorderLayout;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.ListCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class FileListRenderer extends JPanel implements ListCellRenderer {

    private JLabel jlItem;

    /** Creates new form FileListRenderer */
    public FileListRenderer() {
        super();

        setLayout(new BorderLayout());

        jlItem = new JLabel();
        add(jlItem, java.awt.BorderLayout.NORTH);
    }

    public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
        FileItem file = (FileItem) value;

        jlItem.setIcon(getIconForFile(file));

        jlItem.setText("<html>" + file.getLabel() + "<br><i><small>"
                + (index > 0 ? file.getDescription() : " ") + "</small></i></html>");

        if (isSelected) {
            setBackground(list.getSelectionBackground());
            setForeground(list.getSelectionForeground());
        } else {
            setBackground(list.getBackground());
            setForeground(list.getForeground());
        }

        return this;
    }

    protected static ImageIcon getIconForFile(FileItem file) {
        if (file instanceof FolderFileItem) {
            return ICONS_MANAGER.FILE_TYPE_DIRECTORY;
        } else {

            int index = REGISTERED_FILE_ITEMS.getIndexForRegisteredItem(file.getFile());
            if (index >= 0) {
                return ICONS_MANAGER.FILE_TYPES_ICONS[index];
            }
        }

        return ICONS_MANAGER.FILE_TYPE_UNKNOWN;
    }
}
