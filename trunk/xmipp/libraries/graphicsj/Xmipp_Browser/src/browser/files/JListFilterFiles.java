package browser.files;

import javax.swing.JList;
import javax.swing.JTextField;
import javax.swing.ListModel;

public class JListFilterFiles extends JList {

    //private JTextField input;
    public JListFilterFiles() {
        super();
    }

    /**
     * Associates filtering document listener to text
     * component.
     */
    public void installJTextField(JTextField input) {
        if (input != null) {
            FilterFilesModel model = (FilterFilesModel) getModel();
            input.getDocument().addDocumentListener(model);
        }
    }

    /**
     * Disassociates filtering document listener from text
     * component.
     */
    public void uninstallJTextField(JTextField input) {
        if (input != null) {
            FilterFilesModel model = (FilterFilesModel) getModel();
            input.getDocument().removeDocumentListener(model);
        }
    }

    /**
     * Doesn't let model change to non-filtering variety
     */
    @Override
    public void setModel(ListModel model) {
        if (!(model instanceof FilterFilesModel)) {
            throw new IllegalArgumentException();
        } else {
            super.setModel(model);
        }
    }

    /**
     * Adds item to model of list
     */
    public void addElement(Object element) {
        ((FilterFilesModel) getModel()).addElement(element);
    }
}
