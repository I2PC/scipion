package browser.filebrowsers;

import browser.filebrowsers.model.ListModelFilesBrowser;
import javax.swing.JList;
import javax.swing.JTextField;

public class JListFileFilter extends JList {

    //private JTextField input;
    public JListFileFilter() {
        super();
    }

    /**
     * Associates filtering document listener to text
     * component.
     */
    public void installJTextField(JTextField input) {
        if (input != null) {
            ListModelFilesBrowser model = (ListModelFilesBrowser) getModel();
            input.getDocument().addDocumentListener(model);
        }
    }

    /**
     * Disassociates filtering document listener from text
     * component.
     */
    public void uninstallJTextField(JTextField input) {
        if (input != null) {
            ListModelFilesBrowser model = (ListModelFilesBrowser) getModel();
            input.getDocument().removeDocumentListener(model);
        }
    }

    /**
     * Doesn't let model change to non-filtering variety
     */
//    @Override
//    public void setModel(ListModel model) {
//        if (!(model instanceof iListModelBrowser)) {
//            throw new IllegalArgumentException();
//        } else {
//            super.setModel(model);
//        }
//    }
    /**
     * Adds item to model of list
     */
    public void addElement(Object element) {
        ((ListModelFilesBrowser) getModel()).addElement(element);
    }
}
