/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wizards;

import browser.filebrowsers.model.ImageStackListModel;
import javax.swing.DefaultListCellRenderer;
import javax.swing.ListSelectionModel;

/**
 *
 * @author Juanjo Vega
 */
public abstract class JPanelXmippFilterMetadata extends JPanelXmippFilter {

    public JPanelXmippFilterMetadata(String metadata) {
        super("");

        ImageStackListModel model = new ImageStackListModel(metadata);
        jlFileFilter.setModel(model);
        jlFileFilter.setCellRenderer(new DefaultListCellRenderer());

        jlFileFilter.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    }
}
