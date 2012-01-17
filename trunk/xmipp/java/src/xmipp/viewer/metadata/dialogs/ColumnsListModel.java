/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.metadata.dialogs;

import javax.swing.ListModel;
import javax.swing.event.ListDataListener;

/**
 *
 * @author Juanjo Vega
 */
public class ColumnsListModel implements ListModel {

    protected String labels[];
    protected boolean states[];

    public ColumnsListModel(String[] labels) {
        this.labels = labels;

        states = new boolean[labels.length];
        for (int i = 0; i < states.length; i++) {
            states[i] = true;
        }
    }

    public ColumnsListModel(String[] labels, boolean[] states) {
        this.labels = labels;
        this.states = states;
    }

    public int getSize() {
        return labels.length;
    }

    public Object getElementAt(int index) {
        return new Object[]{labels[index], states[index]};
    }

    public void setStates(boolean[] states) {
        this.states = states;
    }

    public boolean[] getStates() {
        return states;
    }

    public void toggleItem(int index) {
        states[index] = !states[index];
    }

    public void addListDataListener(ListDataListener l) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    public void removeListDataListener(ListDataListener l) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }
}
