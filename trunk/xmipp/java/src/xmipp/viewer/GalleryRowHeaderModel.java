/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer;

import javax.swing.ListModel;
import javax.swing.event.ListDataListener;

/**
 *
 * @author Juanjo Vega
 */
public class GalleryRowHeaderModel implements ListModel {
    private int first_index = 0;
    private int n = 0;

    public GalleryRowHeaderModel(int n, int first_index) {        
        this.first_index = first_index;
        this.n = n;
    }

    public GalleryRowHeaderModel(int n) {
       this.n = n;
    }

    public int getSize() {
        return n;
    }

    public void setSize(int n){
    	this.n = n;
    }
    
    public Object getElementAt(int i) {
        return new Integer(first_index + i);
    }

    public void addListDataListener(ListDataListener ll) {
    }

    public void removeListDataListener(ListDataListener ll) {
    }
}
