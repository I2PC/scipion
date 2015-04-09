/***************************************************************************
 * Authors:     Juanjo Vega
 * 				J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
package xmipp.viewer.models;

import javax.swing.ListModel;
import javax.swing.event.ListDataListener;

public class GalleryRowHeaderModel implements ListModel {
    private int first_index = 0;
    private int n = 0;
    private GalleryData data = null;

    public GalleryRowHeaderModel(int n, int first_index) {      
        
        this.first_index = first_index;
        this.n = n;
    }

    public GalleryRowHeaderModel(int n) {
       this.n = n;
    }
    
    public GalleryRowHeaderModel(GalleryData data){
    	this.data = data;
    }

    public int getSize() {
    	if (data != null)
    		return data.labels.size();
        return n;
    }

    public void setSize(int n){
    	this.n = n;
    }
    
    public Object getElementAt(int i) {
    	if (data != null)
    		return data.labels.get(i).labelName;
        return new Integer(first_index + i);
    }

    public void addListDataListener(ListDataListener ll) {
    }

    public void removeListDataListener(ListDataListener ll) {
    }
}
