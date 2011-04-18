/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs;

import browser.table.micrographs.ctf.EllipseCTF;

/**
 *
 * @author Juanjo Vega
 */
public interface iMicrographsGUI {

    public void recalculateCTF(EllipseCTF ellipseCTF, double angle, String PSDFilename);
}
