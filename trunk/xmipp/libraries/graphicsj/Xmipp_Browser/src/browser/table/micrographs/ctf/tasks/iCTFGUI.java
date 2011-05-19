/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.ctf.tasks;

/**
 *
 * @author Juanjo Vega
 */
public interface iCTFGUI {

    public void setRunning(boolean running);

    public void setRowBusy(int row);

    public void setRowIdle(int row);

    public String getFilename();

    public void done();
}
