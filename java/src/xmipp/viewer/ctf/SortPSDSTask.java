/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.ctf;

/**
 *
 * @author Juanjo Vega
 */
public class SortPSDSTask extends CommandTask {

    public SortPSDSTask(String filename) {
        super(getCommand(filename));
    }

    public SortPSDSTask(String filename, iTaskCompletionListener commandsListener) {
        super(getCommand(filename), commandsListener);
    }
    
    public static String getCommand(String filename){
    	return String.format("xmipp_ctf_sort_psds -i %s -o %s", filename, filename);
    }
}
