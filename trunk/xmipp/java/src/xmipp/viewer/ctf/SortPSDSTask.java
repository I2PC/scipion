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

    private final static String XMIPP_CTF_SORT_PSDS = "xmipp_ctf_sort_psds";

    public SortPSDSTask(String filename) {
        super(XMIPP_CTF_SORT_PSDS + " -i " + filename);
    }

    public SortPSDSTask(String filename, iTaskCompletionListener commandsListener) {
        super(XMIPP_CTF_SORT_PSDS + " -i " + filename, commandsListener);
    }
}
