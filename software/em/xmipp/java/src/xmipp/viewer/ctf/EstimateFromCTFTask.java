/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.ctf;

import xmipp.jni.EllipseCTF;


/**
 *
 * @author juanjo
 */
public class EstimateFromCTFTask extends CommandTask {

    private final static String XMIPP_CTF_ESTIMATE_FROM_PSD = "xmipp_ctf_estimate_from_psd";

    public EstimateFromCTFTask(EllipseCTF ellipseCTF, double angle, String PSDFilename, int modelSize,
            iTaskCompletionListener commandsListener, int row, String sortFn) {
        super(XMIPP_CTF_ESTIMATE_FROM_PSD
                + " --sampling_rate " + ellipseCTF.getSamplingRate()
                + " --downSamplingPerformed " + ellipseCTF.getDownsamplingFactor()
                + " --kV " + ellipseCTF.getVoltage()
                + " --Cs " + ellipseCTF.getSphericalAberration()
                + " --defocusU " + ellipseCTF.getDefocusU()
                + " --defocusV " + ellipseCTF.getDefocusV()
                + " --azimuthal_angle " + angle
                + " --psd " + PSDFilename
                + " --ctfmodelSize " + modelSize
                + " --min_freq " + ellipseCTF.getLowFreq()
                + " --max_freq " + ellipseCTF.getHighFreq()
                + " --Q0 " + ellipseCTF.getQ0()
                + " --defocus_range 5000"
                
                ,//+ " ; " + getSortCmd(sortFn),
                row,
                commandsListener);
        
        addCommand(getSortCmd(sortFn));
    }
    
    public static String getSortCmd(String filename){
    	return String.format("xmipp_ctf_sort_psds -i %s -o %s", filename, filename);
    }
}
