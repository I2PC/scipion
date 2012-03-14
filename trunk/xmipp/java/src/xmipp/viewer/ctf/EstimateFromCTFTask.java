/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.ctf;


/**
 *
 * @author juanjo
 */
public class EstimateFromCTFTask extends CommandTask {

    private final static String XMIPP_CTF_ESTIMATE_FROM_PSD = "xmipp_ctf_estimate_from_psd";

    public EstimateFromCTFTask(EllipseCTF ellipseCTF, double angle, String PSDFilename, int modelSize,
            iTaskCompletionListener commandsListener, int row) {
        super(XMIPP_CTF_ESTIMATE_FROM_PSD
                + " --sampling_rate " + ellipseCTF.getSamplingRate()
                + " --kV " + ellipseCTF.getVoltage()
                + " --Cs " + ellipseCTF.getSphericalAberration()
                + " --defocusU " + ellipseCTF.getDefocusU()
                + " --defocusV " + ellipseCTF.getDefocusV()
                + " --azimuthal_angle " + angle
                + " --psd " + PSDFilename
                + " --ctfmodelSize " + modelSize,
                row,
                commandsListener);
    }
}
