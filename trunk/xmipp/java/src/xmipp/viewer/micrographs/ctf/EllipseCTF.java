/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.micrographs.ctf;

import ij.IJ;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class EllipseCTF {

    private double Q0, Cs, D, Ts, kV, lambda;
    double defU, defV;
    private double defocusU, defocusV;

    public EllipseCTF(String CTFFilename, double D) {
        this.D = D;

        loadCTFFilename(CTFFilename);
    }

    private void loadCTFFilename(String CTFFilename) {
        try {
            MetaData ctfMetaData = new MetaData(CTFFilename);

            long firstID = ctfMetaData.findObjects()[0];

            Q0 = ctfMetaData.getValueDouble(MDLabel.MDL_CTF_Q0, firstID);
            Cs = ctfMetaData.getValueDouble(MDLabel.MDL_CTF_CS, firstID);
            Ts = ctfMetaData.getValueDouble(MDLabel.MDL_CTF_SAMPLING_RATE, firstID);
            kV = ctfMetaData.getValueDouble(MDLabel.MDL_CTF_VOLTAGE, firstID);

            defU = ctfMetaData.getValueDouble(MDLabel.MDL_CTF_DEFOCUSU, firstID);
            defV = ctfMetaData.getValueDouble(MDLabel.MDL_CTF_DEFOCUSV, firstID);

            lambda = lambda(kV);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }

    public double getDefocusU() {
        return defocusU;
    }

    public double getDefocusV() {
        return defocusV;
    }

    // Calculates defocus U and V according to ellipse parameters.
    public void calculateDefocus(double minor, double major) {
        defocusU = defocus(minor, Q0, lambda, Cs, D, Ts);
        defocusV = defocus(major, Q0, lambda, Cs, D, Ts);
    }

    private double lambda(double Kv) {
        double local_kV = Kv * 1.0e3;

        return 12.3 / Math.sqrt(local_kV * (1. + local_kV * 1e-6));
    }

    private double defocus(double r, double Q0, double lambda, double Cs, double D, double Ts) {
        double R = r / (D * Ts);

        double a = lambda * R * R;
        double defocus = Math.atan(-Q0) / (Math.PI * a) - 1 / a - 0.5 * (Cs * a * lambda);

        return defocus;
    }

    public double getSamplingRate() {
        return Ts;
    }

    public double getVoltage() {
        return kV;
    }

    public double getSphericalAberration() {
        return Cs;
    }
    /*
    @Override
    public String toString() {
    return "---------------------------------------------------\n"
    + "Q0: " + Q0 + "\n"
    + "Cs: " + Cs + "\n"
    + "D: " + D + "\n"
    + "Ts: " + Ts + "\n"
    + "Kv: " + kV + "\n"
    + "lambda: " + lambda + "\n"
    + "---------------------------------------------------\n"
    + "defocusU: " + defocusU + " > (file: " + defU + ")\n"
    + "defocusV: " + defocusV + " > (file: " + defV + ")\n"
    + "---------------------------------------------------\n";
    }*/
}
