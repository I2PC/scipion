/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.ctf;

import xmipp.jni.CTFParams;

/**
 *
 * @author Juanjo Vega
 */
public class EllipseCTF {

    private double Q0, Cs, D, Ts, kV, lambda, lowFreq, highFreq;
    double defU, defV;
    private double defocusU, defocusV;

   
    
    public EllipseCTF(CTFParams ctfparams, double D)
    {
        this.D = D;
        Q0 = ctfparams.Q0;
        Cs = ctfparams.Cs;
        Ts = ctfparams.Ts;
        kV = ctfparams.kV;
        lambda = ctfparams.lambda;
        defU = ctfparams.defU;
        defV = ctfparams.defV;
    }

    

    public double getDefocusU() {
        return defocusU;
    }

    public double getDefocusV() {
        return defocusV;
    }
    
    /** Set the frequencies range to calculate CTF */
    public void setFreqRange(double lowFreq, double highFreq){
    	this.lowFreq = lowFreq;
    	this.highFreq = highFreq;
    }
    
    public double getLowFreq(){
    	return lowFreq;
    }
    
    public double getHighFreq(){
    	return highFreq;
    }

    public double getQ0(){
    	return Q0;
    }

    // Calculates defocus U and V according to ellipse parameters.
    public void calculateDefocus(double minor, double major) {
        defocusU = defocus(minor, Q0, lambda, Cs, D, Ts);
        defocusV = defocus(major, Q0, lambda, Cs, D, Ts);
    }

    

    private double defocus(double r, double Q0, double lambda, double Cs, double D, double Ts) {
        double R = r / (D * Ts);

        double a = lambda * R * R;
        double defocus = Math.atan(-Q0) / (Math.PI * a) - 1 / a - 0.5 * (Cs * a * lambda);

        return -defocus;
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
