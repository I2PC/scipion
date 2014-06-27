/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.jni;

import ij.process.EllipseFitter;


/**
 *
 * @author Juanjo Vega
 */
public class EllipseCTF {

    private double Q0, Cs, Ts, kV, lambda, lowFreq, highFreq;
    int D;
    double mddefU, mddefV;
    private double defocusU, defocusV;
    private long id;
    private EllipseFitter ellipseFitter;
    private String psd;

    public void setEllipseFitter(EllipseFitter ellipseFitter) {
        this.ellipseFitter = ellipseFitter;
    }

    public EllipseFitter getEllipseFitter() {
        return ellipseFitter;
    }
   
    
    public EllipseCTF(long id, double Q0, double Cs, double Ts, double kV, double mddefU, double mddefV, int D)
    {
        this.id = id;
        this.D = D;
        this.Q0 = Q0;
        this.Cs = Cs;
        this.Ts = Ts;
        this.kV = kV;
        this.lambda = EllipseCTF.lambda(kV);
        this.mddefU = mddefU;
        this.mddefV = mddefV;
    }
    
    public long getId()
    {
        return id;
    }

    public static double lambda(double Kv) {
            double local_kV = Kv * 1.0e3;

            return 12.3 / Math.sqrt(local_kV * (1. + local_kV * 1e-6));
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
    + "defocusU: " + defocusU + " > (file: " + mddefU + ")\n"
    + "defocusV: " + defocusV + " > (file: " + mddefV + ")\n"
    + "---------------------------------------------------\n";
    }
    
    public int getD()
    {
        return D;
    }
   
    public MetaData getCTFMd()
    {
        MetaData md = new MetaData();
        md.setColumnFormat(true);
        long id = md.addObject();
        md.setValueDouble(MDLabel.MDL_CTF_Q0, Q0, id);
        md.setValueDouble(MDLabel.MDL_CTF_CS, Cs, id);
        md.setValueDouble(MDLabel.MDL_CTF_DOWNSAMPLE_PERFORMED, 1, id);
        md.setValueDouble(MDLabel.MDL_SAMPLINGRATE, Ts,id);
        md.setValueDouble(MDLabel.MDL_CTF_VOLTAGE, kV, id);
        md.setValueDouble(MDLabel.MDL_CTF_DEFOCUSU, mddefU, id);
        md.setValueDouble(MDLabel.MDL_CTF_DEFOCUSV, mddefV, id);
        md.print();
        return md;       
    }
    
}
