/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.jni;

 public class CTFParams
    {
        public double Q0;
        public double Cs;
        public double downsampleFactor;
        public double samplingRate; 
        public double Ts;
        public double kV; 
        public double defU; 
        public double defV; 
        public double D;
        public double lambda;
        public String psd;
        public String psdenhanced;

        public CTFParams(String psd, String psdenhanced, 
                double Q0, double Cs, double Ts, double kV, double defU, double defV)
        {
            this.psd = psd;
            this.psdenhanced = psdenhanced;
            this.Q0 = Q0;
            this.Cs = Cs;
            this.Ts = Ts;
            this.kV = kV;
            this.defU = defU;
            this.defV = defV;
            this.lambda = lambda(kV);
            System.out.printf("psd:%s, Q0:%.2f Cs:%.2f\n", psd, Q0, Cs);
        }
        
        
        public static double lambda(double Kv) {
            double local_kV = Kv * 1.0e3;

            return 12.3 / Math.sqrt(local_kV * (1. + local_kV * 1e-6));
        }
        
        
    }
