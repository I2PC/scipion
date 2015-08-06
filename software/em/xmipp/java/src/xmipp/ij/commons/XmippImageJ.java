package xmipp.ij.commons;

import ij.ImageJ;

import java.awt.event.WindowEvent;

public class XmippImageJ extends ImageJ
{
	
    public static final String duplicateFilter = "Duplicate";
    public static final String bandPassFilter = "Bandpass Filter...";
    public static final String anisotropicDiffFilter = "Anisotropic Diffusion...";
    public static final String substractBackgroundFilter = "Subtract Background...";
    public static final String gaussianBlurFilter = "Gaussian Blur...";
    public static final String brightnessContrastFilter = "Brightness/Contrast...";
    public static final String invertLUTFilter = "Invert LUT";
    public static final String enhanceContrastFilter = "Enhance Contrast...";

    
	public void windowClosing(WindowEvent e) {
		setVisible(false);
	}
	
	
}
