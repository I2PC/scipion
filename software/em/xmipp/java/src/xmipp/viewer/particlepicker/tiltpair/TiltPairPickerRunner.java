package xmipp.viewer.particlepicker.tiltpair;


import java.util.logging.Level;

import javax.swing.SwingUtilities;

import xmipp.ij.commons.XmippApplication;
import xmipp.utils.XmippDialog;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.ParticlePickerParams;
import xmipp.viewer.particlepicker.tiltpair.gui.TiltPairPickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;

public class TiltPairPickerRunner implements Runnable {

    
    
    private final ParticlePickerParams params;


    public TiltPairPickerRunner(String[] args) {

        params = new ParticlePickerParams(args);

    }
	// 0 --> input metadata
    // 1 --> output dir
    // 2 --> mode


    

    @Override
    public void run() {
    	try {
        
            TiltPairPicker ppicker = null;
            ppicker = new TiltPairPicker(params.inputfile, params.outputdir, params);
            
            if(params.isScipion())
                XmippApplication.setIsScipion(true);

            new TiltPairPickerJFrame(ppicker);
    	} catch (Exception e) {
    		ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
    		if (!e.getMessage().isEmpty())
    			XmippDialog.showError(null, e.getMessage());
    		
    	}
        
    }

    public static void main(String[] args) {
            TiltPairPickerRunner spr = new TiltPairPickerRunner(args);
            SwingUtilities.invokeLater(spr);

    }

}
