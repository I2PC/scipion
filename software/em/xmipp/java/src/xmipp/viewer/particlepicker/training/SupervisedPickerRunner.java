package xmipp.viewer.particlepicker.training;


import java.util.logging.Level;
import javax.swing.SwingUtilities;
import org.apache.commons.cli.ParseException;
import xmipp.ij.commons.XmippApplication;
import xmipp.utils.XmippDialog;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.ParticlePickerParams;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;

public class SupervisedPickerRunner implements Runnable {
    private ParticlePickerParams params;

    public SupervisedPickerRunner(String[] args) {
        params = new ParticlePickerParams(args);
    }

    @Override
    public void run() {
        
    	try {
            SupervisedParticlePicker ppicker = null;
            
            ppicker = new SupervisedParticlePicker(params.inputfile, params.outputdir, params.threads,
                                                   params.fast, params.incore, params);
            if(params.isScipion())
                XmippApplication.setIsScipion(true);
            new SupervisedPickerJFrame(ppicker);
    	} catch (Exception e) {
    		ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
    		if (!e.getMessage().isEmpty())
    			XmippDialog.showException(null, e);
    	}
    }

    public static void main(String[] args) {
            SupervisedPickerRunner spr = new SupervisedPickerRunner(args);
            SwingUtilities.invokeLater(spr);
    }
}
