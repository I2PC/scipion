package xmipp.viewer.particlepicker.tiltpair;


import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.SwingUtilities;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import xmipp.utils.XmippDialog;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.tiltpair.gui.TiltPairPickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.training.*;
import xmipp.viewer.particlepicker.training.gui.SupervisedParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;

public class TiltPairPickerRunner implements Runnable {

    private Options options;

    public final static String SCIPIONOPT = "scipion";
    public final static String INPUTOPT = "input";
    public final static String OUTPUTOPT = "output";
    public final static String MODEOPT = "mode";
    public final static String THREADSOPT = "threads";
    public final static String FASTOPT = "fast";
    public final static String INCOREOPT = "incore";
    public boolean isscipionsave;
    
    public String python;
    public String script;
    public String inputfile;
    public String outputdir;
    public String protid;
    public String projectid;
    public Mode mode;
    private CommandLine cmdLine;

    public TiltPairPickerRunner(String[] args) throws ParseException {

        defineArgs();
        processArgs(args);

    }
	// 0 --> input metadata
    // 1 --> output dir
    // 2 --> mode

	// On Supervised
    // 3 --> number of threads for supervised mode
    // 4 --> fast mode for supervised mode
    // 5 --> incore for supervised mode
	// On Review
    // 3 -->external auto dir
    public void defineArgs() {
        options = new Options();
        options.addOption(INPUTOPT, true, "");
        options.addOption(OUTPUTOPT, true, "");
        options.addOption(MODEOPT, true, "");

        Option cmdoption = new Option(SCIPIONOPT, "");
        cmdoption.setArgs(4);
        options.addOption(cmdoption);
    }

    public void processArgs(String args[]) throws ParseException {

        String[] cmdargs;
        BasicParser parser = new BasicParser();
        cmdLine = parser.parse(options, args);
        inputfile = cmdLine.getOptionValue(INPUTOPT);
        outputdir = cmdLine.getOptionValue(OUTPUTOPT);
        mode = Mode.getMode(cmdLine.getOptionValue(MODEOPT));
        
        
       
        if (cmdLine.hasOption(SCIPIONOPT)) {
            isscipionsave = true;
            cmdargs = cmdLine.getOptionValues(SCIPIONOPT);
            python = cmdargs[0];
            script = cmdargs[1];
            projectid = cmdargs[2];
            protid = cmdargs[3];
            
        }

    }

    @Override
    public void run() {
        try
        {
            TiltPairPicker ppicker = null;
            ppicker = new TiltPairPicker(inputfile, outputdir, Mode.Manual);
            
            ppicker.setPython(python);
            ppicker.setScipionScript(script);
            ppicker.setProjectId(projectid);
            ppicker.setProtId(protid);

            new TiltPairPickerJFrame(ppicker);
        } catch (Exception e) {
            System.out.println("Error catched on main");
            ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
            XmippDialog.showException(null, e);

        }

    }

    public static void main(String[] args) {
        try {
            TiltPairPickerRunner spr = new TiltPairPickerRunner(args);
            SwingUtilities.invokeLater(spr);
        } catch (Exception e) {
            System.out.println("Error catched on main");
            ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
            XmippDialog.showException(null, e);

        }

    }

}
