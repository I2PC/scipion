package xmipp.viewer.particlepicker.training;


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
import xmipp.viewer.particlepicker.training.gui.SupervisedParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;

public class SupervisedPickerRunner implements Runnable {

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
    public String projectid;
    public String inputid;
    public String inputfile;
    public String outputdir;
    public String protid;
    public String dbpath;
    public Mode mode;
    public int threads;
    public boolean fast;
    public boolean incore;
    private CommandLine cmdLine;

    public SupervisedPickerRunner(String[] args) throws ParseException {

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
        options.addOption(THREADSOPT, true, "");
        options.addOption(FASTOPT, true, "");
        options.addOption(INCOREOPT, true, "");

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
        if (cmdLine.hasOption(THREADSOPT)) 
            threads = Integer.parseInt(cmdLine.getOptionValue(THREADSOPT));
        if (cmdLine.hasOption(FASTOPT)) 
            fast = Boolean.parseBoolean(cmdLine.getOptionValue(FASTOPT));
        if (cmdLine.hasOption(INCOREOPT)) 
            incore = Boolean.parseBoolean(cmdLine.getOptionValue(INCOREOPT));
        
       
        if (cmdLine.hasOption(SCIPIONOPT)) {
            isscipionsave = true;
            cmdargs = cmdLine.getOptionValues(SCIPIONOPT);
            python = cmdargs[0];
            script = cmdargs[1];
            
            
            if(!cmdargs[2].contains(".db"))
            {
                projectid = cmdargs[2];
                inputid = cmdargs[3];
            }
            else
            { 
                dbpath = cmdargs[2];
                protid = cmdargs[3];
            }
            
        }

    }

    @Override
    public void run() {
        try
        {
            SupervisedParticlePicker ppicker = null;
            if (mode == Mode.Manual) 
                ppicker = new SupervisedParticlePicker(inputfile, outputdir, threads, fast, incore);
            else 
                ppicker = new SupervisedParticlePicker(inputfile, outputdir, mode);
            ppicker.setPython(python);
            ppicker.setScipionScript(script);
            ppicker.setProjectId(projectid);
            ppicker.setInputId(inputid);
            ppicker.setDbPath(dbpath);
            ppicker.setProtId(protid);

            new SupervisedParticlePickerJFrame(ppicker);
        } catch (Exception e) {
            System.out.println("Error catched on main");
            ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
            XmippDialog.showException(null, e);

        }

    }

    public static void main(String[] args) {
        try {
            SupervisedPickerRunner spr = new SupervisedPickerRunner(args);
            SwingUtilities.invokeLater(spr);
        } catch (Exception e) {
            System.out.println("Error catched on main");
            ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
            XmippDialog.showException(null, e);

        }

    }

}
