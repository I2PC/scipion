/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;

import java.io.File;
import org.apache.commons.cli.Option;
import xmipp.utils.Params;
import xmipp.utils.XmippWindowUtil;


/**
 *
 * @author airen
 */
public class ScipionParams extends Params {

    
    public String projectid;
    public String inputid;
    public String other;
    public String python;
    public String scripts;
    public String[] objectCommands;

    public final static String PYTHON = "python";
    public final static String PROJECT = "project";
    
    
    public ScipionParams(String args[]) {
        super(args);
    }
    
    public boolean isScipion()
    {
        return cmdLine.hasOption(PROJECT);
    }

    public void defineArgs() {
        super.defineArgs();
        Option opt = new Option(PROJECT, "");
        opt.setArgs(3);
        options.addOption(opt);
        opt = new Option(PYTHON, "");
        opt.setArgs(2);
        options.addOption(opt);
        opt = new Option(OBJECT_CMDS, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
    }

    @Override
    public void processArgs(String args[]) {
        super.processArgs(args);
        
        
        if (cmdLine.hasOption(PROJECT)) {
            
            String[] cmdargs = cmdLine.getOptionValues(PROJECT);
            projectid = cmdargs[0];
            inputid = cmdargs[1];
            other = (cmdargs.length == 3)? cmdargs[2]: "";
        }
        if(cmdLine.hasOption(PYTHON))
            {
                args = cmdLine.getOptionValues(PYTHON);
            	python = args[0];
                scripts = args[1];
            }
    }
    
    public String getCTFScript()
    {
        return scripts + File.separator + "pw_recalculate_ctf.py";
    }
    
    public String getSubsetScript()
    {
        return scripts + File.separator + "pw_create_image_subset.py";
    }

    public String getObjectCmdScript() {
        return scripts + File.separator + "pw_run_obj_cmd.py";
    }

    public String getRegisterMaskScript() {
        return scripts + File.separator + "pw_create_mask.py";
    }

    public String getPlotSqliteScript() {
        return scripts + File.separator + "pw_sqlite_plot.py";
    }
    
   
}
