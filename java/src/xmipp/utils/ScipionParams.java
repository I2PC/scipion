/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;

import java.io.File;
import org.apache.commons.cli.Option;


/**
 *
 * @author airen
 */
public class ScipionParams extends Params {

    
    public int port;
    public String inputid;
    public String other;
    public String scripts;
    public String[] objectCommands;

    public final static String SCIPION = "scipion";
    
    
    public ScipionParams(String args[]) {
        super(args);
    }
    
    public boolean isScipion()
    {
        return cmdLine.hasOption(SCIPION);
    }

    public void defineArgs() {
        super.defineArgs();
        Option opt = new Option(SCIPION, "");
        opt.setArgs(3);
        options.addOption(opt);
        opt = new Option(OBJECT_CMDS, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
    }

    @Override
    public void processArgs(String args[]) {
        super.processArgs(args);
        
        if (cmdLine.hasOption(SCIPION)) {
            
            String[] cmdargs = cmdLine.getOptionValues(SCIPION);
            port = Integer.parseInt(cmdargs[0]);
            inputid = cmdargs[1];
            other = (cmdargs.length == 3)? cmdargs[2]: "";
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
