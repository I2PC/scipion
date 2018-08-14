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

    

    public String inputid;
    public String other;

    public final static String SCIPION = "scipion";
    
    
    public ScipionParams(String args[]) {
        super(args);
    }
    
    public ScipionParams() {
        super();
    }
    
    public ScipionParams(Integer port, String inputid, String other) {
        this.port = port;
        this.inputid = inputid;
        this.other = other;
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
    
    public ScipionParams getScipionParams()
    {
    	return new ScipionParams(port, inputid, other);
    }
    
    
    
    
   
}
