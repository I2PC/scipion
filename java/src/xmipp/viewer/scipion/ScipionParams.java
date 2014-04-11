/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import xmipp.utils.Param;


/**
 *
 * @author airen
 */
public class ScipionParams extends Param {

    public final static String SCIPION = "scipion";
    public String python;
    public String type;
    public String script;
    public String projectid;
    public String inputimagesid;
    public String inputid;
    

    public ScipionParams(String args[]) {
        super(args);
    }

    public void defineArgs() {
        super.defineArgs();
        Option cmdoption = new Option(SCIPION, "");
        cmdoption.setArgs(6);
        options.addOption(cmdoption);
    }

    @Override
    public void processArgs(String args[]) {
        super.processArgs(args);

        if (cmdLine.hasOption(SCIPION)) {
            String[] cmdargs = cmdLine.getOptionValues(SCIPION);
            type = cmdargs[0];
            python = cmdargs[1];
            script = cmdargs[2]; 
            projectid = cmdargs[3];
            inputid = cmdargs[4];
            inputimagesid = cmdargs[5];
        }
    }
}
