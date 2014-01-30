/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import xmipp.utils.Param;
import static xmipp.utils.Param.COMMAND;

/**
 *
 * @author airen
 */
public class ScipionParams extends Param {

    public String cmdname;
    public String cmdscript;
    public String dbpath;
    public String projectid;
    public String protid;

    public ScipionParams(String args[]) {
        super(args);
    }

    public void defineArgs() {
        super.defineArgs();
        Option cmdoption = new Option(COMMAND, "");
        cmdoption.setArgs(5);
        options.addOption(cmdoption);
    }

    @Override
    public void processArgs(String args[]) {
        super.processArgs(args);

        if (cmdLine.hasOption(COMMAND)) {
            String[] cmdargs = cmdLine.getOptionValues(COMMAND);
            cmdname = cmdargs[0];
            cmdscript = cmdargs[1]; 
            projectid = cmdargs[2];
            protid = cmdargs[3];
            dbpath = cmdargs[4];
        }
    }
}
