/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/**
 *
 * @author Juanjo Vega
 */
public class COSS {

    private final static String COMMAND_OPTION_FILE = "file";

    public static void main(String args[]) {
        String FILES[] = processArgs(args);

        if (FILES != null) {
            for (int i = 0; i < FILES.length; i++) {
                if (!FILES[i].isEmpty()) {
                    final String fileName = FILES[i];
                    java.awt.EventQueue.invokeLater(new Runnable() {

                        public void run() {
                            new JFrameCOSS(fileName).setVisible(true);
                        }
                    });
                }
            }
        }

        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                new JFrameCOSS("-- TEST --").setVisible(true);
            }
        });
    }

    public static String[] processArgs(String args[]) {
        Options options = new Options();

        options.addOption(COMMAND_OPTION_FILE, true, "file(s)");

        // It should be able to handle multiple files.
        options.getOption(COMMAND_OPTION_FILE).setOptionalArg(true);
        options.getOption(COMMAND_OPTION_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(COMMAND_OPTION_FILE)) {
                return cmdLine.getOptionValues(COMMAND_OPTION_FILE);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }
}
