package maskstoolbar;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.plugin.PlugIn;
import ij.plugin.frame.PlugInFrame;
import ij.process.ByteProcessor;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.lang.reflect.Method;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JToggleButton;
import maskstoolbar.constants.ICONS;
import maskstoolbar.constants.LABELS;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class MasksToolBarPlugin implements PlugIn {

    @Override
    public void run(String args) {
        MasksToolBar.getInstance().setVisible(true);
    }
}
