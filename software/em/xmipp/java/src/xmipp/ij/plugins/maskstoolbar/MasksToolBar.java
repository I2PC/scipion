/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.ij.plugins.maskstoolbar;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.Toolbar;
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
import xmipp.ij.plugins.maskstoolbar.ICONS;
import xmipp.ij.plugins.maskstoolbar.LABELS;

/**
 *
 * @author Juanjo Vega
 */
public class MasksToolBar extends PlugInFrame implements ActionListener {

    private static MasksToolBar instance = null;    // To have an unique bar.
    private ImagePlus mask = null;
    private JButton jbSetBrushSize = new JButton();
    private JButton jbCreateMask = new JButton();
    private JButton jbCreateSelection = new JButton();
    private JButton jbInvertMask = new JButton();
    private JButton jbInvertSelection = new JButton();
    private JButton jbSpecifySelection = new JButton();
    private JToggleButton jtbBrush = new JToggleButton();
    private JToggleButton jtbEllipse = new JToggleButton();
    private JToggleButton jtbFreehand = new JToggleButton();
    private JToggleButton jtbLockMask = new JToggleButton();
    private JToggleButton jtbOval = new JToggleButton();
    private JToggleButton jtbPolygon = new JToggleButton();
    private JToggleButton jtbRectangle = new JToggleButton();
    private JToggleButton jtbRoundRectangle = new JToggleButton();
    private JToggleButton toolsButtons[] = new JToggleButton[]{
        jtbRectangle,
        jtbRoundRectangle,
        jtbEllipse,
        jtbOval,
        jtbPolygon,
        jtbFreehand,
        jtbBrush};

    @SuppressWarnings("LeakingThisInConstructor")
    public MasksToolBar() {
        super(LABELS.TITLE);

        setLayout(new GridLayout(2, 8));

        jtbRectangle.setIcon(ICONS.TOOL_RECTANGLE);
        jtbRoundRectangle.setIcon(ICONS.TOOL_ROUNDRECTANGLE);
        jtbEllipse.setIcon(ICONS.TOOL_ELLIPSE);
        jtbOval.setIcon(ICONS.TOOL_OVAL);
        jtbPolygon.setIcon(ICONS.TOOL_POLYGON);
        jtbFreehand.setIcon(ICONS.TOOL_FREEHAND);
        jtbBrush.setIcon(ICONS.TOOL_BRUSH);
        jbSetBrushSize.setIcon(ICONS.ACTION_BRUSHSIZE);
        jtbLockMask.setIcon(ICONS.MASK_UNLOCKED);
        jbCreateMask.setIcon(ICONS.ACTION_CREATEMASK);
        jbCreateSelection.setIcon(ICONS.ACTION_CREATESELECTION);
        jbInvertSelection.setIcon(ICONS.ACTION_INVERTSELECTION);
        jbInvertMask.setIcon(ICONS.ACTION_INVERTMASK);
        jbSpecifySelection.setIcon(ICONS.ACTION_SPECIFYSELECTION);

        jtbRectangle.setToolTipText(LABELS.TOOLTIP_RECTANGLE);
        jtbRoundRectangle.setToolTipText(LABELS.TOOLTIP_ROUNDRECTANGLE);
        jtbEllipse.setToolTipText(LABELS.TOOLTIP_ELLIPSE);
        jtbOval.setToolTipText(LABELS.TOOLTIP_OVAL);
        jtbPolygon.setToolTipText(LABELS.TOOLTIP_POLYGON);
        jtbFreehand.setToolTipText(LABELS.TOOLTIP_FREEHAND);
        jtbBrush.setToolTipText(LABELS.TOOLTIP_BRUSH);
        jbSetBrushSize.setToolTipText(LABELS.TOOLTIP_SETBRUSHSIZE);
        jtbLockMask.setToolTipText(LABELS.TOOLTIP_MASK_UNLOCKED);
        jbCreateMask.setToolTipText(LABELS.TOOLTIP_CREATEMASK);
        jbCreateSelection.setToolTipText(LABELS.TOOLTIP_CREATESELECTION);
        jbInvertSelection.setToolTipText(LABELS.TOOLTIP_INVERTSELECTION);
        jbInvertMask.setToolTipText(LABELS.TOOLTIP_INVERTMASK);
        jbSpecifySelection.setToolTipText(LABELS.TOOLTIP_SPECIFYSELECTION);

        add(jtbRectangle);
        add(jtbRoundRectangle);
        add(jtbEllipse);
        add(jtbOval);
        add(jtbPolygon);
        add(jtbFreehand);
        add(jtbBrush);
        add(jbSetBrushSize);
        add(jtbLockMask);
        add(jbCreateMask);
        add(jbCreateSelection);
        add(jbInvertSelection);
        add(jbInvertMask);
        add(jbSpecifySelection);
        add(Box.createVerticalGlue());

        jtbRectangle.addActionListener(this);
        jtbRoundRectangle.addActionListener(this);
        jtbEllipse.addActionListener(this);
        jtbOval.addActionListener(this);
        jtbPolygon.addActionListener(this);
        jtbFreehand.addActionListener(this);
        jtbBrush.addActionListener(this);
        jbSetBrushSize.addActionListener(this);
        jtbLockMask.addActionListener(this);
        jbCreateMask.addActionListener(this);
        jbCreateSelection.addActionListener(this);
        jbInvertSelection.addActionListener(this);
        jbInvertMask.addActionListener(this);
        jbSpecifySelection.addActionListener(this);
        jtbRectangle.setSelected(true);

        pack();
        setResizable(false);
    }

    static void setTool(String tool) {
        IJ.setTool(tool);
    }

    private void uselectOthers(JToggleButton source) {
        for (int i = 0; i < toolsButtons.length; i++) {
            if (toolsButtons[i] != source) {
                toolsButtons[i].setSelected(false);
            }
        }
    }

    /*
     * ij.gui.ToolBar.showBrushDialog() is protected.
     * Using java reflection it becomes accessible.
     */
    void showBrushDialog() {
        try {
            Toolbar toolbar = Toolbar.getInstance();

            Class[] classes = {};
            Method m = toolbar.getClass().getDeclaredMethod("showBrushDialog", classes);
            m.setAccessible(true);

            Object[] arguments = {};
            m.invoke(toolbar, arguments);
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }
    }

    private void createMask() {
        ImagePlus ip = IJ.getImage();

        if (ip != null) {
            Roi roi = ip.getRoi();

            if (roi != null) {
                if (mask == null || !mask.isVisible() || !jtbLockMask.isSelected()) {
                    ByteProcessor processor = new ByteProcessor(ip.getWidth(), ip.getHeight());
                    mask = new ImagePlus("Mask", processor);
                }

                mask.getProcessor().setColor(Color.WHITE);
                mask.getProcessor().fill(roi);
                mask.updateAndDraw();
                mask.show();
            } else {
                IJ.error("Area selection required.");
            }
        } else {
            IJ.error("There are no images open.");
        }
    }

    void createSelection() {
        ImagePlus ip = IJ.getImage();

        if (ip != null) {
            if (ip.getType() == ImagePlus.GRAY8) {
                IJ.run("Create Selection");
            } else {
                IJ.error("Image type is not valid: 8bit gray scale required.");
            }
        } else {
            IJ.error("There are no images open.");
        }
    }

    void invertSelection() {
        ImagePlus ip = IJ.getImage();

        if (ip != null) {
            Roi roi = ip.getRoi();

            if (roi != null) {
                IJ.run("Make Inverse");
            } else {
                IJ.error("Area selection required.");
            }
        } else {
            IJ.error("There are no images open.");
        }
    }

    private void invertMask() {
        if (mask != null) {
            Roi roi = mask.getRoi();    // Stores current roi and...
            mask.setRoi((Roi) null); // ...clears selection...
            mask.getProcessor().invert();
            mask.setRoi(roi);   // ...to restore it later.
            mask.updateAndDraw();
        } else {
            IJ.error("Mask has not been defined.");
        }
    }

    private void specifySelection() {
        ImagePlus ip = IJ.getImage();

        if (ip != null) {
            Roi roi = ip.getRoi();

            if (roi != null) {
                IJ.run("Specify...");
            } else {
                IJ.error("Area selection required.");
            }
        } else {
            IJ.error("There are no images open.");
        }
    }

    public void actionPerformed(ActionEvent ae) {
        if (ae.getSource() == jtbRectangle) {
            setTool("rectangle");
            uselectOthers(jtbRectangle);
        } else if (ae.getSource() == jtbRoundRectangle) {
            setTool("roundrect");
            uselectOthers(jtbRoundRectangle);
        } else if (ae.getSource() == jtbEllipse) {
            setTool("ellipse");
            uselectOthers(jtbEllipse);
        } else if (ae.getSource() == jtbOval) {
            setTool("oval");
            uselectOthers(jtbOval);
        } else if (ae.getSource() == jtbPolygon) {
            setTool("polygon");
            uselectOthers(jtbPolygon);
        } else if (ae.getSource() == jtbFreehand) {
            setTool("freeroi");
            uselectOthers(jtbFreehand);
        } else if (ae.getSource() == jtbBrush) {
            setTool("brush");
            uselectOthers(jtbBrush);
        } else if (ae.getSource() == jbSetBrushSize) {
            showBrushDialog();
        } else if (ae.getSource() == jtbLockMask) {
            jtbLockMask.setIcon(jtbLockMask.isSelected() ? ICONS.MASK_LOCKED : ICONS.MASK_UNLOCKED);
            jtbLockMask.setToolTipText(jtbLockMask.isSelected() ? LABELS.TOOLTIP_MASK_LOCKED : LABELS.TOOLTIP_MASK_UNLOCKED);
        } else if (ae.getSource() == jbCreateMask) {
            createMask();
        } else if (ae.getSource() == jbCreateSelection) {
            createSelection();
        } else if (ae.getSource() == jbInvertSelection) {
            invertSelection();
        } else if (ae.getSource() == jbInvertMask) {
            invertMask();
        } else if (ae.getSource() == jbSpecifySelection) {
            specifySelection();
        }
    }

    @Override
    public void windowClosed(WindowEvent e) {
        instance = null;
    }

    public static MasksToolBar getInstance() {
        if (instance == null) {
            createInstance();
        }

        return instance;
    }

    private static void createInstance() {
        instance = new MasksToolBar();
    }
}
