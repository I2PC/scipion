/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.ctf;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.ImageWindow;
import ij.gui.Roi;
import ij.process.EllipseFitter;
import ij.process.ImageStatistics;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import xmipp.utils.DEBUG;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author Juanjo Vega
 */
public class CTFRecalculateImageWindow extends ImageWindow implements ActionListener, ChangeListener{

    private JButton button;
    protected EllipseFitter ellipseFitter = new EllipseFitter();
    private EllipseCTF ellipseCTF;
    private TasksEngine tasksEngine;
    private String PSDFilename, sortFn;
    private int row;
    private JSpinner spinnerLowFreq;
    private JSpinner spinnerHighFreq;
    private CTFCanvas canvas;
    
    public CTFRecalculateImageWindow(ImagePlus imp, String CTFFilename, String PSDFilename,
            TasksEngine tasksEngine, int row, String sortFn) {
        super(imp, new CTFCanvas(imp));
        

        this.PSDFilename = PSDFilename;
        this.sortFn = sortFn;
        this.tasksEngine = tasksEngine;
        this.row = row;

        ellipseCTF = new EllipseCTF(CTFFilename, imp.getWidth());

        button = XmippWindowUtil.getTextButton("Recalculate CTF", this);
        button.setEnabled(false);

        canvas = (CTFCanvas)imp.getCanvas();
        canvas.addMouseListener(new MouseListener() {

            public void mouseClicked(MouseEvent e) {
            }

            public void mousePressed(MouseEvent e) {
            }

            public void mouseReleased(MouseEvent e) {
                // If there is no ellipse ROI, button is disabled.
                button.setEnabled(fitEllipse());
            }

            public void mouseEntered(MouseEvent e) {
            }

            public void mouseExited(MouseEvent e) {
            }
        });
        
        canvas.setMaster(this);       
        

        Panel previousContent = new Panel();
        previousContent.setLayout(new ImageLayout(ic));

        Component components[] = getComponents();
        for (int i = 0; i < components.length; i++) {
            previousContent.add(components[i]);
        }

        removeAll();
        setLayout(new BorderLayout());

        add(previousContent, BorderLayout.CENTER);

        Panel panel = new Panel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(5, 5, 5, 5);
        
        
        panel.add(new JLabel("Draw an ellipse to fit the first zero"),XmippWindowUtil.getConstraints(gbc, 0, 0, 4));
        
        panel.add(new JLabel("Low freq"), XmippWindowUtil.getConstraints(gbc, 0, 1));
        spinnerLowFreq = new JSpinner(new SpinnerNumberModel(0.0, 0.0, 0.5, 0.01));
        panel.add(spinnerLowFreq, XmippWindowUtil.getConstraints(gbc, 1, 1));
        panel.add(new JLabel("High freq"), XmippWindowUtil.getConstraints(gbc, 2, 1));
        spinnerLowFreq.addChangeListener(this);        
        
        spinnerHighFreq = new JSpinner(new SpinnerNumberModel(0.5, 0.0, 0.5, 0.01));
        panel.add(spinnerHighFreq, XmippWindowUtil.getConstraints(gbc, 3, 1));
        panel.add(button,  XmippWindowUtil.getConstraints(gbc, 2, 2, 2));
        spinnerHighFreq.addChangeListener(this);   
        
        add(panel, BorderLayout.SOUTH);

        setMaximumSize(getPreferredSize());

        pack();
        imp.updateImage();
        IJ.setTool("ellipse");
    }
    
    public double getLowFreq(){
    	return ((Double)spinnerLowFreq.getValue()).doubleValue();
    }
    
    public double getHighFreq(){
    	return ((Double)spinnerHighFreq.getValue()).doubleValue();
    }

    public boolean fitEllipse() {
        boolean fitted = false;
        Roi roi = imp.getRoi();

        if (roi != null && roi.isArea()) {
            IJ.run(imp, "Fit Ellipse", "");
            ImageStatistics is = imp.getStatistics();

            ellipseFitter = new EllipseFitter();
            ellipseFitter.fit(imp.getProcessor(), is);

            // Centers ellipse in image.
            roi = imp.getRoi();
            if (roi != null) {
                Rectangle r = roi.getBounds();

                double xCenter = r.x + r.width / 2;
                double yCenter = r.y + r.height / 2;

                double dx = imp.getWidth() / 2 - xCenter;
                double dy = imp.getHeight() / 2 - yCenter;

                roi.setLocation((int) (r.x + dx), (int) (r.y + dy));
                imp.draw();

                // Store values for later use.
                fitted = true;
                ellipseCTF.calculateDefocus(ellipseFitter.minor / 2, ellipseFitter.major / 2);
                IJ.showStatus(String.format("Defocus(U,V)=%f,%f Angstroms", ellipseCTF.getDefocusU(),ellipseCTF.getDefocusV()));
            }
        }

        return fitted;
    }

    private void recalculateCTF() {
        ellipseCTF.calculateDefocus(ellipseFitter.minor / 2, ellipseFitter.major / 2);
        ellipseCTF.setFreqRange(getLowFreq(), getHighFreq());
        // Add "estimate..." to tasks.
        EstimateFromCTFTask estimateFromCTFTask = new EstimateFromCTFTask(
                ellipseCTF, 90 - ellipseFitter.angle, 
                PSDFilename, imp.getWidth(), tasksEngine, row, sortFn);
        tasksEngine.add(estimateFromCTFTask);
        dispose();
    }

	@Override
	public void actionPerformed(ActionEvent e) {
		recalculateCTF();
	}

	@Override
	public void stateChanged(ChangeEvent arg0) {
		getImagePlus().updateAndRepaintWindow();		
	}
}
