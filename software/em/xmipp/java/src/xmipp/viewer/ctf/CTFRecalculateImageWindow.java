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
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import static javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import xmipp.ij.commons.XmippApplication;
import xmipp.jni.EllipseCTF;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.GalleryData;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author Juanjo Vega
 */
public class CTFRecalculateImageWindow extends ImageWindow implements ActionListener, ChangeListener{

    private JButton recalculatebt, cancelbt;
    protected EllipseFitter ellipseFitter = new EllipseFitter();
    private EllipseCTF ellipseCTF;
    private String PSDFilename, sortFn;
    private int row;
    private JSpinner spinnerLowFreq;
    private JSpinner spinnerHighFreq;
    private CTFCanvas canvas;
    private GalleryData data;
    protected TasksEngine tasksEngine;
    private final boolean[] selection;
    
    public CTFRecalculateImageWindow(GalleryData data, boolean[] selection, ImagePlus imp, String psd, EllipseCTF ellipsectf,
             TasksEngine tasksEngine, int row, String sortFn) {
        super(imp, new CTFCanvas(imp));
        XmippApplication.addInstance(true);
        this.selection = selection;
        this.data = data;
        this.PSDFilename = psd;
        this.sortFn = sortFn;
        this.row = row;
        this.tasksEngine = tasksEngine;
        ellipseCTF = ellipsectf;

        recalculatebt = XmippWindowUtil.getTextButton("Ok", this);
        recalculatebt.setEnabled(false);
        cancelbt = XmippWindowUtil.getTextButton("Cancel", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                exit();
            }
        });

        canvas = (CTFCanvas)imp.getCanvas();
        canvas.addMouseListener(new MouseListener() {

            public void mouseClicked(MouseEvent e) {
            }

            public void mousePressed(MouseEvent e) {
            }

            public void mouseReleased(MouseEvent e) {
                // If there is no ellipse ROI, button is disabled.
                recalculatebt.setEnabled(fitEllipse());
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
        spinnerLowFreq = new JSpinner(new SpinnerNumberModel(0.05, 0.0, 0.5, 0.01));
        panel.add(spinnerLowFreq, XmippWindowUtil.getConstraints(gbc, 1, 1));
        panel.add(new JLabel("High freq"), XmippWindowUtil.getConstraints(gbc, 2, 1));
        spinnerLowFreq.addChangeListener(this);        
        
        spinnerHighFreq = new JSpinner(new SpinnerNumberModel(0.5, 0.0, 0.5, 0.01));
        panel.add(spinnerHighFreq, XmippWindowUtil.getConstraints(gbc, 3, 1));
        JPanel actionspn = new JPanel();
        actionspn.add(cancelbt);
        actionspn.add(recalculatebt);
        panel.add(actionspn,  XmippWindowUtil.getConstraints(gbc, 1, 2, 3, 1, GridBagConstraints.HORIZONTAL));
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
        if(getLowFreq() == 0)
        {
            XmippDialog.showError(null, "Low freq must be above zero");
            return;
        }
        ellipseCTF.calculateDefocus(ellipseFitter.minor / 2, ellipseFitter.major / 2);
        ellipseCTF.setFreqRange(getLowFreq(), getHighFreq());
        ellipseCTF.setEllipseFitter(ellipseFitter);
        // Add "estimate..." to tasks.
        data.recalculateCTF(row, selection, ellipseCTF, sortFn);
        exit();
    }

	@Override
	public void actionPerformed(ActionEvent e) {
		recalculateCTF();
	}

	@Override
	public void stateChanged(ChangeEvent arg0) {
		getImagePlus().updateAndRepaintWindow();		
	}
        

        public void exit()
        {
            setVisible(false);
            dispose();
            XmippApplication.removeInstance(true);
        }
        
        @Override
	public void windowClosing(WindowEvent e) 
        {
            
            super.windowClosing(e);
            XmippApplication.removeInstance(true);
	}

}
