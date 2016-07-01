package xmipp.viewer.particlepicker;

import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.plugin.frame.Recorder;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.net.URL;
import java.net.URLClassLoader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import xmipp.ij.commons.InputFieldsMessageDialog;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippApplication;
import xmipp.ij.commons.XmippImageJ;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.Filename;
import xmipp.utils.ColorIcon;
import xmipp.utils.QuickHelpJDialog;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippQuestionDialog;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.extract.ExtractPickerJFrame;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerCanvas;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.model.Mode;

public abstract class ParticlePickerJFrame extends JFrame implements ActionListener
{

    public static final int TOGGLE_LINEAR_MODE_KEY = KeyEvent.VK_L;
    public static final int TOGGLE_ERASE_MODE_KEY = KeyEvent.VK_E;
    protected ParticlesDialog particlesdialog;

	protected JMenuItem ijmi;
	protected JToggleButton circlechb;
	protected JToggleButton rectanglechb;
	protected JFormattedTextField sizetf;
	protected JToggleButton centerchb;
	protected JPanel shapepn;
	protected JMenuItem savemi;
	protected JMenuItem hcontentsmi;
	protected JMenuItem pmi;
	// protected JMenuItem importffilemi;
	protected JMenuItem exportmi;
	protected JMenu filtersmn;
	protected String activefilter;
	protected JSlider sizesl;
	protected JPanel sizepn;
	private List<JCheckBoxMenuItem> mifilters;
	protected JMenu filemn;
	protected JMenuItem importmi;
	protected JButton colorbt;
	protected Color color;
	protected JPanel colorpn;
	protected JButton resetbt;
	protected JTable micrographstb;
	protected ImportParticlesJDialog importpjd = null;
	protected JMenuItem exitmi;
	protected JLabel positionlb;
	protected JToggleButton usezoombt;
	private JToggleButton eraserbt;

    protected NoneSelectedButtonGroup modesgrp;

    // Linear picking
    private JToggleButton linearPickingbt;
    private JPanel steppn;
    private JLabel steplb;
    private JFormattedTextField steptf;

	private JMenuItem keyassistmi;
	protected JMenu helpmn;
	protected JButton savebt;
	protected JButton saveandexitbt;
	protected JToolBar tb;
    protected ResourceBundle bundle;
    protected String command;
    protected JButton closebt;

	protected JLabel sizelb;


	public ParticlePickerJFrame(ParticlePicker picker)
	{
            try {
                File file = new File(Filename.getXmippPath("resources"));
                URL[] urls = new URL[]{file.toURI().toURL()};
                ClassLoader loader = new URLClassLoader(urls);
                bundle = ResourceBundle.getBundle("Bundle", Locale.getDefault(), loader);
                XmippApplication.addInstance(false);
                setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
                addWindowListener(new WindowAdapter()
                {
                    public void windowClosing(WindowEvent winEvt)
                    {

                        if (canExit()) close();
                        
                    }
                });
                
                initMenuBar(picker);
                
                resetbt = XmippWindowUtil.getTextButton(bundle.getString("resetmic"), new ActionListener()
                {
                    
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        String resetmsg = getResetMsg();
                        XmippQuestionDialog qd = new XmippQuestionDialog(ParticlePickerJFrame.this, resetmsg, false);
                        if (qd.showDialog())
                            resetMicrograph();
                    }
                });
                
                savebt = XmippWindowUtil.getTextIconButton(bundle.getString("save"), "save.gif", new ActionListener()
                {
                    
                    @Override
                    public void actionPerformed(ActionEvent arg0)
                    {
                        saveData(false);

                    }
                });
                
                saveandexitbt = XmippWindowUtil.getTextButton(bundle.getString("saveandexit"), new ActionListener()
                {
                    
                    @Override
                    public void actionPerformed(ActionEvent e)
                            
                    {
                        if (getParticlePicker().getMode() != Mode.ReadOnly)
                            getParticlePicker().saveData();
                        if(getParticlePicker().isScipionSave())
                        {
                            int count = getParticlePicker().getParticlesCount();
                            if(count == 0)
                            {
                                XmippDialog.showInfo(ParticlePickerJFrame.this, XmippMessage.getEmptyFieldMsg("coordinates"));
                                return;
                            }
                            HashMap<String, String> msgfields = new HashMap<String, String>();
                            boolean createprot = getParticlePicker().getPort() == null;
                            if(createprot)
                                msgfields.put("Run name:", "ProtUserCoordinates");
                            
                            String msg = String.format("<html>Are you sure you want to register a new set of Coordinates with <font color=red>%s</font> %s?", count, (count != 1)?"elements":"element");
                            InputFieldsMessageDialog dlg = new InputFieldsMessageDialog(ParticlePickerJFrame.this, "Question", msg);
                            
                            if (dlg.action == InputFieldsMessageDialog.OK_OPTION)
                                executeScipionSaveAndExit();
                            
                        }
                        else
                            close();
                        
                    }
                });
                closebt = XmippWindowUtil.getTextIconButton("Close", "fa-times.png", new ActionListener()
                {
                    
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        if (canExit()) close();
                    }
                });
                if(picker.isScipionSave())
                {
                    saveandexitbt.setText("Coordinates");
                    Image img = Toolkit.getDefaultToolkit().getImage(Filename.getXmippPath("resources" + File.separator + "fa-plus-circle.png"));
                    saveandexitbt.setIcon(new ImageIcon(img));
                    saveandexitbt.setToolTipText("Create Coordinates");
                    Color color = XmippWindowUtil.firebrick; 
                    saveandexitbt.setBackground(color);
                    saveandexitbt.setForeground(Color.WHITE);
                    
                }
                micrographstb = new JTable();
                micrographstb.getSelectionModel().addListSelectionListener(new MicrographsSelectionListener());
                micrographstb.addMouseListener(new MouseListener()
                {
                    
                    @Override
                    public void mouseReleased(MouseEvent arg0)
                    {
                        // TODO Auto-generated method stub
                        
                    }
                    
                    @Override
                    public void mousePressed(MouseEvent arg0)
                    {
                        // TODO Auto-generated method stub
                        
                    }
                    
                    @Override
                    public void mouseExited(MouseEvent arg0)
                    {
                        // TODO Auto-generated method stub
                        
                    }
                    
                    @Override
                    public void mouseEntered(MouseEvent arg0)
                    {
                        // TODO Auto-generated method stub
                        
                    }
                    
                    @Override
                    public void mouseClicked(MouseEvent arg0)
                    {
                        if (micrographstb.getSelectedRow() == -1)
                            return;
                        loadMicrograph();
                    }
                });
                
                XmippWindowUtil.setScipionImageIcon(this);
                
            } catch (Exception ex) {
                Logger.getLogger(ParticlePickerJFrame.class.getName()).log(Level.SEVERE, null, ex);
                throw new IllegalArgumentException(ex);
            }
	}

    protected boolean canExit() {

        if (getParticlePicker().getParams().tmp && (getParticlePicker().isSaved() || getParticlePicker().isChanged()))
        {
            XmippQuestionDialog qd = new XmippQuestionDialog(ParticlePickerJFrame.this, "You are picking in a temporary folder. None of the saved or changed data will remain unless you export \"Coordinates\". Do you really want to exit?", false );
            boolean exit = qd.showDialog();
            if (!exit)
                return false;

        }else if (getParticlePicker().isChanged())
        {
            XmippQuestionDialog qd = new XmippQuestionDialog(ParticlePickerJFrame.this, "Save changes before closing?");
            boolean save = qd.showDialog();
            if (save)
                return saveData(true);
            else if (qd.isCanceled())
                return false;
        }
        return true;
    }

    protected class MicrographsSelectionListener implements ListSelectionListener
        {

			@Override
			public void valueChanged(ListSelectionEvent e)
			{
				if (e.getValueIsAdjusting())
					return;
                                
				if (micrographstb.getSelectedRow() == -1)
					return;// Probably from fireTableDataChanged raised
				loadMicrograph();
			}
        }

	protected abstract void loadMicrograph();
        
    public String getResetMsg()
    {
        return "Are you sure you want to remove all particles from micrograph?";
    }

    private ParticlePickerJFrame getMyself(){
        return this;
    }
	private void initMenuBar(ParticlePicker picker)
	{
		filemn = new JMenu(bundle.getString("file"));
		helpmn = new JMenu(bundle.getString("help"));
		savemi = new JMenuItem("Save", XmippResource.getIcon("save.gif"));
		savemi.setMnemonic('S');
		savemi.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK));
		savemi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
                saveData(true);
            }
		});
		filemn.add(savemi);
		importmi = new JMenuItem("Import coordinates...", XmippResource.getIcon("import_wiz.gif"));
		filemn.add(importmi);
		if (picker.getMode() != Mode.Manual)
			importmi.setEnabled(false);

		importmi.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				if (importpjd == null)
					importpjd = new ImportParticlesJDialog(ParticlePickerJFrame.this);
				importpjd.showDialog();

			}
		});

		exitmi = new JMenuItem("Exit");
		exitmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				close();
			}
		});

		ijmi = new JMenuItem("ImageJ", XmippResource.getIcon("ij.gif"));
		ijmi.setEnabled(picker.getMode() != Mode.ReadOnly);
		ijmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippUtil.showImageJ(Tool.PICKER);
			}
		});

		hcontentsmi = new JMenuItem("Online help", XmippResource.getIcon("online_help.gif"));
		hcontentsmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					openHelpURl();

				}
				catch (Exception ex)
				{
					showException(ex);

				}
			}
		});
		helpmn.add(hcontentsmi);

		keyassistmi = new JMenuItem("Tips...", XmippResource.getIcon("bulb.png"));
		keyassistmi.addActionListener(new ActionListener()
		{

			private QuickHelpJDialog keyassistdlg;

			@Override
			public void actionPerformed(ActionEvent e)
			{
				if (keyassistdlg == null)
					keyassistdlg = new QuickHelpJDialog(ParticlePickerJFrame.this, false, "Tips", getKeyAssist());
				keyassistdlg.setVisible(true);

			}
		});
		helpmn.add(keyassistmi);

		pmi = new JMenuItem("Particles", XmippResource.getIcon("table_view.gif"));
		pmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				loadParticles(false);
			}
		});

		mifilters = new ArrayList<JCheckBoxMenuItem>();
		filtersmn = new JMenu(bundle.getString("filters"));
		filtersmn.addMenuListener(new MenuListener()
		{

			@Override
			public void menuCanceled(MenuEvent arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void menuDeselected(MenuEvent arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void menuSelected(MenuEvent arg0)
			{
				boolean added;
				for (JCheckBoxMenuItem mi : mifilters)
				{
					added = getParticlePicker().isFilterAdded(mi.getText());
					mi.setSelected(added);
				}

			}
		});
                
		addFilterMenuItem(XmippImageJ.gaussianBlurFilter, true, picker);
		addFilterMenuItem(ParticlePicker.xmippsmoothfilter, true, picker);
		addFilterMenuItem(XmippImageJ.bandPassFilter, true, picker);
		addFilterMenuItem(XmippImageJ.enhanceContrastFilter, true, picker);
		addFilterMenuItem(XmippImageJ.brightnessContrastFilter, true, picker);

//		JCheckBoxMenuItem admi = addFilterMenuItem(XmippImageJ.anisotropicDiffFilter, false, picker);
//		admi.addActionListener(new ActionListener()
//		{
//
//			@Override
//			public void actionPerformed(ActionEvent e)
//			{
//				activefilter = "8-bit";
//				IJ.run(activefilter);
//				activefilter = ((JCheckBoxMenuItem) e.getSource()).getText();
//				IJ.run(activefilter);
//			}
//		});
		addFilterMenuItem(XmippImageJ.invertLUTFilter, true, picker);
		addFilterMenuItem(XmippImageJ.substractBackgroundFilter, true, picker);
        addFilterAppliedListener();
	}

    protected boolean saveData(boolean notify) {
        try {
            getParticlePicker().saveData();
            if (notify) showMessage("Data saved successfully");
            setChanged(false);
            return true;
        } catch (Exception e1) {
            XmippDialog.showException(getMyself(), e1);
            return false;
        }
    }

    protected void addFilterAppliedListener() {

        Recorder.record = true;

        // detecting if a command is thrown by ImageJ
        
        ImagePlus.addImageListener(new ImageListener() {

            @Override
            public void imageUpdated(ImagePlus imp) {
                if(command != null)
                {
                    getParticlePicker().updateFilters(command);
                    if(particlesdialog != null)
                        loadParticles(true);
                    
                }
                command = null;
            }

            @Override
            public void imageOpened(ImagePlus arg0) {

            }

            @Override
            public void imageClosed(ImagePlus arg0) {
                // TODO Auto-generated method stub

            }
        });
    }

	protected void openHelpURl()
	{
		XmippWindowUtil.openURI("https://github.com/biocompwebs/scipion/wiki/Picker");
	}

	protected abstract void resetMicrograph();

	protected void enableEdition(boolean enable)
	{
		importmi.setEnabled(enable);
		savemi.setEnabled(enable);
		sizesl.setEnabled(enable);
		sizetf.setEnabled(enable);
		colorbt.setEnabled(enable);
		resetbt.setEnabled(enable);
		savebt.setEnabled(enable);
		eraserbt.setEnabled(enable);

	}

	private JCheckBoxMenuItem addFilterMenuItem(String command, boolean defaultlistener, ParticlePicker picker)
	{
		JCheckBoxMenuItem mi = new JCheckBoxMenuItem(command);
		mifilters.add(mi);
		mi.setSelected(picker.isFilterAdded(command));
		if (defaultlistener)
			mi.addActionListener(this);
		filtersmn.add(mi);
		//mi.setEnabled(picker.getMode() != Mode.ReadOnly);
		return mi;
	}

	@Override
	public void actionPerformed(ActionEvent e)
	{
		try
		{
			JCheckBoxMenuItem item = (JCheckBoxMenuItem) e.getSource();
			activefilter = item.getText();
			if (item.isSelected())// filter added, will be registered by picker
									// with options if needed
				if (activefilter.equals(ParticlePicker.xmippsmoothfilter))
				{
					getParticlePicker().addFilter(ParticlePicker.xmippsmoothfilter, "xmipp");
					reloadImage();
				}
				else
				{
					command = activefilter;
					applyFilter(activefilter);
                                       
				}
			else
			{
				// filter removed
				getParticlePicker().removeFilter(activefilter);
				reloadImage();
                if(particlesdialog != null)
                    loadParticles(true);    
                                
			}
            if(getParticlePicker().getMode() != Mode.ReadOnly)
                getParticlePicker().saveConfig();
		}
		catch (Exception ex)
		{

			ex.printStackTrace();
			showException(ex);
		}

	}
	
	public void applyFilter(String filter)
	{
		IJ.run(getMicrograph().getImagePlus(), filter, "");
	}

	protected void reloadImage()
	{
		getCanvas().updateMicrograph();
                
	}

	public int getSide()
	{
		return 100;
	}

	public abstract ParticlePickerCanvas getCanvas();

	public abstract ParticlesDialog initParticlesJDialog();

	public void loadParticles(boolean reset)
	{
		try
		{
            if(reset)
                for (PickerParticle p : getAvailableParticles())
                    p.resetParticleCanvas();
			if (particlesdialog == null)
				particlesdialog = initParticlesJDialog();
			else
			{

				particlesdialog.loadParticles(false);
				particlesdialog.setVisible(true);
			}
		}
		catch (Exception ex)
		{
			showException(ex);
			if (particlesdialog != null)
				particlesdialog.close();
			particlesdialog = null;
		}

	}

	public void updateMicrographsModel()
	{
		updateMicrographsModel(false);

	}

	public abstract void updateMicrographsModel(boolean all);

	public ParticlesDialog getParticlesJDialog()
	{
		return particlesdialog;
	}

	public abstract Micrograph getMicrograph();

	public abstract List<? extends PickerParticle> getAvailableParticles();

	public boolean isPickingAvailable(MouseEvent e)
	{
		if (getCanvas().getTool() != Tool.PICKER)
			return false;
		if (SwingUtilities.isRightMouseButton(e))
			return false;
		if (getParticlePicker().getMode() == Mode.ReadOnly)
			return false;
		if (e.isControlDown())
			return false;
		return true;
	}

	

	public class ColorActionListener implements ActionListener
	{
		JColorChooser colorChooser;

		@Override
		public void actionPerformed(ActionEvent e)
		{
			// Set up the dialog that the button brings up.
			colorChooser = new JColorChooser();
			JDialog dialog = JColorChooser.createDialog(colorbt, "Pick a Color", true, // modal
			colorChooser, new ActionListener()
			{

				@Override
				public void actionPerformed(ActionEvent e)
				{
					updateColor(colorChooser.getColor());
					getParticlePicker().setColor(colorChooser.getColor());
				}
			}, // OK button handler
			null); // no CANCEL button handler
			XmippWindowUtil.setLocation(0.5f, 0.5f, dialog);
			dialog.setVisible(true);
		}
	}

	public void updateColor(Color color)
	{
		if (colorbt != null)
			colorbt.setIcon(new ColorIcon(color));
		getParticlePicker().setColor(color);
		getCanvas().repaint();
		getParticlePicker().saveConfig();
	}

	protected void initShapePane()
	{

		shapepn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		ShapeItemListener shapelistener = new ShapeItemListener();
        shapepn.add(new JLabel("Shape:"));
        Icon icon = XmippResource.getIcon("circle.png");
		circlechb = new JToggleButton(icon);
		//circlechb.setSelected(true);
		circlechb.addItemListener(shapelistener);
                
		rectanglechb = new JToggleButton(XmippResource.getIcon("square.png"));
        rectanglechb.setPreferredSize(null);
		rectanglechb.setSelected(true);
		rectanglechb.addItemListener(shapelistener);

		centerchb = new JToggleButton(XmippResource.getIcon("plus.png"));
		centerchb.setSelected(true);
		centerchb.addItemListener(shapelistener);

		shapepn.add(circlechb);
		shapepn.add(rectanglechb);
		shapepn.add(centerchb);

	}

	public void initToolBar()
	{
		tb = new JToolBar();

		tb.setFloatable(false);
                
		usezoombt = new JToggleButton("-1", XmippResource.getIcon("zoom.png"));
		usezoombt.setToolTipText("Keep zoom");
		usezoombt.setFocusable(false);
        usezoombt.setSelected(true);
        usezoombt.addActionListener(new ActionListener()
		{
			
			@Override
			public void actionPerformed(ActionEvent e)
			{
				if(!usezoombt.isSelected())
					getParticlePicker().setZoom(-1);
			}
		});
		tb.add(usezoombt);
        initShapePane();
        tb.add(shapepn);
		initSizePane();
		tb.add(sizepn);
		if (!(this instanceof ExtractPickerJFrame))
		{
			initColorPane(getParticlePicker().getColor());
			tb.add(colorpn);
            modesgrp = new NoneSelectedButtonGroup();
			eraserbt = new JToggleButton(bundle.getString("eraser"), XmippResource.getIcon("eraser.png"));
            modesgrp.add(eraserbt);
            eraserbt.setMnemonic(TOGGLE_ERASE_MODE_KEY);
			tb.add(eraserbt);
		}

        if ((this instanceof SupervisedPickerJFrame))
        {
            linearPickingbt = new JToggleButton(bundle.getString("linearPicking"), XmippResource.getIcon("linearPicking.png"));
            linearPickingbt.setMnemonic(TOGGLE_LINEAR_MODE_KEY);
            modesgrp.add(linearPickingbt);
            tb.add(linearPickingbt);
            initStepPane();
            tb.add(steppn);
        }
	}

	public boolean isEraserMode()
	{
		if (!eraseAvailable())
			return false;
		return eraserbt.isSelected();
	}

    protected boolean eraseAvailable() {
        return eraserbt != null;
    }

    public void toogleEraseMode(){
        if (eraseAvailable()){
            eraserbt.setSelected(!eraserbt.isSelected());
            getCanvas().setCustomCursor();
        }
    }

    public boolean isLinearMode() {
        if (linearModeAvailable()){
            return linearPickingbt.isSelected();
        } else {
            return false;
        }
    }

    public void toogleLinearMode(){
        if (linearModeAvailable()){
            linearPickingbt.setSelected(!linearPickingbt.isSelected());
            ((SupervisedPickerCanvas)getCanvas()).resetLinearPicking();
            getCanvas().setCustomCursor();
        }
    }


    protected boolean linearModeAvailable() {
        return linearPickingbt != null;
    }

    public int getStep(){
        if (linearModeAvailable()){
            return Integer.parseInt(steptf.getText());
        } else {
            return 1000;
        }
    }

	protected void displayZoom(double zoom)
	{
		
		usezoombt.setText(String.format(Locale.US, "%.2f", zoom));
		if(usezoombt.isSelected())
			getParticlePicker().setZoom(zoom);
		pack();
	}

	class ShapeItemListener implements ItemListener
	{
		@Override
		public void itemStateChanged(ItemEvent e)
		{
			changeShapes();
		}
	}

	public void changeShapes()
	{
		getCanvas().repaint();

	}

	public boolean isShapeSelected(Shape shape)
	{
		switch (shape)
		{
		case Rectangle:
			return rectanglechb.isSelected();
		case Circle:
			return circlechb.isSelected();
		case Center:
			return centerchb.isSelected();
			// case OnlyLast:
			// return onlylastchb.isSelected();
		}
		return false;
	}

	public abstract ParticlePicker getParticlePicker();

	public abstract void setChanged(boolean changed);

    protected void setTitle(){

        String title = "Xmipp Particle Picker - " + getParticlePicker().getMode();

        if (this.getParticlePicker().getParams().tmp){
            title = title + " - in Temporal folder.";
        }
        setTitle(title);
    }
    public void setSaved(boolean saved){
    };

	protected void initColorPane(Color color)
	{
		colorpn = new JPanel();
		this.color = color;
		colorpn.add(new JLabel("Color:"));
		colorbt = new JButton();
		colorbt.setContentAreaFilled(false);
		colorbt.setFocusPainted(false);
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
                colorbt.setMargin(new Insets(0, 0, 0, 0));
		colorbt.addActionListener(new ColorActionListener());
		colorpn.add(colorbt);
	}

    protected void initStepPane() {
        steppn = new JPanel();

        int size = getParticlePicker().getSize();
        steplb = new JLabel("Step:");
        steppn.add(steplb);

        steptf = new JFormattedTextField(NumberFormat.getIntegerInstance());
        steptf.setColumns(3);
        steptf.setValue(100);
        steppn.add(steptf);
    }


    protected void initSizePane()
	{
		sizepn = new JPanel();

		int size = getParticlePicker().getSize();
		sizelb = new JLabel("Size:");
		sizepn.add(sizelb);
		sizesl = new JSlider(10, ParticlePicker.sizemax, size);
		sizesl.setPaintTicks(true);
		sizesl.setMajorTickSpacing(100);
		int height = (int) sizesl.getPreferredSize().getHeight();
		sizesl.setPreferredSize(new Dimension(50, height));
		sizepn.add(sizesl);
		sizetf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		sizetf.setColumns(3);
		sizetf.setValue(size);
		sizepn.add(sizetf);
		sizetf.addFocusListener(new FocusListener()
		{

                    @Override
                    public void focusGained(FocusEvent fe) {
                    }

                    @Override
                    public void focusLost(FocusEvent fe) {
                        // event from sizes
                        try {
                            if(!fe.isTemporary())
                            {

                                    sizetf.commitEdit();
                                    readSizeFromTextField();
                            }
                        } catch (Exception ex) {
                            XmippDialog.showError(ParticlePickerJFrame.this, XmippMessage.getIllegalValueMsg("size", sizetf.getText()));
                            Logger.getLogger(ParticlePickerJFrame.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
		});
        sizetf.addActionListener(new ActionListener(){

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        // event from sizes
                            readSizeFromTextField();
                    }
		});

		sizesl.addChangeListener(new ChangeListener()
		{

			@Override
			public void stateChanged(ChangeEvent e)
			{
				if (sizesl.getValueIsAdjusting()) 
					return;
                                int size = sizesl.getValue();
				if (size == getParticlePicker().getSize())
					return;
				
				if (!getParticlePicker().isValidSize(ParticlePickerJFrame.this, size))
				{
					sizesl.dispatchEvent(//trick to repaint slider after changing value
                                            new MouseEvent(sizesl, MouseEvent.MOUSE_RELEASED,0,0,0,0,1,false));
					int prevsize = getParticlePicker().getSize();
					sizesl.setValue(prevsize);
					return;
				}
				updateSize(size);
			}
		});

	}
        
    protected void readSizeFromTextField()
    {
        
        int size = ((Number) sizetf.getValue()).intValue();
        if (size == getParticlePicker().getSize())
            return;

        if (!getParticlePicker().isValidSize(ParticlePickerJFrame.this, size))
        {

            int prevsize = getParticlePicker().getSize();
            sizetf.setText(Integer.toString(prevsize));
            return;
        }
        updateSize(size);
       
    }

	public void updateSize(int size)
	{
		
		sizetf.setValue(size);
		sizesl.setValue(size);
		if(getCanvas() != null)
			getCanvas().repaint();
		getParticlePicker().setSize(size);
                //updateMicrographsModel();
                if(particlesdialog != null)
                    loadParticles(true);
		getParticlePicker().saveConfig();
	}

	/** Shortcut function to show messages */
	public boolean showMessage(String message)
	{
		return XmippDialog.showInfo(this, message);
	}

	public boolean showException(Exception e)
	{
		return XmippDialog.showException(this, e);
	}

	public void close()
	{
		setVisible(false);
		dispose();
		if (getCanvas() != null)
			getCanvas().getIw().close();
		XmippApplication.removeInstance(false);
	}

	public abstract String importParticles(Format format, String dir, String preffix, String suffix, float scale, boolean invertx, boolean inverty);

	public Map<Object, Object> getKeyAssist()
	{
		Map<Object, Object> map = Collections.synchronizedMap(new LinkedHashMap<Object, Object>());
		map.put("Shift + scroll up / ctrl + left click", "Zoom in");
		map.put("Shift + scroll down / ctrl + right click", "Zoom out");
		map.put("Right click + mouse move", "Moves image previously expanded");
		map.put("Left click", "Adds or selects a particle. In erase mode, deletes or disables selected particles");
		map.put("Shift + left click", "Deletes or disables selected particle");
		map.put("Left click + mouse move", "Moves selected particle. In erase mode, deletes or disables particles");
		map.put("Left", "Moves selected particle to the left");
		map.put("Right", "Moves selected particle to the right");
		map.put("Up", "Moves selected particle up");
		map.put("Down", "Moves selected particle down");
        map.put((char)TOGGLE_ERASE_MODE_KEY, "Toggle between erase and picking");

        if (linearModeAvailable()) {

            map.put((char)SupervisedPickerCanvas.ERASER_MEDIUM_SIZE_KEY+ "(while erasing)", "Use a Medium size eraser");
            map.put((char)SupervisedPickerCanvas.ERASER_LARGE_SIZE_KEY+ "(while erasing)", "Use a Large size eraser");
            map.put((char)SupervisedPickerCanvas.ERASER_X_LARGE_SIZE_KEY+ "(while erasing)", "Use a X-large size eraser");

            map.put((char)TOGGLE_LINEAR_MODE_KEY, "Toggle between normal picking and linear picking mode");
            map.put((char)SupervisedPickerCanvas.POLIGONAL_MODE_KEY + " (continuous)", "While in linear picking mode, linear picking will be continuous (polygonal mode).");
            map.put("ESCAPE", "Will reset the linear picking canceling the pending line");
        }

		return map;
	}
        
    protected void executeScipionSaveAndExit()
    {
        getCanvas().setEnabled(false);
        XmippWindowUtil.blockGUI(ParticlePickerJFrame.this, "Creating set ...");
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {
                    String cmd = String.format("run function registerCoords '%s'", getParticlePicker().getOutputDir());
                    XmippWindowUtil.runCommand(cmd, getParticlePicker().getParams().port);
                    XmippWindowUtil.releaseGUI(ParticlePickerJFrame.this.getRootPane());
                    getCanvas().setEnabled(true);
                    
                    close();

                } catch (Exception ex) {
                    ex.printStackTrace();
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }

}

class NoneSelectedButtonGroup extends ButtonGroup {

    @Override
    public void setSelected(ButtonModel model, boolean selected) {

        if (selected) {

            super.setSelected(model, selected);

        } else {

            clearSelection();
        }
    }
}
