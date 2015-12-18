package xmipp.viewer.particlepicker.training.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import xmipp.jni.Classifier;
import xmipp.jni.Classifier.Parameter;
import xmipp.utils.ColorIcon;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippQuestionDialog;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.ctf.CTFAnalyzerJFrame;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.ParticlesDialog;
import xmipp.viewer.particlepicker.training.model.ManualParticle;
import xmipp.viewer.particlepicker.training.model.MicrographState;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.ParticleToTemplatesTask;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

public class SupervisedPickerJFrame extends ParticlePickerJFrame {

    protected SupervisedPickerCanvas canvas;
    protected JMenuBar mb;
    protected SupervisedParticlePicker ppicker;
    protected JPanel micrographpn;
    protected MicrographsTableModel micrographsmd;

    protected float positionx;
    protected JButton iconbt;
    protected JLabel manuallb;
    protected JLabel autolb;
    protected JSlider thresholdsl;
    protected JPanel thresholdpn;
    protected JFormattedTextField thresholdtf;

    protected JToggleButton centerparticlebt;
    protected JMenuItem exportmi;
    protected JToggleButton autopickchb;
    protected JPanel sppickerpn;
    
    protected JMenuItem templatesmi;
    TemplatesJDialog templatesdialog;
    protected JLabel checkpercentlb;
    protected JFormattedTextField autopickpercenttf;
    protected JLabel thresholdlb;
    protected JPanel gpickerpn;
	private JButton autopickbt;
	private HashMap<JTextField, Parameter> paramtfs;

    @Override
    public SupervisedParticlePicker getParticlePicker() {
        return ppicker;
    }

    public SupervisedPickerJFrame(SupervisedParticlePicker picker) {

        super(picker);
        try {
            this.ppicker = picker;
            initComponents();
            setChanged(false);
            if(!ppicker.getClassifier().needsTraining())
            	ppicker.autopick(this, getMicrograph());
        } catch (IllegalArgumentException ex) {
            //close();
            throw ex;
        }
    }

    public boolean isCenterParticle() {
        return centerparticlebt.isSelected();
    }

    @Override
    public ParticlesDialog initParticlesJDialog() {
        return new ParticlesDialog(this);
    }

    public SupervisedPickerMicrograph getMicrograph() {
        return ppicker.getMicrograph();
    }

//    @Override
//    protected void openHelpURl() {
//        XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Micrograph_particle_picking_v3");
//
//    }

    public double getThreshold() {
        if (thresholdsl == null) {
            return 0;
        }
        return thresholdsl.getValue() / 100.0;
    }

    @Override
    public List<? extends ManualParticle> getAvailableParticles() {
        return getMicrograph().getAvailableParticles(getThreshold());
    }

    public String importParticlesFromFile(Format format, String file, float scale, boolean invertx, boolean inverty) {
        if (!new File(file).exists()) {
            throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("file", file));
        }
        String result = "";
        if (ppicker.isReviewFile(file)) {
            result = ppicker.importAllParticles(file, scale, invertx, inverty);
            ppicker.saveAllData();
        } else {
            result = importMicrographParticles(format, file, scale, invertx, inverty);
        }
        setChanged(false);
        getCanvas().repaint();
        updateMicrographsModel();
        updateSize(ppicker.getSize());//will also update templates
        canvas.refreshActive(null);

        return result;
    }

    public String importMicrographParticles(Format format, String file, float scale, boolean invertx, boolean inverty) {

        String filename = Micrograph.getName(file, 1);
		// validating you want use this file for this micrograph with different
        // name
        if (!filename.equals(getMicrograph().getName())) {
            String msg = String.format("Are you sure you want to import data from file\n%s to micrograph %s ?", file, getMicrograph().getName());
            boolean importdata = XmippDialog.showQuestion(this, msg);
            if (!importdata) {
                return null;
            }
        }
        ppicker.resetMicrograph(getMicrograph());
        String result = ppicker.importParticlesFromFile(file, format, getMicrograph(), scale, invertx, inverty);
        ppicker.saveData(getMicrograph());
        return result;
    }

    public void updateSize(int size) {
        try {
            
            super.updateSize(size);
            ppicker.updateTemplatesStack(true);
            

        } catch (Exception e) {
        	e.printStackTrace();
            XmippDialog.showError(this, e.getMessage());
        }
    }

    @Override
    public String importParticles(Format format, String dir, String preffix, String suffix, float scale, boolean invertx, boolean inverty) {
        String result = "";

        if (new File(dir).isDirectory()) {
            //System.err.println("JM_DEBUG: ============= import from Folder ============");
            result = ppicker.importParticlesFromFolder(dir, format, preffix, suffix, scale, invertx, inverty);
            boolean resize = ((Integer)sizetf.getValue()).intValue() != ppicker.getSize();
            sizetf.setValue(ppicker.getSize());
            getCanvas().repaint();
            updateMicrographsModel(true);
            getCanvas().refreshActive(null);
            ppicker.updateTemplatesStack(resize);
            

        } else // only can choose file if TrainingPickerJFrame instance
        {
            result = importParticlesFromFile(format, dir, scale, invertx, inverty);
        }
        return result;

    }

    protected void initializeCanvas() {

        if (canvas == null) {
            canvas = new SupervisedPickerCanvas(this);
        } else {
            canvas.updateMicrograph();
        }

    }

    protected void formatMicrographsTable() {
        int width = 515;
        micrographstb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        micrographstb.getColumnModel().getColumn(0).setPreferredWidth(40);
        micrographstb.getColumnModel().getColumn(1).setPreferredWidth(325);
        if(!ppicker.containsPSD())
        {
            micrographstb.getColumnModel().getColumn(1).setPreferredWidth(440);
            width = 630;
        }
        micrographstb.getColumnModel().getColumn(2).setPreferredWidth(70);
        micrographstb.getColumnModel().getColumn(3).setPreferredWidth(80);
        micrographstb.setPreferredScrollableViewportSize(new Dimension(width, 304));
        micrographstb.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        int index = ppicker.getMicrographIndex();
        if (index != -1) {
            micrographstb.setRowSelectionInterval(index, index);
        }
    }

    void updateColor() {
        color = ppicker.getColor();
        colorbt.setIcon(new ColorIcon(color));
        canvas.repaint();
        ppicker.saveConfig();
    }

    public void setChanged(boolean changed) {
        ppicker.setChanged(changed);
        savemi.setEnabled(changed);
        savebt.setEnabled(changed);
        
    }

    public void updateMicrographsModel(boolean all) {

        if (particlesdialog != null) {
            loadParticles(false);
        }

        int index = ppicker.getMicrographIndex();
        if (all) {
            micrographsmd.fireTableRowsUpdated(0, micrographsmd.getRowCount() - 1);
        } else {
            micrographsmd.fireTableRowsUpdated(index, index);
        }

        micrographstb.setRowSelectionInterval(index, index);
        manuallb.setText(Integer.toString(ppicker.getManualParticlesNumber()));
        autolb.setText(Integer.toString(ppicker.getAutomaticParticlesNumber(getThreshold())));
    }

    public ParticlePickerCanvas getCanvas() {
        return canvas;
    }

    private void initComponents() {
        try {

            setResizable(false);
            setTitle("Xmipp Particle Picker - " + ppicker.getMode());
            initMenuBar();
            setJMenuBar(mb);

            GridBagConstraints constraints = new GridBagConstraints();
            constraints.insets = new Insets(0, 5, 0, 5);
            constraints.anchor = GridBagConstraints.WEST;
            setLayout(new GridBagLayout());

            initToolBar();
            centerparticlebt = new JToggleButton(bundle.getString("center"), XmippResource.getIcon("center.png"));
            centerparticlebt.setSelected(true);
            tb.add(centerparticlebt, 0);
            add(tb, XmippWindowUtil.getConstraints(constraints, 0, 0, 2, 1, GridBagConstraints.WEST));

            //add(shapepn, XmippWindowUtil.getConstraints(constraints, 1, 1));

            
            initSupervisedPickerPane();
            if(!ppicker.getClassifier().needsTraining())
            {
            	sppickerpn.setVisible(false);
            	initGenericPickerPane();
            	add(gpickerpn, XmippWindowUtil.getConstraints(constraints, 0, 3, 2, 1, GridBagConstraints.HORIZONTAL));
            }
            else
            	add(sppickerpn, XmippWindowUtil.getConstraints(constraints, 0, 3, 1, 1, GridBagConstraints.HORIZONTAL));
            add(sppickerpn, XmippWindowUtil.getConstraints(constraints, 1, 3, 1, 1, GridBagConstraints.HORIZONTAL));
            enableSupervised(ppicker.getMode() == Mode.Supervised);
            initMicrographsPane();
            add(micrographpn, XmippWindowUtil.getConstraints(constraints, 0, 4, 2, 1, GridBagConstraints.HORIZONTAL));
            JPanel actionspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
            actionspn.add(closebt);
            actionspn.add(savebt);
            actionspn.add(saveandexitbt);
            add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 5, 2, 1, GridBagConstraints.HORIZONTAL));
            if (ppicker.getMode() == Mode.ReadOnly) {
                enableEdition(false);
            }
            pack();
            positionx = 0.9f;
            XmippWindowUtil.setLocation(positionx, 0.2f, this);
            setVisible(true);
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException(e.getMessage());
        }
    }

    void loadTemplates() {

        try {
            if (templatesdialog == null) {
                templatesdialog = new TemplatesJDialog(this);
                ppicker.setTemplatesDialog(templatesdialog);
                ParticleToTemplatesTask.setTemplatesDialog(templatesdialog);
            } else {
                templatesdialog.setVisible(true);
            }
        } catch (Exception e) {
            XmippDialog.showError(this, e.getMessage());
        }

    }

    private void initSupervisedPickerPane() {
    	
        sppickerpn = new JPanel(new FlowLayout(FlowLayout.LEFT));
        sppickerpn.setBorder(BorderFactory.createTitledBorder(bundle.getString("autopick")));
        autopickchb = new JToggleButton("Activate Training", XmippResource.getIcon("pick.png"));
        autopickchb.setSelected(ppicker.isAutopick());
        autopickchb.setText(ppicker.isAutopick()? "Deactivate Training": "Activate Training");
        autopickchb.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    boolean isautopick = autopickchb.isSelected();
                    
                    SupervisedPickerMicrograph trainmic = null;
                    if (isautopick) {
                    	trainmic = getTrainMic();
                        isautopick = trainmic != null;
                    }
                    if (isautopick) {
                        ppicker.setMode(Mode.Supervised);
                        ppicker.trainAndAutopick(SupervisedPickerJFrame.this, trainmic);
                        setTitle("Xmipp Particle Picker - " + ppicker.getMode());

                    } else if (autopickchb.isSelected()) {
                        autopickchb.setSelected(false);
                    } else {
                        boolean ismanual = XmippDialog
                                .showQuestion(SupervisedPickerJFrame.this, "After this operation automatic particles will be converted to manual and classifier training lost. Are you sure you want to continue? ");
                        if (ismanual) {
                            ppicker.setMode(Mode.Manual);
                            resetbt.setEnabled(true);
                            ppicker.saveAllData();
                            updateMicrographsModel(true);
                            if (autopickchb.isSelected()) {
                                autopickchb.setSelected(false);
                            }
                            getMicrograph().resetParticlesRectangle();
                            canvas.refreshActive(null);
                            setTitle("Xmipp Particle Picker - " + ppicker.getMode());
                        } else {
                            autopickchb.setSelected(true);
                        }
                    }

                    ppicker.saveConfig();
                    enableSupervised(autopickchb.isSelected());
                    autopickchb.setText(autopickchb.isSelected()? "Deactivate Training": "Activate Training");

                } catch (Exception ex) {
                    ex.printStackTrace();
                    XmippDialog.showError(SupervisedPickerJFrame.this, ex.getMessage());
                    autopickchb.setSelected(false);
                    return;
                }
            }
        });
        sppickerpn.add(autopickchb);
        initThresholdPane();
        sppickerpn.add(thresholdpn);
        checkpercentlb = new JLabel("Explore (%):");
        sppickerpn.add(checkpercentlb);
        autopickpercenttf = new JFormattedTextField(NumberFormat.getIntegerInstance());
        autopickpercenttf.setColumns(3);
        autopickpercenttf.setValue(ppicker.getAutopickpercent());
        autopickpercenttf.setEnabled(ppicker.getMode() == Mode.Supervised || ppicker.getMode() == Mode.Manual);
        autopickpercenttf.addActionListener(new ActionListener()
        {

                @Override
                public void actionPerformed(ActionEvent arg0)
                {

                        setAutopickPercent();

                }
        });

        sppickerpn.add(autopickpercenttf);

    }

    protected void setAutopickPercent()
    {
            if (autopickpercenttf.getValue() == null || ((Number)autopickpercenttf.getValue()).intValue() <= 0)
            {
                    XmippDialog.showInfo(this, XmippMessage.getEmptyFieldMsg("Check (%)"));
                    autopickpercenttf.setValue(getMicrograph().getAutopickpercent());
                    return;
            }

            int autopickpercent = ((Number) autopickpercenttf.getValue()).intValue();
            getMicrograph().setAutopickpercent(autopickpercent);
            ppicker.setAutopickpercent(autopickpercent);
            ppicker.saveConfig();
    }

    protected void enableSupervised(boolean selected) {

        thresholdsl.setEnabled(selected);
        thresholdlb.setEnabled(selected);
        thresholdtf.setEnabled(selected);
        sizelb.setEnabled(!selected);
        sizesl.setEnabled(!selected);// not really, if there is some micrograph
        // in sup mode size cannot be changed
        sizetf.setEnabled(!selected);
        sizelb.setEnabled(!selected);
        importmi.setEnabled(!selected);
        //autopickpercenttf.setEnabled(selected);
        
        pack();

    }

    public boolean isSupervised() {
        return autopickchb.isSelected();
    }

    protected void enableEdition(boolean isenable) {
        super.enableEdition(isenable);
        centerparticlebt.setEnabled(isenable);
        if (ppicker.getMode() != Mode.Review) {
            autopickchb.setEnabled(isenable);
            thresholdpn.setEnabled(isenable);
        }
        saveandexitbt.setEnabled(isenable);
         if(ppicker.isScipionSave())
        {
            Color color = isenable? XmippWindowUtil.firebrick: XmippWindowUtil.LIGHT_BLUE; 
            Color forecolor = isenable? Color.WHITE: Color.GRAY;
            saveandexitbt.setBackground(color);
            saveandexitbt.setForeground(forecolor);
        }
       
    }

    public void initMenuBar() {
        mb = new JMenuBar();

		// Setting menus
        exportmi = new JMenuItem("Export coordinates...", XmippResource.getIcon("export_wiz.gif"));

        exportmi.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                XmippFileChooser fc = new XmippFileChooser();
                int returnVal = fc.showOpenDialog(SupervisedPickerJFrame.this);

                try {
                    if (returnVal == XmippFileChooser.APPROVE_OPTION) {
                        File file = fc.getSelectedFile();
                        ppicker.exportParticles(file.getAbsolutePath());
                        showMessage("Export successful");
                    }
                } catch (Exception ex) {
                    showException(ex);
                }
            }
        });
        filemn.add(importmi);
        if (ppicker.getMode() != Mode.Manual) {
            importmi.setEnabled(false);
        }
        filemn.add(exportmi);
        filemn.add(exitmi);

        JMenu windowmn = new JMenu(bundle.getString("window"));

        mb.add(filemn);
        mb.add(filtersmn);
        mb.add(windowmn);
        mb.add(helpmn);
        // importffilemi.setText("Import from File...");

        windowmn.add(pmi);
        windowmn.add(ijmi);

       
        templatesmi = new JMenuItem("Templates");
        templatesmi.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent arg0) {
                loadTemplates();
            }
        });
        windowmn.add(templatesmi);

    }

   

    private void initThresholdPane() {
        thresholdpn = new JPanel();
        thresholdlb = new JLabel("Threshold:");
        thresholdpn.add(thresholdlb);
        thresholdsl = new JSlider(0, 100, 0);
        thresholdsl.setPaintTicks(true);
        thresholdsl.setMajorTickSpacing(5);
        int height = (int) thresholdsl.getPreferredSize().getHeight();
        thresholdsl.setPreferredSize(new Dimension(50, height));
        thresholdpn.add(thresholdsl);

        thresholdtf = new JFormattedTextField(NumberFormat.getNumberInstance());
        thresholdtf.setColumns(3);
        thresholdtf.setValue(0);
        thresholdpn.add(thresholdtf);
        thresholdtf.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {

                double threshold = Double.parseDouble(thresholdtf.getText()) * 100;
                setThresholdValue(threshold);

            }
        });

        thresholdsl.addChangeListener(new ChangeListener() {

            @Override
            public void stateChanged(ChangeEvent e) {
                if (thresholdsl.getValueIsAdjusting()) {
                    return;
                }
                double threshold = (double) thresholdsl.getValue() / 100;
                if (getMicrograph().getThreshold() != threshold) {
                    thresholdtf.setText(String.format("%.2f", threshold));
                    setThresholdChanges();
                }
            }
        });
    }

    private void setThresholdValue(double threshold) {
        if (Math.abs(threshold) <= 100) {
            thresholdsl.setValue((int) threshold);
            setThresholdChanges();
        }
    }

    private void setThresholdChanges() {
        // setChanged(true);
        getMicrograph().setThreshold(getThreshold());
        updateMicrographsModel();
        canvas.repaint();
        if (particlesdialog != null) {
            loadParticles(false);
        }

    }

    private void initMicrographsPane() {
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.insets = new Insets(0, 5, 0, 5);
        constraints.anchor = GridBagConstraints.NORTHWEST;
        micrographpn = new JPanel(new GridBagLayout());
        micrographpn.setBorder(BorderFactory.createTitledBorder("Micrographs"));
        
        
        if(ppicker.containsPSD())
        {
            iconbt = new JButton(Micrograph.getNoImageIcon());
            iconbt.setToolTipText("Load CTF Profile");
            iconbt.setBorderPainted(false);
            iconbt.setContentAreaFilled(false);
            iconbt.setFocusPainted(false);
            iconbt.setOpaque(false);
            iconbt.setMargin(new Insets(0, 0, 0, 0));
            iconbt.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent arg0) {
                    String psd = getMicrograph().getPSD();
                    String ctf = getMicrograph().getCTF();
                    if (psd != null && ctf != null) {
                        try {
                            new CTFAnalyzerJFrame(getMicrograph().getPSDImage(), getMicrograph().getCTF(), getMicrograph().getPSD());
                        } catch (Exception ex) {
                            Logger.getLogger(SupervisedPickerJFrame.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }

                }
            });
            micrographpn.add(iconbt, XmippWindowUtil.getConstraints(constraints, 1, 0, 1));
        }
        
        JScrollPane sp = new JScrollPane();
        micrographsmd = new MicrographsTableModel(this);
        micrographstb.setModel(micrographsmd);
        formatMicrographsTable();
        sp.setViewportView(micrographstb);
        micrographpn.add(sp, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
        JPanel infopn = new JPanel();
        manuallb = new JLabel(Integer.toString(ppicker.getManualParticlesNumber()));
        autolb = new JLabel(Integer.toString(ppicker.getAutomaticParticlesNumber(getThreshold())));
        infopn.add(new JLabel("Manual:"));
        infopn.add(manuallb);
        infopn.add(new JLabel("Automatic:"));
        infopn.add(autolb);
        micrographpn.add(infopn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
        JPanel buttonspn = new JPanel(new FlowLayout(FlowLayout.LEFT));

        buttonspn.add(resetbt);
        micrographpn.add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 2));

    }

    protected synchronized void loadMicrograph() {

        if (micrographstb.getSelectedRow() == -1) {
            return;// Probably from fireTableDataChanged raised
        }		// is same micrograph??
        int micindex = ppicker.getMicrographIndex();
        
        int index = micrographstb.getSelectedRow();
        if (micindex == index && canvas != null && canvas.getIw().isVisible()) {
            return;
        }
        

        SupervisedPickerMicrograph next = ppicker.getMicrographs().get(index);
        SupervisedPickerMicrograph current = getMicrograph();
        current.resetParticlesRectangle();
        if (!current.equals(next))// app just started
        {
            int result = tryCorrectAndAutopick(current, next);

            if (result == XmippQuestionDialog.CANCEL_OPTION) {
                micrographstb.setRowSelectionInterval(micindex, micindex);
                return;
            }
            
        }
        if (ppicker.isChanged() && ppicker.getMode() != Mode.Supervised) 
        {
            ppicker.saveData(current);
            setChanged(false);
        }
        ppicker.setMicrograph(next);
        ppicker.saveConfig();
        
        if (ppicker.getMode() == Mode.Supervised) {
            resetbt.setEnabled(next.getState() != MicrographState.Manual);
            double threshold = next.getThreshold();
            if(threshold == 0)
            {
                threshold = current.getThreshold();
                next.setThreshold(threshold);
            }
            thresholdsl.setValue((int) (threshold * 100));
            thresholdtf.setValue(threshold);
        }
        initializeCanvas();
        if(ppicker.containsPSD())
            iconbt.setIcon(next.getCTFIcon());

        pack();

        if (particlesdialog != null) {
            loadParticles(false);
        }

    }

    private int tryCorrectAndAutopick(SupervisedPickerMicrograph current, SupervisedPickerMicrograph next) {
        int result = 3;

        boolean isautopick = ppicker.getMode() == Mode.Supervised && next.getState() == MicrographState.Available;
        if (ppicker.isCorrectPending()) {
            String msg = String.format("Would you like to correct training with added and deleted particles from micrograph %s?", current.getName());
            result = XmippDialog.showQuestionYesNoCancel(this, msg);
            if (result == XmippQuestionDialog.YES_OPTION) {
                ppicker.correctAndAutopick(this, current, next);
                isautopick = false;
            }
            if (result == XmippQuestionDialog.CANCEL_OPTION)
                isautopick = false;
        }
        if (isautopick)// if not done before
        {
            ppicker.autopick(this, next);
        }
        
        return result;
    }

    protected void resetMicrograph() {
        ppicker.resetMicrograph(getMicrograph());
        canvas.refreshActive(null);
        updateMicrographsModel();
        ppicker.updateTemplatesStack(false);
        
        if (ppicker.getMode() == Mode.Supervised) {
            ppicker.autopick(this, getMicrograph());
        }
    }

    public boolean isPickingAvailable(MouseEvent e) {
        if (!super.isPickingAvailable(e)) {
            return false;
        }
        return getMicrograph().getState() != MicrographState.Corrected;
    }

    public SupervisedPickerMicrograph getTrainMic() {
        String trainmsg = "Classifier training for autopick requires that the previous micrographs and the particle's region detected to be fully picked. "
                + "Are you sure you want to continue?";
        if (ppicker.getManualParticlesNumber() < ppicker.getParticlesThreshold()) {
            XmippDialog.showInfo(this, String
                    .format("You should have at least %s particles to go to %s mode", ppicker.getParticlesThreshold(), Mode.Supervised));
            return null;
        }
        boolean isvalidrect = false;
        int index = micrographstb.getSelectedRow();
        SupervisedPickerMicrograph rectmic = null;
        while (index >= 0) 
        {
                rectmic = ppicker.getMicrographs().get(index);
                if (rectmic.hasManualParticles()) 
                {
                    rectmic.initParticlesRectangle(ppicker);
                    isvalidrect = rectmic.isValid(ppicker);
                    if(isvalidrect)
                        break;
                    else
                        rectmic.resetParticlesRectangle();
                }
            index --;
        } 
        if(index < 0)
        {
            index = micrographstb.getSelectedRow() + 1;
            while (index < ppicker.getMicrographs().size()) 
            {
                    rectmic = ppicker.getMicrographs().get(index);
                    if (rectmic.hasManualParticles()) 
                    {
                        rectmic.initParticlesRectangle(ppicker);
                        isvalidrect = rectmic.isValid(ppicker);
                        if(isvalidrect)
                            break;
                    }
                index ++;
            }
        }
        
        if(index < 0 || index == ppicker.getMicrographs().size())
        {
            canvas.repaint();
            XmippDialog.showError(this, "No valid training rectangle could be found in micrographs picked. Particles might be too few.");
            rectmic.resetParticlesRectangle();
            canvas.repaint();
            return null;
        }
        trainmsg = "";
        if(rectmic.equals(getMicrograph()))
        {
            canvas.repaint();
            trainmsg += "Classifier requires that the region detected and the micrographs previously picked to be fully picked. "
                    + "Are you sure you want to continue?";
        }
        else
        	trainmsg += "Classifier requires that the micrographs previously picked to be fully picked. "
                            + "Are you sure you want to continue?";
        int result = XmippDialog
                .showQuestionYesNoCancel(this, trainmsg);
        if (result == XmippQuestionDialog.CANCEL_OPTION) {
            canvas.repaint();
            rectmic.resetParticlesRectangle();
        }
        if(result != XmippQuestionDialog.YES_OPTION)
        	rectmic = null;
        return rectmic;
    }

    

    public void close() {
        if (ppicker.isCorrectPending()) {
            boolean iscorrect = XmippDialog.showQuestion(this, "Would you like to correct training with added and deleted particles?");
            if (iscorrect) {
                ppicker.correct();
            }
        }
        super.close();
    }

    public String getResetMsg() {
        String msg = super.getResetMsg();
        if (ppicker.getMode() == Mode.Supervised) {
            msg += "\nParticles will be automatically picked again";
        }
        return msg;
    }
    
    
    
    protected void initGenericPickerPane()
    {
    	gpickerpn = new JPanel(new FlowLayout(FlowLayout.LEFT));
    	gpickerpn.setBorder(BorderFactory.createTitledBorder("Autopick"));
    	List<Classifier.Parameter> parameters = ppicker.getClassifier().getParams();
    	paramtfs = new HashMap<JTextField, Classifier.Parameter>();
    	JTextField tf;
    	for(Classifier.Parameter param: parameters)
    	{
    		gpickerpn.add(new JLabel(param.label));
    		tf = new JTextField(3);
    		tf.setToolTipText(param.help);
    		tf.setText(param.value);
    		tf.setActionCommand(param.name);
    		tf.addFocusListener(new FocusListener()
			{
				
				@Override
				public void focusLost(FocusEvent e)
				{
					if(e.isTemporary())
						return;
					JTextField tf = (JTextField)e.getComponent();
					updateParam(tf);
				}
				
				@Override
				public void focusGained(FocusEvent e)
				{
					// TODO Auto-generated method stub
					
				}
			});
    		tf.addActionListener(new ActionListener()
			{
				
				@Override
				public void actionPerformed(ActionEvent e)
				{
					JTextField tf = (JTextField)e.getSource();
					updateParam(tf);
					
				}
			});
    		gpickerpn.add(tf);
    		paramtfs.put(tf, param);
    		if(param.name.equals("diameter") || param.name.equals("boxSize"))//classifier diameter and particle size are the same
    			updateSize(Integer.parseInt(tf.getText()));
    	}
    	autopickbt = XmippWindowUtil.getTextButton("Run", new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent arg0)
			{
				for(SupervisedPickerMicrograph m: ppicker.getMicrographs())
					ppicker.resetMicrograph(m);
				ppicker.autopick(SupervisedPickerJFrame.this, getMicrograph());
				
			}});
    	autopickbt.setBackground(XmippWindowUtil.firebrick);
    	autopickbt.setForeground(Color.WHITE);
    	gpickerpn.add(autopickbt);
    }
    
    protected void updateParam(JTextField tf)
    {
		Classifier.Parameter param = paramtfs.get(tf);
		if(param.value.equals(tf.getText()))
			return;
		if(param.name.equals("diameter") || param.name.equals("boxSize"))//classifier diameter and particle size are the same
			updateSize(Integer.parseInt(tf.getText()));
		param.value = tf.getText();
    }

}
