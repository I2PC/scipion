package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import browser.windows.ImagesWindowFactory;

import model.Constants;
import model.Family;
import model.Micrograph;
import model.PPData;
import model.Particle;
import model.XmippJ;

enum Tool {
	IMAGEJ, PICKER
}

enum Shape {
	CIRCLE, RECTANGLE, BOTH
}

public class XmippParticlePickerJFrame extends JFrame implements ActionListener {

	private JButton infobt;
	private JRadioButton pickrbt;
	private JSlider sizesl;
	private JRadioButton circlerbt;
	private JRadioButton rectanglerbt;
	private ButtonGroup shapebg;
	private PPCanvas canvas;
	private JTextField sizetf;
	private Tool tool = Tool.PICKER;
	private JButton adjustbt;
	private JRadioButton bothrbt;
	private Shape shape;
	private JMenuBar mb;
	private JComboBox familiescb;
	private JButton geditbt;
	private PPData ppdata;
	private JLabel colorlb;
	private ColorIcon coloricon;
	private Color color;
	private JPanel familypn;
	private JPanel symbolpn;
	private String activemacro;
	private JPanel micrographpn;
	private JList mglist;
	private ImageWindow iw;
	private boolean changed;

	public Shape getShape() {
		return shape;
	}

	public Tool getTool() {
		return tool;
	}

	

	public PPData getPPData() {
		return ppdata;
	}

	public Family getFamily() {
		return (Family) familiescb.getSelectedItem();
	}

	public XmippParticlePickerJFrame() {

		ppdata = PPData.getInstance();

		initComponents();
		initializeCanvas();
	}
	
	public Micrograph getMicrograph()
	{
		return (Micrograph)mglist.getSelectedValue();
	}



	private void initComponents() {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent winEvt) {
				if(changed)
				{
					int result = JOptionPane.showConfirmDialog(XmippParticlePickerJFrame.this, "Save changes before closing?", "Message", JOptionPane.YES_NO_OPTION);
					if(result == JOptionPane.OK_OPTION)
						XmippParticlePickerJFrame.this.saveChanges();
				}
				System.exit(0);
			}

		});
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Particle Picker");
		mb = getPPMenuBar();
		setJMenuBar(mb);

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;
		setLayout(new GridBagLayout());

		JPanel pickerpn = new JPanel(new GridBagLayout());
		pickerpn.setBorder(BorderFactory.createTitledBorder("Picker"));

		initFamilyPane();
		pickerpn.add(new JLabel("Family:"),
				WindowUtils.updateConstraints(constraints, 0, 0, 1));
		pickerpn.add(familypn,
				WindowUtils.updateConstraints(constraints, 1, 0, 1));

		initSymbolPane();
		pickerpn.add(new JLabel("Symbol:"),
				WindowUtils.updateConstraints(constraints, 0, 1, 1));
		pickerpn.add(symbolpn,
				WindowUtils.updateConstraints(constraints, 1, 1, 1));

		add(pickerpn, WindowUtils.updateConstraints(constraints, 0, 0, 3));

		initMicrographsPane();
		add(micrographpn, WindowUtils.updateConstraints(constraints, 0, 2, 3));

		pack();
		WindowUtils.centerScreen(0.9, this);
		setVisible(true);
	}

	private void initSymbolPane() {

		symbolpn = new JPanel();
		shapebg = new ButtonGroup();
		circlerbt = new JRadioButton("Circle");
		shapebg.add(circlerbt);
		rectanglerbt = new JRadioButton("Rectangle");
		shapebg.add(rectanglerbt);
		bothrbt = new JRadioButton("Both");
		bothrbt.getModel().setSelected(true);
		shape = shape.BOTH;
		shapebg.add(bothrbt);

		symbolpn.add(bothrbt);
		symbolpn.add(circlerbt);
		symbolpn.add(rectanglerbt);

		circlerbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				canvas.repaint();
				shape = Shape.CIRCLE;
			}
		});

		rectanglerbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				canvas.repaint();
				shape = Shape.RECTANGLE;
			}
		});

		bothrbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				canvas.repaint();
				shape = Shape.BOTH;
			}
		});

	}



	private void initFamilyPane() {
		familypn = new JPanel();
		
		// Setting edit button
		geditbt = new JButton("Edit");
		familypn.add(geditbt);

		// Setting combo
		familiescb = new JComboBox(ppdata.getFamilies().toArray());
		familypn.add(familiescb);

		// Setting color
		color = getFamily().getColor();
		familypn.add(new JLabel("Color:"));
		coloricon = new ColorIcon(color);
		colorlb = new JLabel(coloricon);
		familypn.add(colorlb);

		// Setting slider
		int size = getFamily().getSize();
		familypn.add(new JLabel("Size:"));
		sizesl = new JSlider(0, 500, size);
		familypn.add(sizesl);
		sizetf = new JTextField(3);
		sizetf.setText(Integer.toString(size));
		familypn.add(sizetf);

		

		// Setting pane listeners

		geditbt.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				new EditFamiliesJDialog(XmippParticlePickerJFrame.this, true);
				
			}
		});

		sizetf.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				int size = Integer.parseInt(sizetf.getText());
				switchSize(size);

			}
		});

		sizesl.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				int size = sizesl.getValue();
				switchSize(size);
			}
		});

		familiescb.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				color = (getFamily().getColor());
				colorlb.setIcon(new ColorIcon(color));
				sizesl.setValue(getFamily().getSize());
			}
		});
	}

	public JMenuBar getPPMenuBar() {
		JMenuBar mb = new JMenuBar();

		// Setting menus
		JMenu filemn = new JMenu("File");
		JMenu filtersmn = new JMenu("Filters");
		JMenu windowmn = new JMenu("Window");
		JMenu helpmn = new JMenu("Help");
		mb.add(filemn);
		mb.add(filtersmn);
		mb.add(windowmn);
		mb.add(helpmn);

		// Setting menu items
		JMenuItem savemi = new JMenuItem("Save");
		filemn.add(savemi);
		JMenuItem stackmi = new JMenuItem("Generate Stack...");
		filemn.add(stackmi);

		JMenuItem bcmi = new JMenuItem("Brightness/Contrast...");
		filtersmn.add(bcmi);
		bcmi.addActionListener(this);

		JMenuItem fftbpf = new JMenuItem("Bandpass Filter...");
		filtersmn.add(fftbpf);
		fftbpf.addActionListener(this);
		JMenuItem admi = new JMenuItem("Anisotropic Diffusion...");
		filtersmn.add(admi);
		admi.addActionListener(this);
		JMenuItem msmi = new JMenuItem("Mean Shift ...");
		filtersmn.add(msmi);
		msmi.addActionListener(this);
		JMenuItem sbmi = new JMenuItem("Substract Background...");
		filtersmn.add(sbmi);
		sbmi.addActionListener(this);
		JMenuItem gbmi = new JMenuItem("Gaussian Blur...");
		filtersmn.add(gbmi);
		gbmi.addActionListener(this);

		JMenuItem particlesmn = new JMenuItem("Particles");
		windowmn.add(particlesmn);
		JMenuItem ijmi = new JMenuItem("ImageJ");
		windowmn.add(ijmi);

		JMenuItem hcontentsmi = new JMenuItem("Help Contents...");
		helpmn.add(hcontentsmi);

		// Setting menu item listeners

		ijmi.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				if(IJ.getInstance() == null)
					new ImageJ();
				IJ.getInstance().setVisible(true);
			}
		});
		savemi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				saveChanges();
				
			}
		});
		stackmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				generateStack();

			}
		});

		
		

		particlesmn.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				List<ImagePlus> imgs = new ArrayList<ImagePlus>();
				for(Particle p: getMicrograph().getParticles())
					imgs.add(p.getImage(getMicrograph().getImage(), getFamily().getSize()));
				String filename = XmippJ.saveTempImageStack(imgs);
				ImagesWindowFactory.openFileAsImage(filename);
			}
		});

		

		return mb;
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		try
		{
			activemacro = ((JMenuItem)e.getSource()).getText();
			IJ.run(activemacro);
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
			JOptionPane.showMessageDialog(this, ex.getMessage());
		}

	}
	
	private void initMicrographsPane()
	{
		micrographpn = new JPanel();
		micrographpn.setBorder(BorderFactory.createTitledBorder("Micrograph"));
		JScrollPane sp = new JScrollPane();
		
		mglist = new JList(ppdata.getMicrographs().toArray());
		mglist.setSelectedIndex(0);
		mglist.addListSelectionListener(new ListSelectionListener() {
			
			@Override
			public void valueChanged(ListSelectionEvent e) {
				//XmippParticlePickerJFrame.this.iw.close();
				initializeCanvas();
				
			}
		});
		sp.setViewportView(mglist);
		micrographpn.add(sp);
		
	}
	
	private void switchSize(int size) {
		sizetf.setText(Integer.toString(size));
		sizesl.setValue(size);
		canvas.repaint();
		getFamily().setSize(size);
	}

	protected void generateStack() {
		// TODO Auto-generated method stub

	}
	
	

	void initializeCanvas() {
		Micrograph micrograph = getMicrograph();
		if(iw == null)
		{
			canvas = new PPCanvas(this, micrograph);
			iw = new ImageWindow(micrograph.getImage(), canvas);
		}
		else
			canvas.setMicrograph(micrograph);
		iw.setTitle(micrograph.getName());
	}
	


	public void saveChanges() {

		ppdata.saveFamilyData();
		ppdata.saveParticles(getMicrograph());
		changed = false;
	}


	public void updateFamilies() {
		Family item = (Family) familiescb.getSelectedItem();
		DefaultComboBoxModel model = new DefaultComboBoxModel(ppdata
				.getFamilies().toArray());
		familiescb.setModel(model);
		familiescb.setSelectedItem(item);
		color = item.getColor();
		colorlb.setIcon(new ColorIcon(color));
		sizesl.setValue(item.getSize());
		pack();
		canvas.repaint();
		changed = true;
	}

	public void addGroup(Family g) {
		if (ppdata.existsFamilyName(g.getName()))
			throw new IllegalArgumentException(
					Constants.getAlreadyExistsGroupNameMsg(g.getName()));
		ppdata.getFamilies().add(g);
		updateFamilies();
	}

	public void removeFamily(Family family) {
		ppdata.getFamilies().remove(family);
		updateFamilies();
	}

	public void setChanged(boolean changed) {
		this.changed = changed;
		
	}

	
}
