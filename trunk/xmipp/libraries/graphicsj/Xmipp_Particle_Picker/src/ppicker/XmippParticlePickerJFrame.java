package ppicker;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

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
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.ListCellRenderer;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

enum Tool {
	IMAGEJ, PICKER
}

enum Shape {
	CIRCLE, RECTANGLE, BOTH
}

public class XmippParticlePickerJFrame extends JFrame implements ActionListener {

	private JButton infobt;
	private JRadioButton pickrbt;
	private JSlider radiussl;
	private JRadioButton circlerbt;
	private JRadioButton rectanglerbt;
	private ButtonGroup shapebg;
	private PPCanvas canvas;
	private JTextField radiustf;
	private JRadioButton moverbt;
	private JRadioButton pickactionrbt;
	private ButtonGroup actionbg;
	private Tool tool = Tool.PICKER;
	private ImagePlus img;
	private JButton adjustbt;
	private JRadioButton bothrbt;
	private Shape shape;
	private JMenuBar mb;
	private JMenu filemn;
	private JMenuItem loadmi;
	private JMenuItem savemi;
	private JMenuItem stackmi;
	private JComboBox familiescb;
	private JButton geditbt;
	private PPData ppdata;
	private JLabel colorlb;
	private ColorIcon coloricon;
	private Color color;
	private JMenu helpmn;
	private JMenuItem hcontentsmi;

	public XmippParticlePickerJFrame(ImagePlus img) {
		
		this.img = img;
		ppdata = new PPData();

		initComponents();
		initializeImage(img);
	}

	public XmippParticlePickerJFrame() {
		img = WindowManager.getCurrentImage();
		ppdata = new PPData();
		initComponents();
		initializeImage(img);
	}

	private void initComponents() {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Particle Picker");

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(0, 5, 0, 5);
		constraints.anchor = GridBagConstraints.WEST;

		mb = new JMenuBar();
		filemn = new JMenu();
		filemn.setText("File");
		mb.add(filemn);
		setJMenuBar(mb);
		loadmi = new JMenuItem("Load Particles ...");
		filemn.add(loadmi);
		savemi = new JMenuItem("Save Particles ...");
		filemn.add(savemi);
		stackmi = new JMenuItem("Generate Stack ...");
		filemn.add(stackmi);
		helpmn = new JMenu();
		helpmn.setText("Help");
		hcontentsmi = new JMenuItem("Help Contents ...");
		helpmn.add(hcontentsmi);
		mb.add(helpmn);
		setLayout(new GridBagLayout());

		JPanel pickerpn = new JPanel(new GridBagLayout());
		pickerpn.setBorder(BorderFactory.createTitledBorder("Picker"));

		JPanel grouppn = new JPanel();
		pickerpn.add(new JLabel("Family:"), WindowUtils.updateConstraints(constraints, 0, 0, 1));
		familiescb = new JComboBox(ppdata.getFamilies().toArray());
		grouppn.add(familiescb);
		color = Family.getDefaultFamily().getColor();
		grouppn.add(new JLabel("Color:"));
		coloricon = new ColorIcon(color);
		colorlb = new JLabel(coloricon);
		grouppn.add(colorlb);
		
		int radius = getFamily().getRadius();
		grouppn.add(new JLabel("Radius:"));
		radiussl = new JSlider(0, 500, radius);
		grouppn.add(radiussl);
		radiustf = new JTextField(3);
		radiustf.setText(Integer.toString(radius));
		grouppn.add(radiustf);
		add(pickerpn, WindowUtils.updateConstraints(constraints, 0, 0, 3));
		
		geditbt = new JButton("Edit");
		grouppn.add(geditbt);
		
		pickerpn.add(grouppn, WindowUtils.updateConstraints(constraints, 1, 0, 1));
		
		JPanel symbolpn = new JPanel();
		shapebg = new ButtonGroup();
		circlerbt = new JRadioButton("Circle");
		shapebg.add(circlerbt);
		rectanglerbt = new JRadioButton("Rectangle");
		shapebg.add(rectanglerbt);
		bothrbt = new JRadioButton("Both");
		bothrbt.getModel().setSelected(true);
		shape = shape.BOTH;
		shapebg.add(bothrbt);
		pickerpn.add(new JLabel("Symbol:"), WindowUtils.updateConstraints(constraints, 0, 1, 1));
		symbolpn.add(bothrbt);
		symbolpn.add(circlerbt);
		symbolpn.add(rectanglerbt);
		pickerpn.add(symbolpn, WindowUtils.updateConstraints(constraints, 1, 1, 1));
		
		
		geditbt.addActionListener(new PPEditGroupsActionListener());

		pickactionrbt = new JRadioButton("Picker");
		pickactionrbt.setSelected(true);

		actionbg = new ButtonGroup();
		actionbg.add(pickactionrbt);
		actionbg.add(moverbt);

		JPanel ctfpn = new JPanel();
		ctfpn.setBorder(BorderFactory.createTitledBorder("Micrograph"));
		add(ctfpn, WindowUtils.updateConstraints(constraints, 0, 2, 3));

		
		adjustbt = new JButton("Adjust Image");
		JPanel buttons = new JPanel();
		infobt = new JButton("View Info");
		infobt.addActionListener(this);
		buttons.add(infobt);
		buttons.add(adjustbt);
		add(buttons, WindowUtils.updateConstraints(constraints, 0, 3, 3));
		setListeners();
		pack();
		centerScreen();
		setVisible(true);
	}

	private void setListeners() {

		familiescb.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				color = (getFamily().getColor());
				colorlb.setIcon(new ColorIcon(color));
				radiussl.setValue(getFamily().getRadius());
			}
		});

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent winEvt) {
				System.exit(0);
			}

		});

		stackmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				generateStack();

			}
		});

		circlerbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				canvas.repaint();
				shape = Shape.CIRCLE;
			}
		});

		loadmi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				loadParticles();
			}
		});

		savemi.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				saveParticles();

			}
		});


		rectanglerbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				canvas.repaint();
				shape = Shape.RECTANGLE;
			}
		});

		adjustbt.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				IJ.run("Brightness/Contrast...");

			}
		});

		radiustf.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				int radius = Integer.parseInt(radiustf.getText());
				switchRadius(radius);
				
			}
		});

		radiussl.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				int radius = radiussl.getValue();
				switchRadius(radius);
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

	class PPEditGroupsActionListener implements ActionListener {

		@Override
		public void actionPerformed(ActionEvent e) {
			new EditFamiliesJDialog(XmippParticlePickerJFrame.this, true);

		}
	}
	
	private void switchRadius(int radius)
	{
		radiustf.setText(Integer.toString(radius));
		radiussl.setValue(radius);
		canvas.repaint();
		getFamily().setRadius(radius);
	}

	protected void generateStack() {
		// TODO Auto-generated method stub

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		new FamilyParticlesJDialog(this, getFamily().getParticles());

	}

	void initializeImage(ImagePlus image) {
		this.img = image;
		canvas = new PPCanvas(image, this);
		new ImageWindow(img, canvas);
	}

	// centers the dialog within the screen
	public void centerScreen() {
		Dimension dim = getToolkit().getScreenSize();
		Rectangle abounds = getBounds();
		int x = 9 * (dim.width - abounds.width) / 10;
		int y = (dim.height - abounds.height) / 2;
		setLocation(x, y);
	}

	public boolean onPick() {
		return pickrbt.getModel().isSelected();
	}

	public Shape getShape() {
		return shape;
	}

	public Tool getTool() {
		return tool;
	}

	public ImagePlus getImage() {
		return img;
	}

	public void saveParticles() {
		JFileChooser fc = new JFileChooser();
		int returnVal = fc.showSaveDialog(this);

		try {
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				String filename = fc.getSelectedFile().getAbsolutePath();
				ppdata.saveData(filename);
			}
		} catch (Exception e) {
			JOptionPane.showMessageDialog(this, e.getMessage());
		}

	}

	public void loadParticles() {
		JFileChooser fc = new JFileChooser();
		int returnVal = fc.showOpenDialog(this);

		try {
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				String filename = fc.getSelectedFile().getAbsolutePath();
				ppdata.loadData(filename);
				canvas.repaint();
				updateFamilies();
			}
		} catch (Exception e) {
			JOptionPane.showMessageDialog(this, e.getMessage());
		}
	}

	public PPData getPPData() {
		return ppdata;
	}


	public Family getFamily() {
		return (Family) familiescb.getSelectedItem();
	}

	public void updateFamilies() {
		Family item = (Family)familiescb.getSelectedItem();
		DefaultComboBoxModel model = new DefaultComboBoxModel(ppdata
				.getFamilies().toArray());
		familiescb.setModel(model);
		familiescb.setSelectedItem(item);
		color = item.getColor();
		colorlb.setIcon(new ColorIcon(color));
		radiussl.setValue(item.getRadius());
		pack();
		canvas.repaint();
	}

	public void addGroup(Family g) {
		if (ppdata.existsGroupName(g.getName()))
			throw new IllegalArgumentException(Constants.getAlreadyExistsGroupNameMsg(g.getName()));
		ppdata.getFamilies().add(g);
		updateFamilies();
	}

}
