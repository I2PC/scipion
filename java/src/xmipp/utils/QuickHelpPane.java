package xmipp.utils;



import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;




public class QuickHelpPane extends JPanel
{
	private JPanel buttonspn;
	private String title;
	private JTable helptb;
	private Map<String, String> helpmap;
	private String[] keys;
	
	public QuickHelpPane(String title, Map<String, String> helpmap)
	{
		
		this.helpmap = helpmap;
		if(helpmap.size() == 0)
			throw new IllegalArgumentException("There is no help information available");
		keys = helpmap.keySet().toArray(new String[]{});
		this.title = title;
		initComponents();
	}

	private void initComponents()
	{
		setBorder(BorderFactory.createTitledBorder(title));
		setLayout(new GridBagLayout());
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.insets = new Insets(10, 10, 10, 10);
		setLayout(new GridBagLayout());

		JScrollPane sp = new JScrollPane();
		add(sp, XmippWindowUtil.getConstraints(constraints, 0, 0, 3));
		initHelpTable();
		sp.setViewportView(helptb);
		int height = Math.min(600, helptb.getRowHeight() * helpmap.size() + 5);
		sp.setPreferredSize(new Dimension(500, height));
		constraints.gridy = 1;
		buttonspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));

		add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1));
		XmippWindowUtil.setLocation(0.5f, 0.5f, this);
		setVisible(true);
	}

	private void initHelpTable()
	{
		helptb = new JTable();
		helptb.setShowHorizontalLines(false);
		helptb.setTableHeader(null);
		int lines = 2;
	    helptb.setRowHeight( helptb.getRowHeight() * lines);
		helptb.setDefaultRenderer(String.class, new MultilineCellRenderer());
		helptb.setModel(new AbstractTableModel()
		{
			
			@Override
			public Object getValueAt(int row, int column)
			{
				String key = keys[row];
				if(column == 0)
					return key;
				return helpmap.get(key);
						
			}
			
			@Override
			public Class getColumnClass(int column)
			{
				return String.class;
			}
			
			@Override
			public int getRowCount()
			{
				return helpmap.size();
			}
			
			@Override
			public int getColumnCount()
			{
				return 2;
			}
		});
		
	}
	
	void addButton(JButton bt)
	{
		buttonspn.add(bt);
	}


}
