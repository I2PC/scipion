package xmipp.utils;



import java.awt.Dimension;
import java.util.Map;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;




public class QuickHelpPane extends JPanel
{
	protected JPanel buttonspn;
	protected String title;
	protected JTable helptb;
	protected Map<Object, Object> helpmap;
	protected Object[] keys;
	protected final boolean editmap;

    public QuickHelpPane(String title, Map<Object, Object> helpmap)
    {
        this(title, helpmap, false);
    }
	
	public QuickHelpPane(String title, Map<Object, Object> helpmap, boolean editmap)
	{
		
		this.helpmap = helpmap;
        this.editmap = editmap;
		if(helpmap.size() == 0)
			throw new IllegalArgumentException("There is no help information available");
		keys = helpmap.keySet().toArray();
		this.title = title;
		initComponents();
	}

	protected void initComponents()
	{
//		setBorder(BorderFactory.createTitledBorder(title));

		JScrollPane sp = new JScrollPane();
		add(sp);
		sp.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		sp.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		initHelpTable();
		sp.setViewportView(helptb);
		int height = Math.min(600, helptb.getRowHeight() * helpmap.size() + 5);
		sp.setPreferredSize(new Dimension(800, height));
		setVisible(true);
	}

	protected void initHelpTable()
	{
		helptb = new JTable();
		helptb.setShowHorizontalLines(false);
		helptb.setTableHeader(null);
//		int lines = 2;
//        helptb.setRowHeight( helptb.getRowHeight() * lines);
		helptb.setModel(new AbstractTableModel()
		{
			
			@Override
			public Object getValueAt(int row, int column)
			{
				Object key = keys[row];
				if(column == 0)
					return key;
				return helpmap.get(key);
						
			}
                        
            @Override
			public boolean isCellEditable(int row, int column)
			{
				if(column == 0)
                    return false;
                return editmap;
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
                        
                        @Override
			public void setValueAt(Object value, int row, int column)
			{
				helpmap.put(keys[row], value);
			}
                        
                        
		});
		helptb.setDefaultRenderer(String.class, new MultilineCellRenderer());
		helptb.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		helptb.getColumnModel().getColumn(0).setPreferredWidth(250);
		helptb.getColumnModel().getColumn(1).setPreferredWidth(550);
		
	}
	
	void addButton(JButton bt)
	{
		buttonspn.add(bt);
	}


}
