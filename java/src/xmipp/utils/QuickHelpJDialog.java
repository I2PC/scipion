package xmipp.utils;


import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Map;

import javax.swing.JButton;
import javax.swing.JDialog;



public class QuickHelpJDialog extends JDialog
{

		private JButton okbt;
		private QuickHelpPane helppn;
		private String helptitle;
		private Map<String, String> helpmap;

		public QuickHelpJDialog(Frame window, boolean modal, String helptitle, Map<String, String> map)
		{
			super(window, modal);
			this.helpmap = map;
			this.helptitle = helptitle;
			initComponents();
		}

		private void initComponents()
		{
			setResizable(false);
			setDefaultCloseOperation(HIDE_ON_CLOSE);
			setTitle(helptitle);
			GridBagConstraints constraints = new GridBagConstraints();
			constraints.insets = new Insets(10, 10, 10, 10);
			setLayout(new GridBagLayout());
			helppn = new QuickHelpPane( helptitle, helpmap);
			add(helppn, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
			okbt = XmippWindowUtil.getTextButton("Ok", new ActionListener()
			{

				@Override
				public void actionPerformed(ActionEvent e)
				{
					setVisible(false);

				}
			});
			helppn.addButton(okbt);
			pack();
			XmippWindowUtil.setLocation(0.5f, 0.5f, this);
		}

		

}
