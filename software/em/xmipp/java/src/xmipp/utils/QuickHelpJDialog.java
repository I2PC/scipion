package xmipp.utils;


import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Map;

import javax.swing.JButton;
import javax.swing.JDialog;



public class QuickHelpJDialog extends JDialog
{

		protected JButton okbt;
		protected QuickHelpPane helppn;
		protected String helptitle;
		protected Map<Object, Object> helpmap;
		protected final boolean editmap;
        
        public QuickHelpJDialog(Frame window, boolean modal, String helptitle, Map<Object, Object> map)
        {
            this(window, modal, helptitle, map, false);
        }

		public QuickHelpJDialog(Frame window, boolean modal, String helptitle, Map<Object, Object> map, boolean editmap)
		{
			super(window, modal);
			this.helpmap = map;
			this.helptitle = helptitle;
            this.editmap = editmap;
			initComponents();
		}

		private void initComponents()
		{
			setIconImage(XmippResource.getIcon("xmipp_logo.png").getImage());
			setResizable(false);
			setDefaultCloseOperation(HIDE_ON_CLOSE);
			setTitle(helptitle);
			GridBagConstraints constraints = new GridBagConstraints();
			constraints.insets = new Insets(5, 5, 5, 5);
			setLayout(new GridBagLayout());
			helppn = new QuickHelpPane( helptitle, helpmap, editmap);
			add(helppn, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));
			okbt = XmippWindowUtil.getTextButton("Ok", new ActionListener()
			{

				@Override
				public void actionPerformed(ActionEvent e)
				{
					setVisible(false);

				}
			});
			add(okbt, XmippWindowUtil.getConstraints(constraints, 0, 1));
			pack();
			XmippWindowUtil.setLocation(0.5f, 0.5f, this);
		}

		

}
