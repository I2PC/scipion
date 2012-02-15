package xmipp.viewer;

import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

import xmipp.utils.DEBUG;

public class TestRenderer extends DefaultTableCellRenderer {
	private static final long serialVersionUID = 1L;

   public TestRenderer() {
       super();
       setOpaque(true);
       setHorizontalAlignment(JLabel.LEFT);
       setHorizontalTextPosition(JLabel.CENTER);
       setVerticalTextPosition(JLabel.BOTTOM);
   }

   @Override
   public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
       Object item = (Object) object;
       // Calls super class so foreground, background, borders and rest of stuff is set.
       super.getTableCellRendererComponent(table, null,
               item != null && isSelected,
               item != null && hasFocus, row, column);

       if (item != null) {
    	   //DEBUG.printMessage("*** 2. Rendering: " + item.getLabel());
           //AbstractXmippTableModel tableModel = (AbstractXmippTableModel) table.getModel();
    	   String v = item.toString();
    	   try {
    		   float f = Float.parseFloat(v);
    		   v = String.format("%10.4f", f);
    	   } catch (Exception e) {
			// TODO: handle exception
		}
           setText(v);
       }
       return this;
   }
}
