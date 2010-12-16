/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

import java.awt.Adjustable;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.text.DecimalFormat;

import javax.swing.BoundedRangeModel;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.event.EventListenerList;

/**
 * @author jcuenca
 * Horizontal scrollbar with range (min / max) labels, and a current index label (with an optional descriptive prefix)
 * Color is used in labels to highlight the current index
 * extends JPanel
 * implements AdjustmentListener (scrollbar change event management), Adjustable (required by AdjustmentEvent at fireAdjustmentValueChanged() )
 */
public class LabelScrollbar extends JPanel implements AdjustmentListener, Adjustable{

	// extent is required by JScrollBar
	private static int extent=20;
	
	private JLabel textLabel,minLabel,maxLabel;
	private String text=""; //optional text to clarify the value shown in textLabel
	
	// the current value is stored in the scrollbar
	private JScrollBar scrollbar;
	private GridBagLayout gbl;
	
	private static String DEFAULT_TEXT_COLOR="GRAY";
	private static String HIGHLIGHTED_TEXT_COLOR="BLACK";
	
	public LabelScrollbar(int min,int max) {
		// JScrollbar range is minimum..maximum+extent (instead of a simple minimum..maximum)
		scrollbar= new JScrollBar(JScrollBar.HORIZONTAL);
		scrollbar.setVisibleAmount(extent);
		scrollbar.addAdjustmentListener(this);
		
		textLabel=new JLabel("text");
		
		// setup min / max labels
		JPanel minMaxPanel = new JPanel();
		minMaxPanel.setLayout(new BoxLayout(minMaxPanel, BoxLayout.LINE_AXIS));
		minLabel=new JLabel("min");
		minLabel.setHorizontalAlignment(JLabel.LEFT);
		minLabel.setAlignmentX(LEFT_ALIGNMENT);
		maxLabel=new JLabel("max");
		maxLabel.setHorizontalAlignment(JLabel.RIGHT);
		maxLabel.setAlignmentX(RIGHT_ALIGNMENT);
		minMaxPanel.add(minLabel);
		minMaxPanel.add(maxLabel);
		
		gbl= new GridBagLayout();
		setLayout(gbl);
		
		// add each element in a different row
		add(textLabel,1);
		add(scrollbar,2);
		add(minMaxPanel,3);
		
		setMinimum(min);
		setMaximum(max);
		
	}
	
	// row: 1..N
	private void add(JComponent c,int row){
		GridBagConstraints constraints=new GridBagConstraints();
		constraints.fill = GridBagConstraints.BOTH;
		constraints.gridx=0; constraints.gridy = row-1;
		constraints.fill = GridBagConstraints.HORIZONTAL;	
		gbl.setConstraints(c, constraints);
		add(c);
	}
	
	/**
	 * @param n 
	 */
	public void setMaximum(int n){
		// adding the extent is a "requirement" of JScrollbar
		if(scrollbar != null)
			scrollbar.setMaximum(n + extent);
		if(maxLabel != null)
			maxLabel.setText(htmlColor(String.valueOf(n), DEFAULT_TEXT_COLOR, true));
	}
	
	/**
	 * @param n 
	 */
	public void setMinimum(int n){
		if(scrollbar != null)
			scrollbar.setMinimum(n);
		if(minLabel != null)
			minLabel.setText(htmlColor(String.valueOf(n), DEFAULT_TEXT_COLOR, true));
	}

	public void addAdjustmentListener(AdjustmentListener l){
		listenerList.add(AdjustmentListener.class, l);
	}
	
    /**
     * Removes an AdjustmentEvent listener.
     *
     * @param l the AdjustmentLister to remove
     * @see #addAdjustmentListener
     */
    public void removeAdjustmentListener(AdjustmentListener l)  {
        listenerList.remove(AdjustmentListener.class, l);
    }
	
    /**
     * (from JScrollBar)
     * Notify listeners that the scrollbar's model has changed.
     * 
     * @see #addAdjustmentListener
     * @see EventListenerList
     */
    private void fireAdjustmentValueChanged(int id, int type, int value,
					    boolean isAdjusting) {
        Object[] listeners = listenerList.getListenerList();
        AdjustmentEvent e = null;
        for (int i = listeners.length - 2; i >= 0; i -= 2) {
            if (listeners[i]==AdjustmentListener.class) {
                if (e == null) {
                    e = new AdjustmentEvent(this, id, type, value, isAdjusting);
                }
                ((AdjustmentListener)listeners[i+1]).adjustmentValueChanged(e);
            }          
        }
    }   
	
	public void adjustmentValueChanged(AdjustmentEvent e){
		if(e.getSource()==scrollbar){
			updateTextLabel(); // projection number
			
			// notify listeners
			int id = AdjustmentEvent.ADJUSTMENT_VALUE_CHANGED;
			int type = AdjustmentEvent.TRACK;
			BoundedRangeModel model = scrollbar.getModel();
			int value = model.getValue();
			boolean isAdjusting = model.getValueIsAdjusting();
			fireAdjustmentValueChanged(id, type, value, isAdjusting);
		}
		// notify();
	}
	
	public void setValue(int v){
		scrollbar.setValue(v);
	}
	
	public int getValue(){
		return scrollbar.getValue();
	}
	
	/**
	 * 
	 * @return value as string with filling zeroes (so the component does not need to resize when changing tenths or hundreds)
	 */
	public String getStringValue(){
		DecimalFormat df=new DecimalFormat();
		df.setMinimumIntegerDigits(String.valueOf(getMaximum()).length());
		return df.format(getValue());
	}
	
	private void updateTextLabel(){
		textLabel.setText("<html>" + htmlColor(getText(), DEFAULT_TEXT_COLOR, false) + htmlColor(getStringValue(), HIGHLIGHTED_TEXT_COLOR, false)+ "</html>");
	}
	 

	/**
	 * @return the text
	 */
	private String getText() {
		return text;
	}

	/**
	 * @param text the text to set
	 */
	public void setText(String text) {
		this.text = text;
		updateTextLabel();
	}
	
	private String htmlColor(String text, String color, boolean htmlHeader){
		String res= "<font color=" + color + ">"+ text + "</font>";
		if(htmlHeader)
			res="<html>" + res + "</html>";
		return res;
	}
	
	/******************** Methods required by Adjustable interface *********************************/
	
    /**
     * Returns the minimum value supported by the scrollbar 
     * (usually zero).
     *
     * @return the value of the model's minimum property
     * @see #setMinimum
     */
    public int getMinimum() { 
        return scrollbar.getMinimum(); 
    }
    
    /**
     * The maximum value of the scrollbar is maximum - extent.
     *
     * @return the value of the model's maximum property
     * @see #setMaximum
     */
    public int getMaximum() { 
        return scrollbar.getMaximum(); 
    }

    /**
     * For backwards compatibility with java.awt.Scrollbar.
     * @see Adjustable#getUnitIncrement
     * @see #getUnitIncrement(int)
     */
    public int getUnitIncrement() {
        return scrollbar.getUnitIncrement();
    }
    
    /**
     * Sets the unitIncrement property.
     * <p>
     * Note, that if the argument is equal to the value of Integer.MIN_VALUE,
     * the most look and feels will not provide the scrolling to the right/down.
     * @see #getUnitIncrement
     * @beaninfo
     *   preferred: true
     *       bound: true
     * description: The scrollbar's unit increment.
     */
    public void setUnitIncrement(int unitIncrement) { 
    	scrollbar.setUnitIncrement(unitIncrement);
    }
    
    /**
     * For backwards compatibility with java.awt.Scrollbar.
     * @see Adjustable#getBlockIncrement
     * @see #getBlockIncrement(int)
     */
    public int getBlockIncrement() {
        return scrollbar.getBlockIncrement();
    }

    /**
     * Sets the blockIncrement property.
     * <p>
     * Note, that if the argument is equal to the value of Integer.MIN_VALUE,
     * the most look and feels will not provide the scrolling to the right/down.
     * @see #getBlockIncrement()
     * @beaninfo
     *   preferred: true
     *       bound: true
     * description: The scrollbar's block increment.
     */
    public void setBlockIncrement(int blockIncrement) { 
    	scrollbar.setBlockIncrement(blockIncrement);
    }
    
    /**
     * Returns the scrollbar's extent, aka its "visibleAmount".  In many 
     * scrollbar look and feel implementations the size of the 
     * scrollbar "knob" or "thumb" is proportional to the extent.
     * 
     * @return the value of the model's extent property
     * @see #setVisibleAmount
     */
    public int getVisibleAmount() { 
        return scrollbar.getVisibleAmount();
    }

    
    /**
     * Set the model's extent property.
     * 
     * @see #getVisibleAmount
     * @see BoundedRangeModel#setExtent
     * @beaninfo
     *   preferred: true
     * description: The amount of the view that is currently visible.
     */
    public void setVisibleAmount(int extent) { 
    	scrollbar.setVisibleAmount(extent);
    }
    
    /**
     * Returns the component's orientation (horizontal or vertical). 
     *                     
     * @return VERTICAL or HORIZONTAL
     * @see #setOrientation
     * @see java.awt.Adjustable#getOrientation
     */
    public int getOrientation() { 
        return scrollbar.getOrientation();
    }
    
    /********************** Main - component testing **************************/

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	    JFrame frame = new JFrame("Test Frame");
	    frame.setSize(400, 400);
	    frame.setLocationRelativeTo(null);
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	    LabelScrollbar ls = new LabelScrollbar(1,10);  // Initialize the component.
	    ls.setText("Projection #");
	    frame.getContentPane().add(ls);      // Place the component on the application
	                                             // window such that it fills the whole
	                                             // window frame.
	    frame.setVisible(true);

	}

	
}
