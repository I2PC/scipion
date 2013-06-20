package xmipp.viewer.particlepicker;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JLabel;
import javax.swing.JPanel;


import xmipp.jni.MetaData;
import xmipp.utils.ColorIcon;
import xmipp.utils.XmippWindowUtil;

public class ColorHelper
{
	int id;
	String name;
	private MetaData md;
	private double max;
	private double min;
	private double range;

	private final static int LOW = 0;
	private final static int HIGH = 255;
	private final static int HALF = (HIGH + 1) / 2;
 
	private final static Map<Integer, Color> map = initNumberToColorMap();
	private static int factor;
	
	
	
	public ColorHelper(int id, String name, Color color, MetaData md)
	{
		this.id = id;
		this.name = name;
		this.md = md;
		max = Math.rint(md.getColumnMax(id) * 100)/100;
		min = Math.rint(md.getColumnMin(id) * 100)/100;
		range = max - min;
	}
	
	public ColorHelper(int id, String name, MetaData md)
	{
		this(id, name, getNextColor(), md);
	}


	private static Color[] colors = new Color[] { Color.BLUE, Color.CYAN, Color.GREEN, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.YELLOW };
	private static int nextcolor;

	public static Color getNextColor()
	{
		Color next = colors[nextcolor];
		nextcolor++;
		if (nextcolor == colors.length)
			nextcolor = 0;
		return next;
	}
	
	public static Color[] getSampleColors() {
		return colors;
	}
	
	public String toString()
	{
		return name;
	}

	

	public String getName()
	{
		return name;
	}

	public int getId()
	{
		return id;
	}
	
	public Color getColor(double score)
	{
		
		double percent = 1;
		if(range > 0)
			percent = (score - min)/range;
		Color scorecolor = numberToColorPercentage(percent);
		return scorecolor;
	}
	
	
 
	/**
	 * 
	 * @param value
	 *            should be from 0 unti 100
	 */
	public static Color numberToColor(final double value) {
		if (value < 0 || value > 100) {
			return null;
		}
		return numberToColorPercentage(value / 100);
	}
 
	/**
	 * @param value
	 *            should be from 0 unti 1
	 * @return
	 */
	public static Color numberToColorPercentage(final double value) {
		if (value < 0 || value > 1) {
			return null;
		}
		Double d = value * factor;
		int index = d.intValue();
		if (index == factor) {
			index--;
		}
		return map.get(index);
	}
	
	public static JPanel getColorMap()
	{
		JPanel mappn = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
		for(int i = 0; i < map.size(); i += 10)
			mappn.add(new JLabel(new ColorIcon(map.get(i), 1, 20, 0, false, false)));
		return mappn;
	}
	
	
	/**
	 * @return
	 */
	private static Map<Integer, Color> initNumberToColorMap() {
		HashMap<Integer, Color> localMap = new HashMap<Integer, Color>();
		int r = LOW;
		int g = LOW;
		int b = HALF;
 
		// factor (increment or decrement)
		int rF = 0;
		int gF = 0;
		int bF = 1;
 
		int count = 0;
		// 1276 steps
		while (true) {
			localMap.put(count++, new Color(r, g, b));
			if (b == HIGH) {
				gF = 1; // increment green
			}
			if (g == HIGH) {
				bF = -1; // decrement blue
				// rF = +1; // increment red
			}
			if (b == LOW) {
				rF = +1; // increment red
			}
			if (r == HIGH) {
				gF = -1; // decrement green
			}
			if (g == LOW && b == LOW) {
				rF = -1; // decrement red
			}
			if (r < HALF && g == LOW && b == LOW) {
				break; // finish
			}
			r += rF;
			g += gF;
			b += bF;
			r = rangeCheck(r);
			g = rangeCheck(g);
			b = rangeCheck(b);
		}
		initList(localMap);
		return localMap;
	}
 
	/**
	 * @param localMap
	 */
	private static void initList(final HashMap<Integer, Color> localMap) {
		List<Integer> list = new ArrayList<Integer>(localMap.keySet());
		Collections.sort(list);
		Integer min = list.get(0);
		Integer max = list.get(list.size() - 1);
		factor = max + 1;
	}
 
	/**
	 * @param value
	 * @return
	 */
	private static int rangeCheck(final int value) {
		if (value > HIGH) {
			return HIGH;
		} else if (value < LOW) {
			return LOW;
		}
		return value;
	}

	public double getMin()
	{
		return min;
	}
	
	public double getMax()
	{
		return max;
	}
 
	

}
