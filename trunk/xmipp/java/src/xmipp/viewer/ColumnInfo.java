package xmipp.viewer;

import xmipp.jni.MetaData;
import java.util.ArrayList;

/** This class will store info about how to display label on gallery */
public class ColumnInfo {
	protected int label;
	protected String labelName; //This is redundant, but for avoid a function call
	public boolean visible; 
	public boolean render;
	public boolean allowRender = false;
	
	/** Constructors */
	public ColumnInfo(int label, String name, boolean visible, boolean render){
		this.label = label;
		this.labelName = name;
		this.visible = visible;
		this.render = render;
		try {
			this.allowRender = MetaData.isImage(label);
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public ColumnInfo(int label){
		this(label, MetaData.getLabelName(label), true, false);
	}
	
	public ColumnInfo(int label, boolean visible){
		this(label, MetaData.getLabelName(label), visible, false);
	}
	
	public ColumnInfo(int label, boolean visible, boolean render){
		this(label, MetaData.getLabelName(label), visible, render);
	}	
	
	public int getLabel(){
		return label;
	}
	
	public String getLabelName(){
		return labelName;
	}
	
	public ColumnInfo clone(){
		return new ColumnInfo(label, labelName, visible, render);
	}
	
	@Override
	public boolean equals(Object o){
		ColumnInfo col = (ColumnInfo)o;
		return visible == col.visible && render == col.render;
	}
	
	/** Create an array list of column info given the md */
	public static ArrayList<ColumnInfo> createListFromMd(MetaData md){
		ArrayList<ColumnInfo> labels = null;
		try {
		int [] lab = md.getActiveLabels();
		labels = new ArrayList<ColumnInfo>(lab.length);
		for (int i = 0; i < lab.length; ++i) 
			labels.add(new ColumnInfo(lab[i]));
		} catch (Exception e){
			e.printStackTrace();
		}
		return labels; 
	}//function createListFromMd
	
}
