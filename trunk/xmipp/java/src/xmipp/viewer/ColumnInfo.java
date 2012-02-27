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
	public boolean allowEdit = true;
	
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
	
	/** Update the column information with the provided one
	 * return true if some field has changed
	 */
	public boolean updateInfo(ColumnInfo ci){
		boolean result = visible != ci.visible 
				|| render != ci.render || allowEdit != ci.allowEdit;
		visible = ci.visible;
		render = ci.render;
		allowEdit = ci.allowEdit;
		return result;
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
	
}//class ColumnInfo
