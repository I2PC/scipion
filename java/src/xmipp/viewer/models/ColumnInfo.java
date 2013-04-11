/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.viewer.models;

import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

import java.util.ArrayList;

/** This class will store info about how to display label on gallery */
public class ColumnInfo {
	protected int label;
	public String labelName; //This is redundant, but for avoid a function call
	public boolean visible; 
	public boolean render;
	public boolean allowRender = false;
	public boolean allowEdit = true;
	public String comment;
	
	/** Constructors */
	public ColumnInfo(int label, String name, boolean visible, boolean render){
		this.label = label;
		this.labelName = name;
		this.visible = visible;
		this.render = render;
		try {
			this.allowRender = MetaData.isImage(label);
			comment = MetaData.getLabelComment(label);
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
		boolean result = labelName != ci.labelName || visible != ci.visible 
				|| render != ci.render || allowEdit != ci.allowEdit;
		visible = ci.visible;
		render = ci.render;
		allowEdit = ci.allowEdit;
		labelName = ci.labelName;
		return result;
	}
	
	public int getLabel(){
		return label;
	}
	
	public int getType(){
		try {
			return MetaData.getLabelType(label);
		}catch(Exception e){
			e.printStackTrace();
		}
		return -1;
	}
	
	public String getLabelTypeString(){
		try {
			int type = MetaData.getLabelType(label);
			return MetaData.getLabelTypeString(type);
		}catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	/** Return the name of the label */
	public String getLabelName(){
		return labelName;
	}
	
	/** Return the real name of the label
	 * This is read from metadata 
	 * */
	public String getLabelRealName(){
		return MetaData.getLabelName(label);
	}
	
	/** Change display name of the label */
	public void changeLabelName(String newName){
		labelName = newName;
	}
	
	public ColumnInfo clone(){
		return new ColumnInfo(label, labelName, visible, render);
	}
	
	@Override
	public boolean equals(Object o){
		ColumnInfo col = (ColumnInfo)o;
		return label == col.label && 
				visible == col.visible && 
				render == col.render;
	}
	
	public class ColumnExtraInfo {
		public String color;
		public String marker;
		public String linestyle;
		
		public ColumnExtraInfo(String c){
			color = c;
			linestyle = "solid";
			marker = "none";
		}
	}
	
	public String toString()
	{
		return getLabelName();
	}
	
}//class ColumnInfo
