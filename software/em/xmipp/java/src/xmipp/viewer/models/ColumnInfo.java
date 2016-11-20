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

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

/** This class will store info about how to display label on gallery */
public class ColumnInfo {
    

	public int label;
	public String labelName; //This is redundant, but for avoid a function call
	public boolean visible; 
	public boolean render;
	public boolean allowRender = false;
	public boolean allowEdit = false;
	public String comment;
        public int type;
	
	/** Constructors */
	
	public ColumnInfo(int label, String name, String comment, int type, boolean allowRender, boolean visible, boolean render, boolean allowEdit){
		this.label = label;
		this.labelName = name;
		this.visible = visible;
		this.render = render;
		this.allowRender = allowRender;
                this.comment = comment;
                this.type = type;
                this.allowEdit = allowEdit;
	}
        
	
        public ColumnInfo(int label, String name, String comment, int type, boolean allowRender,boolean allowEdit){
            this(label, name, comment, type, allowRender, true, allowRender, allowEdit);
        }
	
	public ColumnInfo(int label){
		this(label, MetaData.getLabelName(label), MetaData.getLabelComment(label), MetaData.getLabelType(label), MetaData.isImage(label), true, false, true);
	}
        
	
	public ColumnInfo(int label, boolean visible){
		this(label, MetaData.getLabelName(label), MetaData.getLabelComment(label), MetaData.getLabelType(label), MetaData.isImage(label), visible, false, true);
	}
	
	public ColumnInfo(int label, boolean visible, boolean render){
		this(label, MetaData.getLabelName(label), MetaData.getLabelComment(label), MetaData.getLabelType(label), MetaData.isImage(label), visible, render, true);
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
	
	public String getLabelTypeString(){
		try {
			
			return MetaData.getLabelTypeString(type);
		}catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}

	
	public ColumnInfo clone(){
		return new ColumnInfo(label, labelName, comment, type, allowRender, visible, render, allowEdit);
	}
	
	@Override
	public boolean equals(Object o){

        if (o == null) return false;
        if (!getClass().equals(o.getClass())) return false;

        ColumnInfo col = (ColumnInfo)o;
		return label == col.label && 
				visible == col.visible && 
				render == col.render;
	}

    public boolean isEnable() {
        return label == MDLabel.MDL_ENABLED || labelName.equals("enabled");
    }
	
	public class ColumnExtraInfo {
		public String color;
		public String marker;
		public String linestyle;
                public boolean plot;
		
		public ColumnExtraInfo(String c){
			color = c;
			linestyle = "solid";
			marker = "none";
                        plot = false;
		}
	}
	
	public String toString()
	{
		return labelName;
	}

        
}//class ColumnInfo
