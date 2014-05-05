/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.util.List;
import javax.swing.table.AbstractTableModel;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.scipion.ScipionMetaData.EMObject;

/**
 *
 * @author airen
 */
public class ScipionViewerTableModel extends AbstractTableModel{
    private ScipionMetaData md;
    private List<ScipionMetaData.EMObject> emobjects;
    private List<ColumnInfo> columns;
    
    public ScipionViewerTableModel(ScipionMetaData md)
    {
        this.md = md;
        this.emobjects = md.getEMObjects();
        this.columns = md.getColumnsInfo();

    }
    
    @Override
    public String getColumnName(int c)
    {
        
        return columns.get(c).labelName;
    }
    
    @Override
    public Class getColumnClass(int c)
    {
        if(getRowCount() == 0)
            return String.class;
        Object value = emobjects.get(0).getValue(c);
        if(value == null)
            return String.class;
        
        return value.getClass();
    }

    @Override
    public int getRowCount() {
        
        return emobjects.size();
        
    }

    @Override
    public int getColumnCount() {
        return columns.size();
    }

    @Override
    public Object getValueAt(int row, int col) {
        EMObject emo = emobjects.get(row);
        return emo.getValue(col);
        
    }
    
}
