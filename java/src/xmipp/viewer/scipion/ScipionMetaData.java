/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.swing.SwingUtilities;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.StopWatch;
import xmipp.viewer.models.ColumnInfo;

/**
 *
 * @author airen
 */
public class ScipionMetaData extends MetaData{
    private String dbfile;
    private ArrayList<ColumnInfo> columns;
    private String self;
    private int cellwidth = 100, cellheight = 100;
    private ArrayList<EMObject> emobjects;
    private static Map<String, Integer> labels = new HashMap<String, Integer>();

    
    public ScipionMetaData(String dbfile)
    {
        
        this.dbfile = dbfile;
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        Connection c = null;
        Statement stmt = null;
        try {
            String name, alias, clsname; 
            int type;
            boolean allowRender;
            ColumnInfo ci;
            int label = 0;
            Class.forName("org.sqlite.JDBC");
            c = DriverManager.getConnection("jdbc:sqlite:" + dbfile);
            stmt = c.createStatement();
            ResultSet rs = stmt.executeQuery( "SELECT * FROM CLASSES;" );
            while ( rs.next() ) {
               name = rs.getString("label_property");
               alias = rs.getString("column_name");
               clsname = rs.getString("class_name");
               if(isPointer(clsname))
               {
                   if(name.equals("self"))
                       self = clsname;
                   continue;
               }
               if (clsname.equals("Integer"))
               {
                        type = MetaData.LABEL_INT;
               } 
               else if (clsname.equals("Float"))
               {
                        type = MetaData.LABEL_DOUBLE;
               }
               else
               {
                        type = MetaData.LABEL_STRING;
               }
               
               labels.put(name, labels.size());
               label = labels.get(name);
               allowRender = isImage(self, name);
               ci = new ColumnInfo(label, name, alias, type, allowRender);
               columns.add(ci);
            }
            
            Object value = null;
           
            EMObject emo;
            int id;
            ci = getColumnInfo("_index");
            rs = stmt.executeQuery( "SELECT * FROM OBJECTS;" );
            while ( rs.next() ) {
                id = rs.getInt("Id");
                emo = new EMObject(id);
                for(ColumnInfo column: columns)
                {
                    name = column.labelName;
                    alias = column.comment;
                    switch(column.type)
                    {
                        
                        case MetaData.LABEL_INT:
                            value = rs.getInt(alias);
                            break;
                        case MetaData.LABEL_DOUBLE:
                            value = rs.getFloat(alias);
                            break;
                        case MetaData.LABEL_STRING:
                            value = rs.getString(alias);
                            if(name.equals("_filename"))
                            {
                                
                                if(column.render && ci != null)
                                    value = rs.getInt(ci.comment) + "@" + value;
                            }
                            break;
                        default:
                            value = rs.getString(alias);
                       
                    }
                    emo.addValue(value);
                    
                }
                emobjects.add(emo);
            }
            rs.close();
            stmt.close();
            c.close();
            
        } 
        catch ( Exception e ) {
            e.printStackTrace();
            
        }
      
    }
    
    public static void main( String[] args )
    {
        final String dbfile = args[1];
        
        SwingUtilities.invokeLater(new Runnable() {

            @Override
            public void run() {
                ScipionMetaData md = new ScipionMetaData(dbfile);
                new ScipionViewerJFrame(md);
            }
        });
        
    }
    
   
    public String getDBFile()
    {
        return dbfile;
    }

    public List<EMObject> getEMObjects() {
        return emobjects;
    }

    public EMObject getEMObject(long id) {
        for(EMObject emo: emobjects)
            if(emo.id == id)
                return emo;
        return null;
    }

    public String getValueFromLabel(int index, int label) {
        EMObject emo = emobjects.get(index);
        return emo.getValue(getColumnIndex(label)).toString();
    }
    
    
  
    public int getColumnIndex(int label)
    {
        int index = 0;
        for(ColumnInfo ci: columns)
        {
            if(ci.label == label)
                return index;
            index ++;
        }
        return -1;
    }

    public static boolean isPointer(String name) {
        if(name.equals("Particle") || name.equals("Coordinate") 
                || name.equals("Micrograph") || name.equals("CTFModel") 
                || name.equals("Acquisition") || name.equals("Class2D"))
            return true;
        return false;
    }

    private ColumnInfo getColumnInfo(String labelName) {
        for(ColumnInfo ci: columns)
            if(ci.labelName.equals(labelName))
                return ci;
        return null;
    }
    public class EMObject
    {
        long id;
        List<Object> values;
       
        boolean isenabled;
        
        public EMObject(long id)
        {
            this.id = id;
            this.isenabled = true;
            values = new ArrayList<Object>();
        }
        
        public void addValue(Object value)
        {
            values.add(value);
        }

        List<Object> getValues() {
            return values;
        }    

        Object getValue(int c) {
            return values.get(c);
        }

       
    }

    public int getCellHeight()
    {
        return cellheight;
    }
    
    public int getCellWidth()
    {
        return cellwidth;
    }
    
    public long[] findObjects()
    {
        long[] ids = new long[emobjects.size()];
        for(int i = 0; i < ids.length; i++)
            ids[i] = emobjects.get(i).id;
        return ids;
            
       
    }
    
    public int[] getActiveLabels()
    {
        int[] labels = new int[columns.size()];
        for(int i = 0; i < labels.length; i ++)
            labels[i] = columns.get(i).label;
        return labels;
    }
    
    public boolean containsLabel(int label)
    {
        for(int i = 0; i < columns.size(); i ++)
            if(columns.get(i).label == label)
                return true;
        return false;
    }
    
    public boolean setEnabled(boolean isenabled, long id)
    {
        EMObject emo = getEMObject(id);
        emo.isenabled = isenabled;
        return true;
    }
    
    public static String getLabelName(int label)
    {
        for(Map.Entry<String, Integer> e: labels.entrySet())
            if( e.getValue().intValue() == label)
                return e.getKey();
        return null;
    }
    
    public ArrayList<ColumnInfo> getColumnsInfo()
    {
        return columns;
    }
    
    public static boolean isImage(String self, String label)
    {
        if(self.equals("Micrograph") || self.equals("Particle"))
        {
            if(label.equals("_filename"))
                return true;
        }
        else if (self.equals("CTFModel"))
        {
            if(label.equals("_micFile"))
                return true;
        }
        else if (self.equals("Class2D"))
        {
            if(label.equals("_representative._filename"))
                return true;
        }
        return false;
    }
    
    public static boolean hasIndex(String self)
    {
        if(self.equals("Particle"))
            return true;
        return false;
    }
    
    public boolean containsGeometryInfo()
    {
        return false;
    }
    
    public boolean isColumnFormat()
    {
        return true;
    }
    
    public Object getValueObject(int label, long id)
    {
        EMObject emo = getEMObject(id);
        int c = getColumnIndex(label);
        if(c != -1)
            return emo.getValue(c);
        return null;
    }
    
    public int getValueInt(int label, long id)
    {
        Object value = getValueObject(label, id);
        return (Integer)value;
    }
    
    public double getValueDouble(int label, long id)
    {
        Object value = getValueObject(label, id);
        return (Float)value;
    }
    
    public String getValueString(int label, long id)
    {
        Object value = getValueObject(label, id);
        if(value == null)
            return "";
        return value.toString();
    }
    
    public String getSelf()
    {
        return self;
    }
    
    public int size()
    {
        return emobjects.size();
    }
    


            
}
