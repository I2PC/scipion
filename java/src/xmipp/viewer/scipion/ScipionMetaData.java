/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Color;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import xmipp.jni.MetaData;
import xmipp.viewer.models.ClassInfo;
import xmipp.viewer.models.ColumnInfo;

/**
 *
 * @author airen
 */
public class ScipionMetaData extends MetaData{
    private String dbfile;
    private List<ColumnInfo> columns;
    private String self;
    private List<EMObject> emobjects;
    private Map<String, Integer> labels = new HashMap<String, Integer>();
    private List<ScipionMetaData> childmds;
    private ScipionMetaData parent;
    private List<ClassInfo> classes;
    private String id;
    
    public ScipionMetaData()
    {
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        
    }
    
    public ScipionMetaData(List<ColumnInfo> columns, List<EMObject> emobjects)
    {
        this.columns = columns;
        this.emobjects = emobjects;
    }
    
    public ScipionMetaData(String dbfile)
    {
        this.dbfile = dbfile;
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        
        loadData("Classes", "Objects");
        if(self.equals("Class2D"))
        {
            childmds = new ArrayList<ScipionMetaData>();
            classes = new ArrayList<ClassInfo>();
            String id;
            String classestb = "%s_Classes";
            String objectstb = "%s_Objects";
            
            for(EMObject emo: emobjects)
            {
                id = String.format("Class00%s", emo.id);
                emo.childmd = new ScipionMetaData(dbfile, String.format(classestb, id), String.format(objectstb, id), id);
               
                emo.childmd.setParent(this);
                childmds.add(emo.childmd);
                classes.add(new ClassInfo(emo.childmd.getBlock(), Color.red, emo.childmd.size()));
            }
            id = "Representatives";
            ScipionMetaData childmd = new ScipionMetaData(dbfile, String.format(classestb, id), String.format(objectstb, id), id);
            childmd.setParent(this);
            childmds.add(childmd);
            
        }
    }
    
    public ScipionMetaData(String dbfile, String classestb, String objectstb, String id)
    {
        this(dbfile, classestb, objectstb);
        this.id = id;
    }
    
    public ScipionMetaData(String dbfile, String classestb, String objectstb)
    {
        
        this.dbfile = dbfile;
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        loadData(classestb, objectstb);
        
    }
    
    public void loadData(String classestb, String objectstb)
    {
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
            ResultSet rs;
            
            String query = String.format( "SELECT * FROM %s;", classestb);
            rs = stmt.executeQuery(query );
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
            int id, index;
            ColumnInfo indexci;
            query =  String.format("SELECT * FROM %s;", objectstb);
            rs = stmt.executeQuery( query );
            while ( rs.next() ) {
                id = rs.getInt("Id");
                emo = new EMObject(id, this);
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
                            if(name.endsWith("_filename"))
                            {
                                indexci = getColumnInfo(name.replace("_filename", "_index"));
                                if(column.render && indexci != null)
                                {
                                    index = rs.getInt(indexci.comment);
                                    if(index > 0)
                                        value = index + "@" + value;
                                }
                            }
                            break;
                        default:
                            value = rs.getString(alias);
                       
                    }
                    //System.out.printf("%s: %s render: %s type: %s \n", column.labelName, value, column.render, MetaData.getLabelTypeString(column.type));
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
    
    public void setParent(ScipionMetaData md)
    {
        this.parent = md;
    }
    
    public String[] getBlocks()
    {
        if(childmds == null)
            return new String[]{getBlock()};
        String[] blocks = new String[childmds.size() + 1];
        blocks[0] = getBlock();
        for(int i = 1; i < blocks.length; i  ++)
            blocks[i] = childmds.get(i - 1).getBlock();
        return blocks;
            
    }
    
    public String getBlock()
    {
        if (id == null)
            return self + "s";
        return id + self + "s";
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
        if(index >= size())
            return null;
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

    public ClassInfo[] getClasses() {
        if(classes == null)
            return null;
        return classes.toArray(new ClassInfo[]{});
    }

    public List<ScipionMetaData> getChilds() {
        return childmds;
    }

    
    public class EMObject
    {
        long id;
        List<Object> values;
       
        boolean isenabled;
        ScipionMetaData childmd;
        ScipionMetaData md;
        
        public EMObject(long id, ScipionMetaData md)
        {
            this.id = id;
            this.isenabled = true;
            values = new ArrayList<Object>();
            this.md = md;
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

        public boolean setValue(int label, Object value) {
            int i = 0;
            for(ColumnInfo ci: columns)
            {
                if(ci.label == label)
                {
                    values.set(i, value);
                    return true;
                }
                i ++;
            }       
            return false;
        }

       
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
    
    public boolean getEnabled(long id)
    {
        EMObject emo = getEMObject(id);
        return emo.isenabled;
    
    }
    
    
    
    public List<ColumnInfo> getColumnsInfo()
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
    
    public boolean setValueString(int label, String value, long id)
    {
        EMObject emo = getEMObject(id);
        return emo.setValue(label, value);
    }

    public  void write(String path)
    {
        //save data on sqlite db path
    }
    
    public  void writeBlock(String path)
    {
        //save block data on sqlite db path
    }
    
    public void writeBlock(String block, String path, long[] ids)
    {
        for(Long id: ids)
            System.out.println(id);
    }
         
    public void importObjects(MetaData md, long [] ids)
    {
        ScipionMetaData smd = (ScipionMetaData)md;
        EMObject emo;
        for(int i = 0; i < ids.length; i++)
        {
            emo = smd.getEMObject(ids[i]);
            emobjects.add(emo);
            System.out.println(emo.id);
        }
    }
    
    public List<EMObject> getChilds(long[] ids)
    {
        ArrayList<EMObject> childs = new ArrayList<EMObject>();
        EMObject emo;
        for(int i = 0; i < ids.length; i++)
        {
            emo = getEMObject(ids[i]);
            if(emo.childmd != null)
                childs.addAll(emo.childmd.getEMObjects());
        }
        return childs;
    }
    
    public void unionAll(MetaData md)
    {
        emobjects.addAll(((ScipionMetaData)md).getEMObjects());
    }
    
    public void destroy()
    {}
    
    public ScipionMetaData getParent() {
        return parent;
    }
    
    public ScipionMetaData getChild(String block)
    {
        if(childmds == null)
            return null;
        for(ScipionMetaData child: childmds)
            if(child.getBlock().equals(block))//from parent to child
                return child;
        return null;
        
    }
}
