/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Color;
import java.io.File;
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
    private String self, selfalias;
    private List<EMObject> emobjects;
    private Map<String, Integer> labels = new HashMap<String, Integer>();
    private ScipionMetaData parent;
    private List<ClassInfo> classes;
    private String id;
    private String classestb, objectstb;
    private ScipionMetaData representativesmd;
    private boolean haschilds;
    
    public ScipionMetaData(String classestb, String objectstb, String self, String selfalias, List<ColumnInfo> columns)
    {
        this(classestb, objectstb, self, selfalias, columns, new ArrayList<EMObject>());
    }
    
    public ScipionMetaData(String classestb, String objectstb, String self, String selfalias, List<ColumnInfo> columns, ScipionMetaData representativesstruct)
    {
        this(classestb, objectstb, self, selfalias, columns, new ArrayList<EMObject>());
        representativesmd = representativesstruct;
    }
    
    public ScipionMetaData(String classestb, String objectstb, String self, String selfalias, List<ColumnInfo> columns, List<EMObject> emobjects)
    {
        this.classestb = classestb;
        this.objectstb = objectstb;
        this.self = self;
        this.selfalias = selfalias;
        this.columns = columns;
        this.emobjects = emobjects;

    }
    
    public ScipionMetaData(String dbfile)
    {
        this.dbfile = dbfile;
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        classestb = "Classes";
        objectstb = "Objects";
        loadData(classestb, objectstb);
        if(self.equals("Class2D"))
        {
            haschilds = true;
            classes = new ArrayList<ClassInfo>();
            String id;
            String classestb = "%s_Classes";
            String objectstb = "%s_Objects";
            id = "Representatives";
            representativesmd = new ScipionMetaData(dbfile, String.format(classestb, id), String.format(objectstb, id));
            int i = 0;
            for(EMObject emo: emobjects)
            {
                id = String.format("Class00%s", emo.id);
                emo.childmd = new ScipionMetaData(dbfile, String.format(classestb, id), String.format(objectstb, id), id);
               
                emo.childmd.setParent(this);
                classes.add(new ClassInfo(emo.childmd.getBlock(), Color.red, emo.childmd.size()));
                emo.representative = representativesmd.getEMObjects().get(i);
                i ++;
            }
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
        this.classestb = classestb;
        this.objectstb = objectstb;
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
               if(name.equals("self"))
               {
                    self = clsname;
                    selfalias = alias;
                   continue;
               }
               type = getTypeMapping(clsname);
               
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
                    emo.setValue(column, value);
                    
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
        if(!haschilds)
            return new String[]{getBlock()};
        List<String> childblocks = new ArrayList<String>();
        for(EMObject emo: emobjects)
            if(emo.childmd != null)
                childblocks.add(emo.childmd.getBlock());
        String[] blocks = new String[childblocks.size() + 1];
        blocks[0] = getBlock();
        for(int i = 0; i < childblocks.size(); i  ++)
            blocks[i + 1] = childblocks.get(i);
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
        return emo.getValue(getColumnInfo(label)).toString();
    }
    
    public ColumnInfo getColumnInfo(int label)
    {
        
        for(ColumnInfo ci: columns)
        {
            if(ci.label == label)
                return ci;
        }
        return null;
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

   
    ScipionMetaData getSelectionMd(long[] selIds) {
        ScipionMetaData selmd;
        if(representativesmd != null)
            selmd = new ScipionMetaData(classestb, objectstb, self, selfalias, columns, representativesmd.getStructure(representativesmd.classestb, representativesmd.objectstb));
        else
            selmd = new ScipionMetaData(classestb, objectstb, self, selfalias, columns);
        for(long id: selIds)
            selmd.add(getEMObject(id));

        return selmd;
    }
    
    public void add(EMObject emo)
    {
        emobjects.add(emo);
        if(emo.childmd != null)
            haschilds = true;
        if(representativesmd != null)
            representativesmd.add(emo.representative);

            
    }

    public ScipionMetaData getStructure(String classestb, String objectstb) {
        if(representativesmd == null)
            return new ScipionMetaData(classestb, objectstb, self, selfalias, columns);
        else
            return new ScipionMetaData(classestb, objectstb, self, selfalias, columns, representativesmd.getStructure(representativesmd.classestb, representativesmd.objectstb));
    }
    
    public class EMObject
    {
        long id;
        Map<ColumnInfo, Object> values;
       
        boolean isenabled;
        ScipionMetaData childmd;
        ScipionMetaData md;
        EMObject representative;
        
        public EMObject(long id, ScipionMetaData md)
        {
            this.id = id;
            this.isenabled = true;
            values = new HashMap<ColumnInfo, Object>();
            this.md = md;
        }
        
        
        

        Map<ColumnInfo, Object> getValues() {
            return values;
        }    

        Object getValue(ColumnInfo c) {
            return values.get(c);
        }

        public boolean setValue(ColumnInfo ci, Object value) {
            values.put(ci, value);
            return true;
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
        if(self.equals("Micrograph") || self.equals("Particle") || self.equals("Volume"))
        {
            if(label.equals("_filename"))
                return true;
        }
        else if (self.equals("CTFModel"))
        {
            if(label.equals("_micFile"))
                return true;
        }
        else if (self.equals("Class2D") || self.equals("Class3D"))
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
        ColumnInfo c = getColumnInfo(label);
        if(c != null)
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
        return setValueObject(label, value, id);
    }
    
    public boolean setValueDouble(int label, double value, long id)
    {
        return setValueObject(label, new Double(value).floatValue(), id);
    }
    
    public boolean setValueInt(int label, int value, long id)
    {
        return setValueObject(label, value, id);
    }
    
     public boolean setValueObject(int label, Object value, long id)
    {
        EMObject emo = getEMObject(id);
        ColumnInfo ci = getColumnInfo(label);
        return emo.setValue(ci, value);
    }

    public void importObjects(MetaData md, long [] ids)
    {
        ScipionMetaData smd = (ScipionMetaData)md;
        EMObject emo;
        for(int i = 0; i < ids.length; i++)
        {
            emo = smd.getEMObject(ids[i]);
            emobjects.add(emo);
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
        if(!haschilds)
            return null;
        for(EMObject emo: emobjects)
            if(emo.childmd != null && emo.childmd.getBlock().equals(block))//from parent to child
                return emo.childmd;
        return null;
        
    }
    
     public  void write(String path)
    {
        new File(path).delete();//overwrite file
        if (parent != null)
        {
            parent.write(path);//this will write parent and siblings
            return;
        }
        writeBlock(path);
        if(haschilds)
            for(EMObject emo: emobjects)
                if(emo.childmd != null)
                    emo.childmd.writeBlock(path);
        if(representativesmd != null)
            representativesmd.writeBlock(path);
        
    }
     
    public  void writeBlock(String path)
    {//overwrites file if exists, does not save recursively metadata childs
        
        Connection c = null;
        Statement stmt = null;
        try {
            
            Class.forName("org.sqlite.JDBC");
            c = DriverManager.getConnection("jdbc:sqlite:" + path);
            String sql = String.format("DROP TABLE IF EXISTS %1$s; CREATE TABLE %1$s(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n" +
"                      label_property      TEXT UNIQUE,\n" +
"                      column_name TEXT UNIQUE,\n" +
"                      class_name TEXT DEFAULT NULL\n" +
"                      )", classestb);
                
            stmt = c.createStatement();
            stmt.executeUpdate(sql);
            //System.out.println(sql);
            sql = String.format("INSERT INTO %s(id, label_property, column_name, class_name) values (1, \'self\', \'%s\', \'%s\')", classestb, selfalias, self);
            String line = ", (%s, \'%s\', \'%s\', \'%s\')";
            ColumnInfo ci;
            String createcols = "", cols = "id, label, comment", type;
            for(int i = 0; i < columns.size(); i ++)
            {
                ci = columns.get(i);
                type = getTypeMapping(ci.type);
                sql += String.format(line, i + 2, ci.labelName, ci.comment, type);
                createcols += String.format(",\n%s %s DEFAULT NULL", ci.comment, type);
                cols += String.format(", %s", ci.comment);
            }
            stmt.executeUpdate(sql);
            //System.out.println(sql);
            sql = String.format("DROP TABLE IF EXISTS %1$s; CREATE TABLE %1$s(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n" +
"                      label      TEXT DEFAULT NULL,\n" +
"                      comment TEXT DEFAULT NULL" +
"                      %2$s)", objectstb, createcols);
                
            
            stmt.executeUpdate(sql);
            //System.out.println(sql);
            sql = String.format("INSERT INTO %s(%s) VALUES ", objectstb, cols);
            Object value;
            for(EMObject emo: emobjects)
            {
                sql += String.format("(%s, '', ''", emo.id);
                for(int i = 0; i < columns.size(); i ++)
                {
                    ci = columns.get(i);
                    value = emo.values.get(i);
                    if(ci.type == MetaData.LABEL_STRING)
                    {
                        if(value != null)
                        {
                            String str = (String)value;
                            if( str.contains("@"))
                                str = str.substring(str.lastIndexOf("@") + 1);
                            value = str;
                        }
                        sql += String.format(", '%s'", value);
                    }
                    else
                        sql += String.format(", %s", value);
                }
                sql += "),";
            }
            sql = sql.substring(0, sql.length() - 1);//remove first comma
            stmt.executeUpdate(sql);
            //System.out.println(sql);
            stmt.close();
            c.close();
            
        } 
        catch ( Exception e ) {
            e.printStackTrace();
            
        }
    }
    
    public static String getTypeMapping(int type)
    {
         if (type == MetaData.LABEL_INT)
            return "Integer";
        if (type == MetaData.LABEL_DOUBLE)
            return "Float";
        return "String";
    }
    
    public static int getTypeMapping(String type)
    {
        int result;
         if (type.equals("Integer"))
            return MetaData.LABEL_INT;
        if (type.equals("Float"))
            return MetaData.LABEL_DOUBLE;
        return MetaData.LABEL_STRING;
    }
    
    public String getColumnsTable()
    {
        return classestb;
    }
    
    public String getObjectsTable()
    {
        return objectstb;
    }
}
