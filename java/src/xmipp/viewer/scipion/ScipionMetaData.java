/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import xmipp.jni.MetaData;
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
    private ScipionMetaData parent;
    private String id;
    private String classestb, objectstb;
    private boolean haschilds;
    private static int labelscount = 0;
    private ColumnInfo enabledci;
    
    public ScipionMetaData(String classestb, String objectstb, String self, String selfalias, List<ColumnInfo> columns)
    {
        this(classestb, objectstb, self, selfalias, columns, new ArrayList<EMObject>());
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
            String id;
            haschilds = true;
            String classestb = "%s_Classes";
            String objectstb = "%s_Objects";
            
            int i = 0;
            for(EMObject emo: emobjects)
            {
                id = String.format("Class%03d", emo.id);
                emo.childmd = new ScipionMetaData(dbfile, String.format(classestb, id), String.format(objectstb, id), id);
                emo.childmd.setParent(this);
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
            name = alias = "enabled";
            labelscount ++;
            enabledci = new ColumnInfo(labelscount, name, alias, MetaData.LABEL_INT, false);
            columns.add(enabledci);
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
               
               labelscount ++;
               allowRender = isImage(self, name);
               ci = new ColumnInfo(labelscount, name, alias, type, allowRender);
               columns.add(ci);
            }
            
            Object value = null;
           
            EMObject emo;
            int id, index;
            String label, comment;
            boolean enabled;
            ColumnInfo indexci;
            query =  String.format("SELECT * FROM %s;", objectstb);
            rs = stmt.executeQuery( query );
            while ( rs.next() ) {
                id = rs.getInt("id");
                label = rs.getString("label");
                comment = rs.getString("comment");
                enabled = (rs.getInt("enabled") == 1)? true: false;
                emo = new EMObject(id, label, comment, enabled, this);
                for(ColumnInfo column: columns)
                {
                    if(!column.isEnable())
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
                }
                emo.setValue(enabledci, 1);
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
    
    public class EMObject
    {
        long id;
        String label, comment;
        protected Map<ColumnInfo, Object> values;
        boolean changed;
        ScipionMetaData childmd;
        ScipionMetaData md;
        
        public EMObject(long id, String label, String comment, boolean enabled, ScipionMetaData md)
        {
            this.id = id;
            this.label = label;
            this.comment = comment;
            values = new HashMap<ColumnInfo, Object>();
            setEnabled(enabled);
            this.md = md;
            changed = false;
        }
        
        Object getValue(ColumnInfo c) {
            Object result = values.get(c);
            if(result != null)
                return result;
            for(ColumnInfo ci: columns)
                if(ci.labelName.equals(c.labelName))
                    return values.get(ci);//when mergin metadatas required
            return null;
        }

        public void setValue(ColumnInfo ci, Object value) {
                if(values.containsKey(ci))
                    changed = true;
                values.put(ci, value);
                
        }
        
        public boolean isEnabled()
        {
            Integer value = (Integer)getValue(enabledci);
            return value.intValue() == 1;
        }
        
        public void setEnabled(boolean isenabled)
        {
            int value = isenabled? 1: 0;
            setValue(enabledci, value);
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
        Object value = emo.getValue(getColumnInfo(label));
        if(value == null)
            return null;
        return value.toString();
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
   
    public ScipionMetaData getMd(List<EMObject> emos) {
        ScipionMetaData selmd;
        //either selection of classes or particles main selection will be saved on Classes and Objects table
        selmd = new ScipionMetaData("Classes", "Objects", self, selfalias, columns);
        for(EMObject emo: emos)
            selmd.add(emo);
        
        return selmd;
    }
    
    public void add(EMObject emo)
    {
        emobjects.add(emo);
        if(emo.childmd != null)
            haschilds = true;
    }

    public ScipionMetaData getStructure(String classestb, String objectstb) {
        return new ScipionMetaData(classestb, objectstb, self, selfalias, columns);
    }

    public List<EMObject> getEnabledObjects() {
        List<EMObject> emos = new ArrayList<EMObject>();
        for(EMObject emo: emobjects)
            if(emo.isEnabled())
                emos.add(emo);
        return emos;
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
        emo.setEnabled(isenabled);
        return true;
    }
    
    public boolean getEnabled(long id)
    {
        EMObject emo = getEMObject(id);
        return emo.isEnabled();
    
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
            if(label.equals("_micFile") || label.equals("_psdFile"))
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
        emo.setValue(ci, value);
        return true;
    }

    public void importObjects(MetaData md, long [] ids)
    {
        ScipionMetaData smd = (ScipionMetaData)md;
        EMObject emo;
        for(int i = 0; i < ids.length; i++)
        {
            emo = smd.getEMObject(ids[i]);
            add(emo);
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
        System.out.println("writing file in " + path);
        new File(path).delete();//overwrite file
        if(path.contains("@"))
        {
            writeBlock(path.substring(path.lastIndexOf("@") + 1));
            return;
        }
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
        
    }
     
    public String toString()
    {
            String sql = String.format("DROP TABLE IF EXISTS %1$s; CREATE TABLE %1$s(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n" +
"                      label_property      TEXT UNIQUE,\n" +
"                      column_name TEXT UNIQUE,\n" +
"                      class_name TEXT DEFAULT NULL\n" +
"                      )", classestb);
                
            sql += String.format(";INSERT INTO %s(id, label_property, column_name, class_name) values (1, \'self\', \'%s\', \'%s\')", classestb, selfalias, self);
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
            sql += String.format(";DROP TABLE IF EXISTS %1$s; CREATE TABLE %1$s(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n" +
"                      label      TEXT DEFAULT NULL,\n" +
"                      comment TEXT DEFAULT NULL" +
"                      %2$s)", objectstb, createcols);
            if(size() == 0)    
                return sql;
            sql += String.format(";INSERT INTO %s(%s) VALUES ", objectstb, cols);
            Object value;
            for(EMObject emo: emobjects)
            {
                sql += String.format("(%s, '', ''", emo.id);
                for (ColumnInfo column : columns) {
                        value = emo.getValue(column);
                        if(value != null)
                        {
                            if(column.type == MetaData.LABEL_STRING)
                            {

                                    String str = (String)value;
                                    if( str.contains("@"))
                                        str = str.substring(str.lastIndexOf("@") + 1);
                                    value = str;
                                    sql += String.format(", '%s'", value);


                            }
                            else
                                sql += String.format(", %s", value);
                        }
                        else
                            sql += ", NULL";
                }
                sql += "),";
            }
            sql = sql.substring(0, sql.length() - 1);//remove first comma
            System.out.println(sql);
            return sql;
    }
    public  void writeBlock(String path)
    {//overwrites file if exists, does not save recursively metadata childs
        
        Connection c = null;
        Statement stmt = null;
        try {
            
            Class.forName("org.sqlite.JDBC");
            c = DriverManager.getConnection("jdbc:sqlite:" + path);
            stmt = c.createStatement();
               String sql = String.format("DROP TABLE IF EXISTS %1$s; CREATE TABLE %1$s(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n" +
"                      label_property      TEXT UNIQUE,\n" +
"                      column_name TEXT UNIQUE,\n" +
"                      class_name TEXT DEFAULT NULL\n" +
"                      )", classestb);
            stmt.executeUpdate(sql);
            
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
            
            sql = String.format("DROP TABLE IF EXISTS %1$s; CREATE TABLE %1$s(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n" +
"                      label      TEXT DEFAULT NULL,\n" +
"                      comment TEXT DEFAULT NULL" +
"                      %2$s)", objectstb, createcols);
            stmt.executeUpdate(sql);
            if(size() == 0)    
                return;
            sql = String.format("INSERT INTO %s(%s) VALUES ", objectstb, cols);
            Object value;
            int count = 0;
            for(EMObject emo: emobjects)
            {
                sql += String.format("(%s, '', ''", emo.id);
                for (ColumnInfo column : columns) {
                        value = emo.getValue(column);
                        if(value != null)
                        {
                            if(column.type == MetaData.LABEL_STRING)
                            {

                                    String str = (String)value;
                                    if( str.contains("@"))
                                        str = str.substring(str.lastIndexOf("@") + 1);
                                    value = str;
                                    sql += String.format(", '%s'", value);


                            }
                            else
                                sql += String.format(", %s", value);
                        }
                        else
                            sql += ", NULL";
                }
                sql += "),";
                count ++;
                if(count == 200)
                {
                     
                    sql = sql.substring(0, sql.length() - 1);//remove comma
                    stmt.executeUpdate(sql);
                    sql = String.format("INSERT INTO %s(%s) VALUES ", objectstb, cols);
                    count = 0;
                }
            }
            if(count > 0)
            {
                sql = sql.substring(0, sql.length() - 1);//remove comma
                stmt.executeUpdate(sql);
            }
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
    
    public void print()
    {
        System.out.println(toString());
    }
    
    public double[] getStatistics(boolean applyGeo)
    {
        return null;
    }
    

    public void sort(int sortLabel, boolean ascending)
    {
        ColumnInfo ci = getColumnInfo(sortLabel);
        Comparable valuei, valuej;
        EMObject emo;
        boolean bigger, exchange;
        for(int i = 0; i < emobjects.size() - 1; i ++)
        {
            valuei = (Comparable)emobjects.get(i).getValue(ci);
            for(int j = i + 1; j < emobjects.size(); j ++)
            {
                valuej = (Comparable)emobjects.get(j).getValue(ci);
                bigger = valuei.compareTo(valuej) > 0;
                exchange = ascending? bigger: !bigger;
                if(exchange) 
                {
                    emo = emobjects.get(i);
                    emobjects.set(i, emobjects.get(j));
                    emobjects.set(j, emo);
                    valuei = valuej;
                }
            }
        }   
    }
    
    public  void overwrite(String src, String path)
    {
        try {
            Files.copy(Paths.get(src), Paths.get(path), StandardCopyOption.REPLACE_EXISTING);
            if (parent != null)
            {
                parent.overwrite(src, path);//this will write parent and siblings
                return;
            }
            overwriteBlock(path);
            if(haschilds)
                for(EMObject emo: emobjects)
                    if(emo.childmd != null)
                        emo.childmd.overwriteBlock(path);
        } catch (IOException ex) {
            Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    public  void overwriteBlock(String path)
    {//overwrites enabled column in existing sqlite objects table
        Connection c = null;
        Statement stmt = null;
        try {
            
            Class.forName("org.sqlite.JDBC");
            c = DriverManager.getConnection("jdbc:sqlite:" + path);
            stmt = c.createStatement();
            String sql = String.format("ALTER TABLE %s ADD COLUMN _enabled INT NOT NULL DEFAULT 1", objectstb);
            stmt.executeUpdate(sql);
            String updatesql = "UPDATE %s SET %s=%s WHERE id=%s;";
            sql = "";
            int count = 0;
            for(EMObject emo: emobjects)
            {
                if(emo.changed)
                {
                    sql+= String.format(updatesql, objectstb, enabledci.labelName, (emo.isEnabled())? 1: 0, emo.id);
                    count ++;
                    if(count == 200)
                    {
                        stmt.executeUpdate(sql);
                        sql = "";
                        count = 0;
                    }
                }
                
            }
            if(count > 0)
                stmt.executeUpdate(sql);
            stmt.close();
            c.close();
            
        } 
        catch ( Exception e ) {
            e.printStackTrace();
            
        }
    }
    
   
   
}
