/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import ij.IJ;
import ij.ImagePlus;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.sqlite.SQLiteConfig;
import xmipp.ij.commons.ImagePlusFromFile;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.CTFDescription;
import xmipp.jni.EllipseCTF;
import xmipp.jni.Filename;
import xmipp.jni.MetaData;
import xmipp.utils.StopWatch;
import xmipp.viewer.models.ColumnInfo;

/**
 *
 * @author airen
 */
public class ScipionMetaData extends MetaData {

    protected List<ColumnInfo> columns;
    protected List<EMObject> emobjects;
    protected ScipionMetaData parent;
    protected boolean haschilds;
    protected static int labelscount = 0;
    protected ColumnInfo idci, labelci, commentci, enabledci;
    protected String preffix = "";
    protected String[] blocks;
    protected int enableds;
    protected Boolean checkTmp;
    protected HashMap<String, String> properties;
    protected HashMap <Long, EMObject> idsmap;
    protected String self;
    

    public ScipionMetaData(String dbfile) {
       
        this.filename = dbfile;
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        loadData();
        if (isClassificationMd()) {
            String preffix;
            haschilds = true;
            List<String> childblocks = new ArrayList<String>();
            for (EMObject emo : emobjects) {
                if((Integer)emo.getValue("_size") != 0)
                {
                    preffix = String.format("Class%03d_", emo.getId());
                    emo.childmd = new ScipionMetaData(dbfile, preffix);
                    emo.childmd.setParent(this);
                    childblocks.add(emo.childmd.getBlock());
                }
            }
            blocks = new String[childblocks.size() + 1];
            blocks[0] = getBlock();
            for (int i = 0; i < childblocks.size(); i++) {
                blocks[i + 1] = childblocks.get(i);

            }
        }
        else
            blocks = new String[]{getBlock()};

    }

    public ScipionMetaData(String dbfile, String preffix) {

        this.filename = dbfile;
        columns = new ArrayList<ColumnInfo>();
        emobjects = new ArrayList<EMObject>();
        this.preffix = preffix;
        loadData();
        blocks = new String[]{getBlock()};

    }

    public void loadData() {
        Connection c = null;
        Statement stmt = null;
        try {
        	String query;
//        	 c = etConnection(filename);
             c = getConnectionReadOnly(filename);
             stmt = c.createStatement();
             ResultSet rs;
        	if (preffix == null || preffix.equals(""))
            {
        		query = "SELECT name FROM sqlite_master WHERE type='table' AND name='Properties';";
                rs = stmt.executeQuery(query);
        		boolean exists = rs.next();
        		if(exists)
        		{
	                properties = new HashMap<String, String>();
	                String key, value;
	                query = "SELECT * FROM Properties;";
	                rs = stmt.executeQuery(query);
	
	                while (rs.next()) {
	                    key = rs.getString("key");
	                    value = rs.getString("value");
	                    properties.put(key, value);
	                }
        		}
            }
            String name, alias, clsname;
            int type;
            boolean allowRender;
            ColumnInfo ci;

            name = alias = "enabled";
            labelscount++;
            enabledci = new ColumnInfo(labelscount, name, alias, MetaData.LABEL_INT, false, true);
            columns.add(enabledci);

            name = alias = "id";
            labelscount++;
            idci = new ColumnInfo(labelscount, name, alias, MetaData.LABEL_INT, false, false);
            columns.add(idci);

            name = alias = "label";
            labelscount++;
            labelci = new ColumnInfo(labelscount, name, alias, MetaData.LABEL_STRING, false, false);
            columns.add(labelci);

            name = alias = "comment";
            labelscount++;
            commentci = new ColumnInfo(labelscount, name, alias, MetaData.LABEL_STRING, false, false);
            columns.add(commentci);

           

            loadValues(c);
           
            
            query = String.format("SELECT * FROM %sClasses;", preffix);
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                name = rs.getString("label_property");
                alias = rs.getString("column_name");
                clsname = rs.getString("class_name");
                if(name.equals("self"))
                {
                	self = clsname;
                	continue;
                }
                
                type = getTypeMapping(clsname);
                if(type == -1)
                	continue;//complex object is empty, unless is matrix
                labelscount++;
                allowRender = isImage(name);
                ci = new ColumnInfo(labelscount, name, alias, type, allowRender, false);
                columns.add(ci);
            }
            

            rs.close();
            stmt.close();
            c.close();
            StopWatch.getInstance().printElapsedTime("loaded data");
        } catch (Exception e) {
            throw new IllegalArgumentException(e);
            

        }
    }

    protected void loadValues(Connection c)
    {
        try {
            idsmap = new HashMap<Long, EMObject>();
            Statement stmt = c.createStatement();
            ResultSet rs;
            EMObject emo;
            String alias;
            
            String query = String.format("SELECT id, label, comment, enabled FROM %sObjects;", preffix);
            rs = stmt.executeQuery(query);
            Object value;
            int index = 0;
            while (rs.next()) {
                emo = new EMObject(index, this);
                for (ColumnInfo column : columns) {
                    alias = column.comment;
                    switch (column.type) {
                        
                        case MetaData.LABEL_INT:
                            value = rs.getInt(alias);
                            break;
                        case MetaData.LABEL_DOUBLE:
                            value = rs.getFloat(alias);
                            break;
                        default:
                            value = rs.getString(alias);
                    }
                    emo.values.put(column, value);
                }
                emobjects.add(emo);
                idsmap.put(emo.getId(), emo);
                if(emo.isEnabled())
                    enableds ++;
                else
                    emo.changed = true;
                index ++;
            }
        } catch (SQLException ex) {
            throw new IllegalArgumentException(ex);
            
        }
    }
    
    public void setParent(ScipionMetaData md) {
        this.parent = md;
    }

    public String[] getBlocks() {
        return blocks;

    }
    

    public String getBlock() {
    	String setType = getSetType();
    	if(setType.isEmpty())
    		setType = "noname";
    	String block = setType.replace("SetOf", "");
    	
        if (preffix == null) {
            return block;
        }
        return preffix + block;
    }

   

    public List<EMObject> getEMObjects() {
        return emobjects;
    }

    public EMObject getEMObject(long id) {
        return idsmap.get(id);
    }

    public synchronized String getValueFromLabel(int index, int label) {
        if (index >= size()) {
            return null;
        }
        EMObject emo = emobjects.get(index);
        Object value = emo.getValue(getColumnInfo(label));
        if (value == null) {
            return null;
        }
        return value.toString();
    }

    public ColumnInfo getColumnInfo(int label) {

        for (ColumnInfo ci : columns) {
            if (ci.label == label) {
                return ci;
            }
        }
        return null;
    }

    public ColumnInfo getColumnInfo(String labelName) {
        for (ColumnInfo ci : columns) {
            if (ci.labelName.equals(labelName)) {
                return ci;
            }
        }
        return null;
    }

   
    public void add(EMObject emo) {
        emobjects.add(emo);
        if (emo.childmd != null) {
            haschilds = true;
        }
    }

   
    

    public long[] findObjects() {
        long[] ids = new long[emobjects.size()];
        for (int i = 0; i < ids.length; i++) {
            ids[i] = emobjects.get(i).getId();
        }
        return ids;

    }

    public int[] getActiveLabels() {
        int[] labels = new int[columns.size()];
        for (int i = 0; i < labels.length; i++) {
            labels[i] = columns.get(i).label;
        }
        return labels;
    }

    @Override
    public boolean containsLabel(int label) {
        for (int i = 0; i < columns.size(); i++) {
            if (columns.get(i).label == label) {
                return true;
            }
        }
        return false;
    }

    @Override
    public boolean setEnabled(boolean isenabled, long id) {
        EMObject emo = getEMObject(id);
        emo.setEnabled(isenabled);
        return true;
    }

    @Override
    public boolean getEnabled(long id) {
        EMObject emo = getEMObject(id);
        return emo.isEnabled();

    }

    public List<ColumnInfo> getColumnsInfo() {
        return columns;
    }

    public boolean isImage(String label) {
    	String setType = getSetType();
    	if (setType == null)
    		return false;
        if(label.contains("_filename"))
        {
            String indexLabel = label.replace("_filename", "_index");
            if(getColumnInfo(indexLabel) != null)
                return true;
        }
        if (setType.equals("SetOfMicrographs") || setType.equals("SetOfParticles") || setType.equals("SetOfVolumes")) {
            if (label.equals("_filename")) {
                return true;
            }
        } else if (setType.equals("SetOfCTF")) {
            if (label.equals("_micObj._filename") || label.equals("_psdFile") || label.equals("_xmipp_enhanced_psd") 
                    ||label.equals("_xmipp_ctfmodel_quadrant") ||label.equals("_xmipp_ctfmodel_halfplane")) {
                return true;
            }
        } else if (setType.equals("SetOfClasses2D") || setType.equals("SetOfClasses3D")) {
            if (label.equals("_representative._filename")) {
                return true;
            }
        }
        return false;
    }

    public static boolean hasIndex(String self) {
        if (self.equals("Particle")) {
            return true;
        }
        return false;
    }

    public boolean isColumnFormat() {
        return true;
    }

    public synchronized Object getValueObject(int label, long id) {
        EMObject emo = getEMObject(id);
        ColumnInfo c = getColumnInfo(label);
        if (c != null) {
            return emo.getValue(c);
        }
        return null;
    }

    public int getValueInt(int label, long id) {
        Object value = getValueObject(label, id);
        if (value != null) {
            return (Integer) value;
        }
        return Integer.MIN_VALUE;
    }

    public double getValueDouble(int label, long id) {
        Object value = getValueObject(label, id);
        if (value != null) {
            return (Float) value;
        }
        return Double.NaN;
    }

    public String getValueString(int label, long id) {
        Object value = getValueObject(label, id);
        if (value == null) 
            return "";
        if (value instanceof Float)
            return String.format("%.2f", value);
        return value.toString();
    }

    
    @Override
    public int size() {
        return emobjects.size();
    }

    @Override
    public boolean setValueString(int label, String value, long id) {
        return setValueObject(label, value, id);
    }

    @Override
    public boolean setValueDouble(int label, double value, long id) {
        return setValueObject(label, new Double(value).floatValue(), id);
    }

    @Override
    public boolean setValueInt(int label, int value, long id) {
        return setValueObject(label, value, id);
    }

    public boolean setValueObject(int label, Object value, long id) {
        EMObject emo = getEMObject(id);
        ColumnInfo ci = getColumnInfo(label);
        emo.setValue(ci, value);
        return true;
    }

    @Override
    public void importObjects(MetaData md, long[] ids) {
        ScipionMetaData smd = (ScipionMetaData) md;
        EMObject emo;
        for (int i = 0; i < ids.length; i++) {
            emo = smd.getEMObject(ids[i]);
            add(emo);
        }
    }

    public void unionAll(MetaData md) {
        emobjects.addAll(((ScipionMetaData) md).getEMObjects());
    }

    public void destroy() {
    }

    public ScipionMetaData getParent() {
        return parent;
    }

    public ScipionMetaData getChild(String block) {
        if (!haschilds) {
            return null;
        }
        for (EMObject emo : emobjects) {
            if (emo.childmd != null && emo.childmd.getBlock().equals(block))//from parent to child
            {
                return emo.childmd;
            }
        }
        return null;
    }

    public void write(String path) {
        new File(path).delete();//overwrite file
        if (path.contains("@")) {
            writeBlock(path.substring(path.lastIndexOf("@") + 1));
            return;
        }
        if (parent != null) {
            parent.write(path);//this will write parent and siblings
            return;
        }
        writeBlock(path);
        if (haschilds) {
            for (EMObject emo : emobjects) {
                if (emo.childmd != null) {
                    emo.childmd.writeBlock(path);
                }
            }
        }
    }

   

    public void writeBlock(String path) {
        //Might fail if some fixed column was added

        Connection c = null;
        Statement stmt = null;
        try {

            c = getConnection(path);
            c.setAutoCommit(false);
            stmt = c.createStatement();

            String sql = String.format("DROP TABLE IF EXISTS %1$sClasses; CREATE TABLE %1$sClasses(\n"
                    + "id        INTEGER PRIMARY KEY AUTOINCREMENT,\n"
                    + "                      label_property      TEXT UNIQUE,\n"
                    + "                      column_name TEXT UNIQUE,\n"
                    + "                      class_name TEXT DEFAULT NULL\n"
                    + "                      )", preffix);
            stmt.executeUpdate(sql);

            String line = ", (%s, \'%s\', \'%s\', \'%s\')";
            ColumnInfo ci;
            String createcols = "", cols = "", type;
            for (int i = 0; i < columns.size(); i++) {
                ci = columns.get(i);
                type = getTypeMapping(ci.type);
                sql += String.format(line, i + 2, ci.labelName, ci.comment, type);
                createcols += String.format("%s %s DEFAULT NULL,", ci.comment, type);
                cols += String.format("%s,", ci.comment);
            }
            stmt.executeUpdate(sql);
            createcols = createcols.substring(0, createcols.length() - 1);
            cols = cols.substring(0, cols.length() - 1);
            sql = String.format("DROP TABLE IF EXISTS %1$sObjects; CREATE TABLE %1$sObjects(\n"
                    + "                      %2$s)", preffix, createcols);

            stmt.executeUpdate(sql);
            if (size() == 0) {
                return;
            }
            sql = String.format("INSERT INTO %sObjects(%s) VALUES ", preffix, cols);
            Object value;
            for (EMObject emo : emobjects) {
                sql += "(";
                for (ColumnInfo column : columns) {
                    value = emo.getValue(column);
                    if (value != null) {
                        if (column.type == MetaData.LABEL_STRING) {

                            String str = (String) value;
                            if (str.contains("@")) {
                                str = str.substring(str.lastIndexOf("@") + 1);
                            }
                            value = str;
                            sql += String.format("'%s',", value);
                        } else {
                            sql += String.format("%s,", value);
                        }
                    } else {
                        sql += "NULL,";
                    }
                }
                sql = sql.substring(0, sql.length() - 1) + "),";
            }
            sql = sql.substring(0, sql.length() - 1);
            stmt.executeUpdate(sql);
            c.commit();
            stmt.close();
            c.close();

        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public static String getTypeMapping(int type) {
        if (type == MetaData.LABEL_INT) {
            return "Integer";
        }
        if (type == MetaData.LABEL_DOUBLE) {
            return "Float";
        }
        return "String";
    }

    public static int getTypeMapping(String type) {
        if (type.equals("Integer")) {
            return MetaData.LABEL_INT;
        }
        if (type.equals("Float")) {
            return MetaData.LABEL_DOUBLE;
        }
        if(type.equals("Boolean"))
            return MetaData.LABEL_INT;
        
        if(type.equals("String") || type.equals("Matrix") || type.equals("CsvList"))
        	return MetaData.LABEL_STRING;
        return -1;
    }

    public String getPreffix() {
        return preffix;
    }

    public void print() {
        System.out.println(toString());
    }

    public double[] getStatistics(boolean applyGeo) {
        return null;
    }

    public void sort(int sortLabel, boolean ascending) {
        ColumnInfo ci = getColumnInfo(sortLabel);
        try {
            loadNeighborhoodValues(0, ci, size());
        } catch (SQLException ex) {
            Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
        }
        Comparator<EMObject> comparator = new EMObjectComparator(ci);
        if(!ascending)
            comparator = Collections.reverseOrder(comparator);
        Collections.sort(emobjects, comparator);
        for (int i = 0; i < emobjects.size(); i++) 
            emobjects.get(i).setIndex(i);

    }

    public void overwrite(String src, String path, boolean[] selection) throws SQLException {
        try {
            XmippUtil.copyFile(src, path);
            // Do not overwrite the parent. This looses the selection. issue #498
            //            if (parent != null) {
            //                parent.overwrite(src, path, null);//this will write parent and siblings
            //                return;
            //            }
            overwriteBlock(path, selection);
            if (haschilds) {
                for (EMObject emo : emobjects) {
                    if (emo.childmd != null) {
                        emo.childmd.overwriteBlock(path, null);
                    }
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void overwriteBlock(String path, boolean[] selection) throws SQLException {//overwrites enabled column in existing sqlite objects table
        Connection c = null;
        PreparedStatement stmt = null;
        try {

//            Class.forName("org.sqlite.JDBC");
//            c = DriverManager.getConnection("jdbc:sqlite:" + path);
            c = getConnection(path);

            c.setAutoCommit(false);
            stmt = c.prepareStatement(String.format("UPDATE %sObjects SET %s=?, %s=?, %s=? WHERE id=?;", preffix, enabledci.labelName, labelci.labelName, commentci.labelName));
            boolean enabled;
            for (EMObject emo : emobjects) 
            {
                if(selection != null)
                {
                	enabled = selection[emo.index];
                    stmt.setInt(1, enabled ? 1 : 0);
                    stmt.setString(2, emo.getLabel());
                    stmt.setString(3, emo.getComment());
                    stmt.setInt(4, emo.getId().intValue());
                    stmt.executeUpdate();
                 
                }
                else if (emo.changed) {
                    stmt.setInt(1, emo.isEnabled() ? 1 : 0);
                    stmt.setString(2, emo.getLabel());
                    stmt.setString(3, emo.getComment());
                    stmt.setInt(4, emo.getId().intValue());
                    stmt.executeUpdate();
                    
                }
            }
            c.commit();
            
        } 
        catch (Exception e) {
            e.printStackTrace();
            if (c != null) 
                try {
                    System.err.print("Transaction is being rolled back");
                    c.rollback();
                } 
                catch(SQLException excep) {
                    excep.printStackTrace();
                }

        }
        finally {
            if (stmt != null)
                stmt.close();
            
            c.close();
        }
    }

    @Override
    public boolean isCTFMd() {
        return getSetType().equals("SetOfCTF");
    }

    @Override
    public EllipseCTF getEllipseCTF(long id, int D) {
        if (!isCTFMd()) {
            return null;
        }
        try {
            EMObject emo = getEMObject(id);

            double Q0, Cs = 0, Ts = 0, kV = 0, defU = 0, defV = 0, defAngle;

            Q0 = emo.getValueDouble("_micObj._acquisition._amplitudeContrast");
            Cs = emo.getValueDouble("_micObj._acquisition._sphericalAberration");
            Ts = emo.getValueDouble("_micObj._samplingRate");
            kV = emo.getValueDouble("_micObj._acquisition._voltage");

            defU = emo.getValueDouble("_defocusU");
            defV = emo.getValueDouble("_defocusV");
            defAngle = emo.getValueDouble("_defocusAngle");
            Double downsampleFactor = emo.getValueDouble("_xmipp_CtfDownsampleFactor");
            if(downsampleFactor == null)
                downsampleFactor = 1.;
                //read params from sqlite
            
            return new EllipseCTF(id, Q0, Cs, downsampleFactor, Ts, kV, defU, defV, defAngle, D);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            throw new IllegalArgumentException(ex);
        }

    }

    @Override
    public String getCTFFile(long id) {
                String psde = getPSDEnhanced(id);
                if(psde != null)
                {
                    String ctf = psde.replace("_enhanced_psd.xmp", ".ctfparam");
                    if(new File(ctf).exists())
                        return ctf;
                }
        return null;
    }
    
    @Override
    public String getPSDFile(long id) {
        return (String) getEMObject(id).getValue("_psdFile");
    }
    
    @Override
    public String getPSDEnhanced(long id) {
        return (String) getEMObject(id).getValue("_xmipp_enhanced_psd");
    }
    
    public CTFDescription getCTFDescription(long id)
    {
        try {
            MetaData md = getEllipseCTF(id).getCTFMd();
            String file = "ctf.xmd";
            md.write(file);
            return new CTFDescription(file);
        } catch (Exception ex) {
            Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
            throw new IllegalArgumentException(ex);
        }
    }
    
    public ColumnInfo getGeoMatrixColumn(ColumnInfo ci)
    {
    	String column = "_transform._matrix";
    	if(ci.labelName.contains("."))
    		column = ci.labelName.substring(0, ci.labelName.lastIndexOf(".") + 1) + column;
        return getColumnInfo(column);
    }
    
    
    // Check if the underlying data has geometrical information
    public boolean containsGeometryInfo(String type) {
        String column = "_alignment";
        String value = null;
        if(properties == null)
        {
            if (parent == null) 
                return false; 
            else
                for(EMObject emo: parent.emobjects)
                    if(preffix.contains(emo.getId().toString()))
                        value = (String)emo.getValue(column);
        }
        else
            value = properties.get(column);

        if (value == null)
            return false;
            
        return value.equals(type);
    }
    
    
    
    public int getEnabledCount()
    {
        return enableds;
    }
    
    synchronized void loadNeighborhoodValues(int index, ColumnInfo column) throws SQLException
    {
        int neighborhood = 25;
        loadNeighborhoodValues(index, column, neighborhood);
    }
    
    synchronized void loadNeighborhoodValues(int index, ColumnInfo column, int neighborhood) throws SQLException
        {
            Connection c = null;
            PreparedStatement stmt = null;
            ResultSet rs = null;
            try {
                String name, alias;
                Object value;
//                c = getConnection(filename);
                c = getConnectionReadOnly(filename);
                String columns = column.comment;
                ColumnInfo indexci = null;
                if(column.labelName.endsWith("_filename"))
                {
                    indexci = getColumnInfo(column.labelName.replace("_filename", "_index"));
                    if(indexci != null)
                        columns += ", " + indexci.comment;
                }
                stmt = c.prepareStatement(String.format("SELECT %s FROM %sObjects where Id=?;", columns, preffix));
                EMObject emo;
                int fnIndex;
                for(int i = index; i <= index + neighborhood && i < size(); i ++)
                {
                    
                    emo = emobjects.get(i);
                    if(!emo.values.containsKey(column))
                    {
                        StopWatch.getInstance().printElapsedTime("reading " + columns + " from " + emo.getId());
                        stmt.setInt(1, emo.getId().intValue());
                        rs = stmt.executeQuery();
                        name = column.labelName;
                        alias = column.comment;
                        switch (column.type) {

                            case MetaData.LABEL_INT:
                                value = rs.getInt(alias);
                                break;
                            case MetaData.LABEL_DOUBLE:
                                value = rs.getFloat(alias);
                                break;
                            case MetaData.LABEL_STRING:
                                value = rs.getString(alias);
                                if (indexci != null) {
                                        fnIndex = rs.getInt(indexci.comment);
                                        if (fnIndex > 0) 
                                            value = fnIndex + "@" + value;
                                    
                                }
                                break;
                            default:

                                value = rs.getString(alias);
                        }
                        emo.values.put(column, value);
                    }
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
            }
            finally {
                if(rs != null)
                    rs.close();
                if (stmt != null)
                    stmt.close();
                c.close();
            }
        }
   
    public boolean isChanged()
    {
        
        for(EMObject emo: emobjects)
            if(emo.changed)
                return true;
        return false;
    }

    boolean isClassificationMd() {
    	String setType = getSetType();
    	if(setType == null)
    		return false;
        return setType.startsWith("SetOfClasses");
    }
    
    
    public class EMObject {

        int index;
        Map<ColumnInfo, Object> values;
        boolean changed;
        ScipionMetaData childmd;
        ScipionMetaData md;
        

        public EMObject(int index, ScipionMetaData md) {
            this.index = index;
            values = new HashMap<ColumnInfo, Object>();
            this.md = md;
            changed = false;
        }

        Object getValue(String column) {
            ColumnInfo ci = getColumnInfo(column);
            return getValue(ci);
        }

        Object getValue(ColumnInfo c) {
            if(c == null)
                return null;
            if(!values.containsKey(c))
                try {
                    loadNeighborhoodValues(index, c);
            } catch (SQLException ex) {
                Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
            }
            Object result = values.get(c);

            return result;
        }
        
        

        
        Double getValueDouble(String column) {
            ColumnInfo ci = getColumnInfo(column);
            return getValueDouble(ci);
        }

        Double getValueDouble(ColumnInfo c) {
            Object value = getValue(c);
            if (value != null) {
                return ((Float) value).doubleValue();//in sqlite double is float
            }
            return null;
        }

        public void setValue(ColumnInfo ci, Object value) {
            
            changed = true;
            values.put(ci, value);
            
        }
        
        public Boolean getValueBoolean(String column) {
            ColumnInfo ci = getColumnInfo(column);
            if(column == null)
                return null;
            Object value = getValue(ci);
            if(value == null)
                return null;
            return (Integer)value == 1;
        }

        public boolean isEnabled() {
            if(!values.containsKey(enabledci))
                return true;
            Integer value = (Integer) values.get(enabledci);
            return value == 1;
        }

        public void setEnabled(boolean isenabled) {
            int value = isenabled ? 1 : 0;
            if(isEnabled() && !isenabled)
                enableds --;
            else if( !isEnabled() && isenabled)
                enableds ++;
            setValue(enabledci, value);
        }

        public Long getId() {
            return ((Integer)values.get(idci)).longValue();
        }

        public String getLabel() {
            Object value = getValue(labelci);
            if (value != null) {
                return ((String) getValue(labelci));
            }
            return null;
        }

        public String getComment() {
            Object value = getValue(commentci);
            if (value != null) {
                return ((String) getValue(commentci));
            }
            return null;
        }

        public void setComment(String comment) {
            setValue(commentci, comment);
        }

        String getValueString(ColumnInfo ci) {
            Object value = getValue(ci);
            if (value != null) {
                return  value.toString();//in sqlite double is float
            }
            return "";
        }

        void setIndex(int index) {
            this.index = index;
        }

        void setLabel(String label) {
            setValue(labelci, label);
        }
        
        public EllipseCTF getEllipseCTF()
        {
            String comment = getComment();
            if(comment == null || comment.isEmpty())
                return null;
            String[] params = comment.trim().split("\\s+");
            double defU = Double.parseDouble(params[0]);
            double defV = Double.parseDouble(params[1]);
            double angle = Double.parseDouble(params[2]);
            double lowFreq = Double.parseDouble(params[3]);
            double highFreq = Double.parseDouble(params[4]);
            EllipseCTF ctf = md.getEllipseCTF(getId());
            ctf.setDefocus(defU, defV, angle);
            ctf.setFreqRange(lowFreq, highFreq);
            return ctf;
        }
    }
    
    
        

    public static Connection getConnection(String filename) throws ClassNotFoundException, SQLException
    {
        return getConnection(filename, null);
    }
    public static Connection getConnectionReadOnly(String filename) throws ClassNotFoundException, SQLException
    {

        SQLiteConfig config = new SQLiteConfig();
        config.setReadOnly(true);
        Connection c = getConnection(filename, config.toProperties());
        c.setReadOnly(true);
        return c;
    }

    private static Connection getConnection(String filename, Properties properties) throws ClassNotFoundException, SQLException {

        Class.forName("org.sqlite.JDBC");
        
        Connection c = null ;

        if (properties == null) {
            c = DriverManager.getConnection("jdbc:sqlite:" + filename);
        } else {
            c = DriverManager.getConnection("jdbc:sqlite:" + filename, properties);
        }
        return c;

    }


    void loadSelection(String selectedPath) {
        try {

            Connection c = getConnectionReadOnly(selectedPath);
            if(parent != null)
                parent.loadSelection(selectedPath, c);
            else
                loadSelection(selectedPath, c);
            c.close();
        } catch (Exception ex) {
            Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    
    void loadSelection(String selectedPath, Connection c) {
        
        try {
            
            Statement stmt = c.createStatement();
            ResultSet rs;
            EMObject emo;
            String alias;

            String query = String.format("SELECT id, enabled, label, comment FROM %sObjects;", preffix);
            rs = stmt.executeQuery(query);
            Object enabled, label, comment;
            long id;
            enableds = 0;
            boolean isctfmd = isCTFMd();
            EllipseCTF ctf;
            while (rs.next()) {
                id = rs.getInt("id");
                emo = getEMObject(id);
                alias = enabledci.comment;
                enabled = rs.getInt(alias);
                emo.setValue(enabledci, enabled);
                alias = labelci.comment;
                label = rs.getString(alias);
                emo.setValue(labelci, label);
                alias = commentci.comment;
                comment = rs.getString(alias);
                emo.setValue(commentci, comment);
                if(emo.isEnabled())
                    enableds ++;
                if(isctfmd)
                {
                    ctf = emo.getEllipseCTF();
                    if(ctf != null)
                        ctfs.put(id, ctf);
                }
                
            }
            
            rs.close();
            stmt.close();
            
            if (haschilds) {
                for (EMObject emobject : emobjects) {
                    if(emobject.childmd != null)
                        emobject.childmd.loadSelection(selectedPath, c);
                }
            }
        } catch (Exception ex) {
            Logger.getLogger(ScipionMetaData.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
   
    
    class EMObjectComparator implements Comparator<EMObject> {
        private final ColumnInfo ci;
        
        public EMObjectComparator(ColumnInfo ci)
        {
            this.ci = ci;
        }

        @Override
        public int compare(EMObject o1, EMObject o2) {
            Comparable value1 = (Comparable)o1.getValue(ci);
            Comparable value2 = (Comparable)o2.getValue(ci);
            return value1.compareTo(value2);
          }

    }
    
    public String getSetType() {
        if(properties == null)
            return self + "s";//child md from classes set
        if(properties.isEmpty())
        	return "";
        return properties.get("self");
    }
    
    public double[] getStatistics(boolean applyGeo, int label)
    {
        double[] minMax = new double [2];
        ColumnInfo ci = getColumnInfo(label);
        if(!ci.render)
            throw new IllegalArgumentException("No images to process");
        String imageFn, mddir = this.getBaseDir();
        ImagePlus imp;
        minMax[0] = Double.MAX_VALUE;
        minMax[1] = -Double.MAX_VALUE;
        double aux;
        for (EMObject emo: emobjects)
        {
            imageFn = emo.getValueString(ci);
            if(imageFn != null)
            {
                imageFn = Filename.findImagePath(imageFn, mddir, true);
                if (imageFn != null && Filename.exists(imageFn))
                {
                    imp = new ImagePlusFromFile(imageFn).getImagePlus();
                    aux = imp.getProcessor().getMin();
                    if(aux < minMax[0])
                        minMax[0] = aux;
                    aux = imp.getProcessor().getMax();
                    if(aux > minMax[1])
                        minMax[1] = aux;
                }
            }
        }
        
        return minMax;
    }
    
    public boolean containsMicrographsInfo() 
    {
    	return self.equals("Micrograph") || self.equals("Movie");
    }

	public String getSelf()
	{
		// TODO Auto-generated method stub
		return self;
	}
	
	public static boolean isScipionMetaData(String filename)
	{
		if (!filename.endsWith(".sqlite") && !filename.endsWith(".db"))
				return false; 
		//Now it will check if contains the required tables to be read as scipion metadata, otherwise will be read as xmipp metadata
		Connection c;
		Statement stmt;
		ResultSet rs;
		boolean result = true;
		try
		{
//			c = getConnection(filename);
            c = getConnectionReadOnly(filename);
			stmt = c.createStatement();
			String query = "SELECT name FROM sqlite_master WHERE type='table' AND name='Classes';";
            rs = stmt.executeQuery(query);
    		boolean exists = rs.next();
    		if (!exists)
    			result = false;
    		else
    		{
	    		query = "SELECT name FROM sqlite_master WHERE type='table' AND name='Objects';";
	            rs = stmt.executeQuery(query);
	    		exists = rs.next();
	    		if (!exists)
	    			result = false;
    		}
            rs.close();
            stmt.close();
            c.close();
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
		
        return result;
	}
   
}
