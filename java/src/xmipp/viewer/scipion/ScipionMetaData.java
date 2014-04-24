/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.sql.*;

/**
 *
 * @author airen
 */
public class ScipionMetaData {
    private String dbfile;
    
    public ScipionMetaData(String dbfile)
    {
        this.dbfile = dbfile;
        Connection c = null;
        Statement stmt = null;
        try {
            Class.forName("org.sqlite.JDBC");
            c = DriverManager.getConnection("jdbc:sqlite:" + dbfile);
            stmt = c.createStatement();
            ResultSet rs = stmt.executeQuery( "SELECT * FROM CLASSES;" );
            String label, column;
            while ( rs.next() ) {
               label = rs.getString("label_property");
               column = rs.getString("column_name");
               System.out.println( "label = " + label );
               System.out.println( "column = " + column );
               System.out.println();
            }
            rs.close();
            stmt.close();
            c.close();
        } 
        catch ( Exception e ) {
            System.err.println( e.getClass().getName() + ": " + e.getMessage() );
            
        }
        System.out.println("Opened database successfully");
    }
    
    public static void main( String args[] )
    {
        String dbfile = "/home/airen/pyworkflow-code/data/tests/xmipp_tutorial/gold/micrographs_gold.sqlite";
        ScipionMetaData md = new ScipionMetaData(dbfile);
    }
    
}
