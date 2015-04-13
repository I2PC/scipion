/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package m;

import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 *
 * @author roberto
 */
@WebServlet(name = "ReportInstallation", urlPatterns = {"/ReportInstallation"})
public class ReportInstallation extends HttpServlet {

    private static final boolean glassfishprod = true; //ifdef equivalent for java

    /**
     * Processes requests for both HTTP
     * <code>GET</code> and
     * <code>POST</code> methods.
     *
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    protected void processRequest(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        response.setContentType("text/html;charset=UTF-8");
        PrintWriter out = response.getWriter();
        try {
            String remoteAddr;
            if (glassfishprod) {
                //remoteAddr is added by apache (on i2pc) before redirecting
                remoteAddr = request.getParameter("remoteAddr");
            } else {
                //if you want to check this in your machine set debug to false. Note: glassdishdev does not add remoteAddr
                remoteAddr = request.getRemoteAddr();
            }
            
            String osystem = request.getParameter("osystem");
            Connection c = null;
            
            // doNotStore allows you to check that the server is up without filling the database
            if (osystem.equals("doNotStore")) {
                out.print("doNotStore");
            } else {
                try {
                    //load sqlite driver
                    Class.forName("org.sqlite.JDBC");
                    Connection connection;
                    if (glassfishprod) {
                        connection = DriverManager.getConnection("jdbc:sqlite:/home/databases/xmipp/xmippInstall.sqlite");
                    } else {
                        connection = DriverManager.getConnection("jdbc:sqlite:/tmp/xmippInstall.sqlite");
                    }
                    // prepare statement so SQL injection is harder
                    String insertStatement = "INSERT INTO xmipp_installations (operating_system, remoteAddr) "
                            + "VALUES (?,?);";
                    PreparedStatement prepStmt = connection.prepareStatement(insertStatement);
                    
                    // no buffer overflows (sanity check)
                    osystem = osystem.length() > 199 ? osystem.substring(0, 198) : osystem;
                    remoteAddr = remoteAddr.length() > 16 ? remoteAddr.substring(0, 15) : remoteAddr;
                    
                    //parameter are now sanitized, set them
                    prepStmt.setString(1, osystem);
                    prepStmt.setString(2, remoteAddr);
                    prepStmt.executeUpdate();

                    prepStmt.close();
                    connection.close();
                } catch (Exception e) {
                    System.err.println(e.getMessage());
                }
            }
        } finally {
            out.close();
        }
    }

    // <editor-fold defaultstate="collapsed" desc="HttpServlet methods. Click on the + sign on the left to edit the code.">
    /**
     * Handles the HTTP
     * <code>GET</code> method.
     *
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        processRequest(request, response);
    }

    /**
     * Handles the HTTP
     * <code>POST</code> method.
     *
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        processRequest(request, response);
    }

    /**
     * Returns a short description of the servlet.
     *
     * @return a String containing servlet description
     */
    @Override
    public String getServletInfo() {
        return "Short description";
    }// </editor-fold>
}
