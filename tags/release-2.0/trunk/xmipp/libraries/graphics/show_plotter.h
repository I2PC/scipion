/***************************************************************************
 *
 * Author:     Javier Rodriguez Fernandez (javrodri@gmail.com)
 *             Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Universidad San Pablo CEU (Monteprincipe, Madrid)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef SHOWPLOTTER_H
#define SHOWPLOTTER_H

#include <vector>

#include <qevent.h>
#include <qlabel.h>
#include <qmainwindow.h>
#include <qpixmap.h>
#include <qtoolbutton.h>
#include <qwidget.h>

#include <data/matrix2d.h>

/**@defgroup ShowPlotter Show Plotter
   @ingroup GraphicsLibrary */
//@{
/** Class to store the Settings of the plotter.
    It has been implemented to let the plotter
    to have different settings, used for the zoom.
*/
class PlotSettings
{
public:
    /// Minimum value showed in X axis
    double minX;
    /// Maximum value showed in X axis
    double maxX;
    /** Number of Xticks.
        By convention numXTicks are off by one,
        if numXTicks is 5, plotter will draw 6 ticks marks */
    int numXTicks;

    /// Minimum value showed in Y axis
    double minY;
    /// Maximum value showed in Y axis
    double maxY;
    /** Number of Yticks.
        By convention numXTicks are off by one,
        if numXTicks is 5, plotter will draw 6 ticks marks */
    int numYTicks;

    /// X position of the zoom
    double posX;
    /// Y position of the zoom
    double posY;
public:
    /** Empty Constructor */
    PlotSettings();
public:
    // Internal functions
    // Function to scroll settings, when we move through the plotter
    void scroll(int dx, int dy);

    // Function to set nice values at the axis
    void adjust();

    // X range
    double spanX() const
    {
        return maxX - minX;
    }

    // Y range
    double spanY() const
    {
        return maxY - minY;
    }

    // Procedure to take nice values for the different Ticks of the axis and
    // for the min and max values of the axis.
    void adjustAxis(double &min, double &max);
};

/** Class to represent pairs of values.
    It may receive a Matrix2D object or two matrices 1D;
    It has some utilities like ZoomIn & Out, move through it, save plotter as
    .png image.

    An example of use is provided by:
    @code
         #include <XmippGraphics/showPlotter.hh>
         #include <qapplication.h>

         int main(int argc, char *argv[]) {
            try {
               // Data to plot
               Matrix1D<double> y;
               y.initLinear(1,100);

               // Setup the interface
               QApplication app(argc,argv);
               Plotter plotter;
               plotter.setCurveData(0,y);
               plotter.setCurveData(1,y/2);
               app.setMainWidget(&plotter);
               return app.exec();
            } catch (Xmipp_error XE) {cout << XE;}
            return 0;
         }
    @endcode
*/
class Plotter : public QMainWindow
{
    Q_OBJECT
public:
    /** Empty constructor */
    Plotter(QWidget *parent = 0, const char *name = 0);

    /** Destructor */
    ~Plotter();

    /** Assign a new curve to the plotter.
        It needs an ID for the curve and a Matrix1D with the curve.
        it will be assigned to CurveMap and shown in the plotter.*/
    void setCurveData(int id, const Matrix1D<double> &Y);

    /** Assign a new curve to the plotter. */
    void setCurveData(int id, const Matrix1D<double> &X,
                      const Matrix1D<double> &Y);

    /** Assign a new curve to the plotter.
        The first column is assumed to be the X axis and the second column
        the Y axis.
    */
    void setCurveData(int id, const Matrix2D<double> &data);

    /** Delete a curve from our plotter. */
    void deleteCurve(int id);

    /** Refresh pixmap */
    void refreshPixmap();

    /** Set active state of the given curve */
    void setActiveState(int id, bool new_state)
    {
        curveActive[id] = new_state;
    }

    /** Set plotsettings */
    void setPlotSettings(const PlotSettings &new_settings);

    /** Copy from another Plotter */
    void copy(const Plotter &plotter);
signals:
    void resizeDone();

public slots:
    void zoomIn();
    void zoomOut();
    void saveToFile();

public:
    // Zoom related
    std::vector<PlotSettings> zoomStack;
    int curZoom;

    // Events
    void paintEvent(QPaintEvent *event);
    void resizeEvent(QResizeEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);
    void mouseMoveStatus(QMouseEvent *event);

    // Status bar related
    void updateCellIndicators(double x, double y);
    void updateModLabel(const QString & message);
    QLabel *locationLabelX;
    QLabel *locationLabelY;
    QLabel *formulaLabel;
    QLabel *modLabel;

    // Buttons and related actions
    QToolButton *zoomInButton;
    QToolButton *zoomOutButton;
    QToolButton *saveButton;

    // Drawing
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void drawGrid(QPainter *painter);
    void updateRubberBandRegion();
    void refreshCurves();
    void drawCurves(QPainter *painter);
    void createStatusBar();
    QPixmap pixmap;
    static const int Margin = 50;
    bool rubberBandIsShown;
    QRect rubberBandRect;

    // Data to plot
    std::map<int, Matrix2D<double> > curveMap;
    std::map<int, bool > curveActive;
};
//@}
#endif // SHOWPLOTTER_H
