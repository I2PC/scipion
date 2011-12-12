package sphere;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Geometry {

    public final static int MAX_ROT = 360;
    public final static int MAX_TILT = 180;
    private final static double FLT_EPSILON = 1.19209e-07;

    /**
     * Returns the 3D point associated to rot and tilt.
     * @param rot
     * @param tilt
     * @return
     */
    public static Point3d getSphereCoordinates(double rot, double tilt) {
        double rotAsRads = Math.toRadians(rot);
        double tiltAsRads = Math.toRadians(tilt);

        return new Point3d(
                Math.sin(tiltAsRads) * Math.cos(rotAsRads),
                Math.sin(tiltAsRads) * Math.sin(rotAsRads),
                Math.cos(tiltAsRads));
    }
    
    public static double[] getAngles(double matrix[]) {
        // Gets alpha, beta...
        Point3d point = new Point3d(matrix[2], -matrix[6], -matrix[10]);
        double angles[] = getAngles(point);

        // ...and psi.
        
        
        return new double[]{angles[0], angles[1], 0};
    }

    public static double[] getAngles(Point3d point) {
        double alpha, beta;
        double abs_ca, abs_sa, sb, cb;
        double aux_alpha, aux_beta;
        double error, newerror;
        Vector3d v_aux = new Vector3d();
        Vector3d v = new Vector3d(point);

        //if not normalized do it so
        v.normalize();

        cb = v.z;

        /*        if (Math.abs((cb)) > 0.999847695) {//one degree
        System.out.println(
        "\tWARNING: Routine Euler_direction2angles is not reliable"
        + "for small tilt angles. Up to 0.001 deg it should be OK\n"
        + "for most applications but you never know");
        }*/

        if (Math.abs((cb - 1.)) < FLT_EPSILON) {
            alpha = 0.;
            beta = 0.;
        } else {//1

            aux_beta = Math.acos(cb); // beta between 0 and PI

            sb = Math.sin(aux_beta);

            abs_ca = Math.abs(v.x) / sb;
            if (Math.abs((abs_ca - 1.)) < FLT_EPSILON) {
                aux_alpha = 0.;
            } else {
                aux_alpha = Math.acos(abs_ca);
            }

            v_aux.x = Math.sin(aux_beta) * Math.cos(aux_alpha);
            v_aux.y = Math.sin(aux_beta) * Math.sin(aux_alpha);
            v_aux.z = Math.cos(aux_beta);

            error = Math.abs(v.dot(v_aux) - 1.);
            alpha = aux_alpha;
            beta = aux_beta;

            v_aux.x = Math.sin(aux_beta) * Math.cos(-1. * aux_alpha);
            v_aux.y = Math.sin(aux_beta) * Math.sin(-1. * aux_alpha);
            v_aux.z = Math.cos(aux_beta);
            newerror = Math.abs(v.dot(v_aux) - 1.);
            if (error > newerror) {
                alpha = -1. * aux_alpha;
                beta = aux_beta;
                error = newerror;
            }

            v_aux.x = Math.sin(-aux_beta) * Math.cos(-1. * aux_alpha);
            v_aux.y = Math.sin(-aux_beta) * Math.sin(-1. * aux_alpha);
            v_aux.z = Math.cos(-aux_beta);
            newerror = Math.abs(v.dot(v_aux) - 1.);
            if (error > newerror) {
                alpha = -1. * aux_alpha;
                beta = -1. * aux_beta;
                error = newerror;
            }

            v_aux.x = Math.sin(-aux_beta) * Math.cos(aux_alpha);
            v_aux.y = Math.sin(-aux_beta) * Math.sin(aux_alpha);
            v_aux.z = Math.cos(-aux_beta);
            newerror = Math.abs(v.dot(v_aux) - 1.);

            if (error > newerror) {
                alpha = aux_alpha;
                beta = -1. * aux_beta;
                error = newerror;
            }
        }//else 1 end

        beta = Math.toDegrees(beta);
        alpha = Math.toDegrees(alpha);

        return new double[]{alpha, beta};
    }

    /**
     * rot, tilt -> 3D point -> abs -> rot', tilt'
     * @param rot
     * @param tilt
     * @return
     */
    public static double[] normalizeAngles(double rot, double tilt) throws Exception {
        //double angles[] = absAngles(rot, tilt); // Avoids negative values
        Point3d point = getSphereCoordinates(rot, tilt); // Associated point in sphere 3D space.
        //double direction[] = new double[]{point.x, point.y, point.z};

        //double angles[] = Projection.eulerDirection2Angles(direction);//getAngles(point);  // Previous operation inversed: from 3D point back to angles.
        double angles[] = getAngles(point);  // Previous operation inversed: from 3D point back to angles.
        
        angles = absAngles(angles[0], angles[1]);   // Avoids negative values again.

        return new double[]{angles[0], angles[1]};
    }

    /**
     * Computes angles absolute values.
     * @param point
     */
    public static double[] absAngles(double rot, double tilt) {
        double abs_tilt = tilt;
        double abs_rot = rot;

        if (abs_tilt < 0) { // tilt <0
            abs_tilt *= -1; // -tilt
            abs_rot += 180; // 180 +- rot
        } else {    // tilt > 0
            if (abs_rot < 0) {  // rot < 0
                abs_rot += 360; // 360 +- rot
            }
        }

//        System.out.println(rot + ", " + tilt + " -> " + (abs_rot % 360) + ", " + abs_tilt);

        return new double[]{abs_rot % 360, abs_tilt};
    }

    /**
     * Fixes point so rot (x) is circular [0-360) ant tilt (y) is [0-180]
     * @param point
     */
    public static void fixPoint(Point2d point) {
        if (point.x < 0) {  // Circular: <- Rot ->
            point.x += MAX_ROT;
        } else if (point.x >= MAX_ROT) {
            point.x -= MAX_ROT;
        }

        if (point.y < 0) {  // Not circular: |<- Tilt ->|
            point.y = 0;
        } else if (point.y > MAX_TILT) {
            point.y = MAX_TILT;
        }

        point.x %= MAX_ROT;
    }

    /**
     * Returns length of point's associated vector
     * @param point
     * @return
     */
    public static double length(Point3d point) {
        Vector3d v = new Vector3d(point);

        return v.length();
    }
}
