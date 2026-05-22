package edu.cnu.mdi.chimera.util;


/**
 * A class representing a plane defined by the equation
 * Ax + By + Cz + D = 0.
 * <p>
 * The plane is defined by three points (p0, p1, p2) and maintains
 * a rotation matrix (rmat) used to align the plane's normal with the z-axis.
 * </p>
 */
public class ChimeraPlane {

    // Constants used for tolerance and limiting values.
    private static final double TOL = 1e-8;

    // Plane parameters for the equation Ax + By + Cz + D = 0.
    public double A = 0.0;
    public double B = 0.0;
    public double C = 0.0;
    public double D = 0.0;

    // Points used to create the plane.
    public Point3D.Double p0 = new Point3D.Double();
    public Point3D.Double p1 = new Point3D.Double();
    public Point3D.Double p2 = new Point3D.Double();

    // The rotation matrix (or rotation representation) used to align the normal with the z-axis.
    public ChimeraRotation rmat = new ChimeraRotation();

    /**
     * Constructs a plane using three points.
     *
     * @param p0 One point of the plane.
     * @param p1 A second point of the plane.
     * @param p2 A third point of the plane.
     */
    public ChimeraPlane(Point3D.Double p0, Point3D.Double p1, Point3D.Double p2) {
        setPlane(p0, p1, p2);
    }

    /**
     * Sets the plane using three points. This method computes the plane's
     * parameters, stores the defining points, and computes the rotation matrix.
     *
     * @param p0 One point of the plane.
     * @param p1 A second point of the plane.
     * @param p2 A third point of the plane.
     */
    public void setPlane(Point3D.Double p0, Point3D.Double p1, Point3D.Double p2) {
        getPlaneParameters(p0, p1, p2);
        // Make copies of the points.
        this.p0 = new Point3D.Double(p0.x, p0.y, p0.z);
        this.p1 = new Point3D.Double(p1.x, p1.y, p1.z);
        this.p2 = new Point3D.Double(p2.x, p2.y, p2.z);
        getRotationMatrix();
    }

    /**
     * Computes the plane parameters (A, B, C, D) from three points.
     *
     * @param p0 One point of the plane.
     * @param p1 A second point of the plane.
     * @param p2 A third point of the plane.
     */
    private void getPlaneParameters(Point3D.Double p0, Point3D.Double p1, Point3D.Double p2) {
        double dx1 = p1.x - p0.x;
        double dy1 = p1.y - p0.y;
        double dz1 = p1.z - p0.z;

        double dx2 = p2.x - p0.x;
        double dy2 = p2.y - p0.y;
        double dz2 = p2.z - p0.z;

        // Compute cross-product components to get the normal.
        double dyz = dy1 * dz2 - dy2 * dz1;
        double dzx = dz1 * dx2 - dz2 * dx1;
        double dxy = dx1 * dy2 - dx2 * dy1;

        A = dyz;
        B = dzx;
        C = dxy;

        D = -(A * p0.x + B * p0.y + C * p0.z);
    }

    /**
     * Computes the rotation matrix that aligns the plane's normal with the z-axis.
     */
    public void getRotationMatrix() {
        double len = Math.sqrt(A * A + B * B + C * C);
        // Clamp the ratio to the interval [-1, 1] to avoid numerical issues.
        double ratio = Math.max(-1, Math.min(1, C / len));
        double theta = Math.acos(ratio);

        double phi;
        if ((Math.abs(B) < TOL) && (Math.abs(A) < TOL)) {
            phi = 0;
        } else {
            phi = Math.atan2(B, A);
        }

        // Assumes ChimeraRotation.rotationFromAngles sets up the rotation matrix appropriately.
        ChimeraRotation.rotationFromAngles(rmat, phi, ChimeraRotation.CR_Z_AXIS, theta, ChimeraRotation.CR_Y_AXIS);
    }

    /**
     * Returns a pair of matrices used for rotating points so that the new z-axis
     * aligns with the plane's normal.
     * <p>
     * The returned 3D array contains:
     * <ul>
     *   <li>index 0: the rotation matrix (double[3][3])</li>
     *   <li>index 1: the inverse rotation matrix (double[3][3])</li>
     * </ul>
     * </p>
     *
     * @return a 3D array containing the rotation matrix and its inverse.
     */
    public double[][][] getRotationMatrices() {
        // Copy the matrices so that the internal state is not exposed.
        double[][] rotMat = ChimeraRotation.copyMatrix(rmat.matrix);
        double[][] rotMatInv = ChimeraRotation.copyMatrix(rmat.invMatrix);
        return new double[][][] { rotMat, rotMatInv };
    }

    /**
     * Returns a unit normal to this plane.
     *
     * @return a Point3D.Double representing the unit normal vector.
     * @throws IllegalStateException if the normal vector has zero length.
     */
    public Point3D.Double getNormal() {
        double len = Math.sqrt(A * A + B * B + C * C);
        if (len == 0) {
            throw new IllegalStateException("Cannot compute unit normal for a zero-length normal vector.");
        }
        return new Point3D.Double(A / len, B / len, C / len);
    }
    /**
     * Checks whether a given point lies in the plane (within a small tolerance).
     *
     * @param p The point to test.
     * @return True if the point lies in the plane; otherwise, false.
     */
    public boolean isIn(Point3D.Double p) {
        return Math.abs(deviation(p)) < TOL;
    }

    /**
     * Computes the deviation of a point from the plane.
     *
     * @param p The point to test.
     * @return The deviation (should be near 0 if the point is in the plane).
     */
    public double deviation(Point3D.Double p) {
        return (A * p.x + B * p.y + C * p.z + D);
    }


    /**
     * Computes the "length" of the plane, defined as sqrt(A^2 + B^2 + C^2).
     *
     * @return The computed length.
     */
    public double planeLength() {
        return Math.sqrt(A * A + B * B + C * C);
    }

    /**
     * Computes the value of z' (the constant z value in the rotated frame).
     *
     * @return The computed zprime value.
     */
    public double zprime() {
        return -D / planeLength();
    }



}

