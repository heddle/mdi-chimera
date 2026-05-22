package edu.cnu.mdi.chimera.util;

/**
 * Implements a rotation via a 3x3 matrix and its inverse.
 * <p>
 * The rotation can be built from up to three single-axis rotations.
 * The axes are specified using the constants CR_X_AXIS, CR_Y_AXIS, and CR_Z_AXIS.
 * </p>
 */
public class ChimeraRotation {

    // Constants to specify the axes.
    public static final int CR_X_AXIS = 0;
    public static final int CR_Y_AXIS = 1;
    public static final int CR_Z_AXIS = 2;

    // A tiny angle for numerical comparisons.
    public static final double TINY_ANG = 1.0e-8;

    // The rotation matrix (3x3)
    public double[][] matrix = new double[3][3];
    // The inverse rotation matrix (3x3)
    public double[][] invMatrix = new double[3][3];
    // A flag indicating if this rotation is the identity rotation.
    public boolean identity = false;

    /**
     * Default constructor. Creates an identity rotation.
     */
    public ChimeraRotation() {
        setIdentity();
    }

    // --------------------------
    // Rotation creation methods that return new objects
    // --------------------------

    /**
     * Creates a rotation from a single axis rotation.
     *
     * @param theta1 The rotation angle (radians) about the given axis.
     * @param axis1  The axis of rotation (use CR_X_AXIS, CR_Y_AXIS, or CR_Z_AXIS).
     * @return A new ChimeraRotation corresponding to the given rotation.
     */
    public static ChimeraRotation rotationFromAngles(double theta1, int axis1) {
        return rotationFromAngles(theta1, axis1, 0.0, CR_X_AXIS);
    }

    /**
     * Creates a rotation from two sequential rotations.
     *
     * @param theta1 The first rotation angle (radians) about axis1.
     * @param axis1  The first rotation axis.
     * @param theta2 The second rotation angle (radians) about axis2.
     * @param axis2  The second rotation axis.
     * @return A new ChimeraRotation corresponding to the two rotations.
     */
    public static ChimeraRotation rotationFromAngles(double theta1, int axis1,
                                                       double theta2, int axis2) {
        return rotationFromAngles(theta1, axis1, theta2, axis2, 0.0, CR_X_AXIS);
    }

    /**
     * Creates a rotation from up to three sequential rotations.
     * <p>
     * If theta2 (or theta3) is effectively zero (within TINY_ANG), that rotation is skipped.
     * </p>
     *
     * @param theta1 The first rotation angle (radians) about axis1.
     * @param axis1  The first rotation axis.
     * @param theta2 The second rotation angle (radians) about axis2.
     * @param axis2  The second rotation axis.
     * @param theta3 The third rotation angle (radians) about axis3.
     * @param axis3  The third rotation axis.
     * @return A new ChimeraRotation representing the combined rotation.
     */
    public static ChimeraRotation rotationFromAngles(double theta1, int axis1,
                                                       double theta2, int axis2,
                                                       double theta3, int axis3) {
        ChimeraRotation rotation = new ChimeraRotation();

        // Check for near-zero rotations.
        boolean t1Zero = Math.abs(theta1) < TINY_ANG;
        boolean t2Zero = Math.abs(theta2) < TINY_ANG;
        boolean t3Zero = Math.abs(theta3) < TINY_ANG;

        if (t1Zero && t2Zero && t3Zero) {
            rotation.setIdentity();
            return rotation;
        }
        else {
            rotation.identity = false;
            // First, create the single-axis rotation for theta1.
            rotation.matrix = singleAxisRotation(theta1, axis1);

            // If theta2 is nonzero, combine that rotation.
            if (!t2Zero) {
                double[][] m2 = singleAxisRotation(theta2, axis2);
                // If theta3 is provided and nonzero, combine that as well.
                if (!t3Zero) {
                    double[][] m3 = singleAxisRotation(theta3, axis3);
                    m2 = matrixMultiply(m3, m2);
                }
                // Combine m2 with the first rotation: newMatrix = m2 * (theta1 rotation).
                rotation.matrix = matrixMultiply(m2, rotation.matrix);
            }
            // Compute the inverse matrix as the transpose of the rotation matrix.
            rotation.invMatrix = transpose(rotation.matrix);
        }
        return rotation;
    }

    // --------------------------
    // Overloaded method to set an existing rotation
    // --------------------------

    /**
     * Sets the provided ChimeraRotation instance according to two sequential rotations.
     * <p>
     * This method is provided to match the call signature from Plane.
     * </p>
     *
     * @param r       The ChimeraRotation instance to set.
     * @param theta1  The first rotation angle (radians).
     * @param axis1   The first rotation axis.
     * @param theta2  The second rotation angle (radians).
     * @param axis2   The second rotation axis.
     */
    public static void rotationFromAngles(ChimeraRotation r, double theta1, int axis1,
                                          double theta2, int axis2) {
        // Use the two-rotation version that returns a new instance.
        ChimeraRotation temp = rotationFromAngles(theta1, axis1, theta2, axis2);
        r.identity = temp.identity;
        r.matrix = temp.matrix;
        r.invMatrix = temp.invMatrix;
    }

    // --------------------------
    // Utility methods
    // --------------------------

    /**
     * Sets this rotation to the identity rotation.
     */
    public void setIdentity() {
        identity = true;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                matrix[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        // The inverse of the identity is itself.
        invMatrix = copyMatrix(matrix);
    }

    /**
     * Creates a single-axis rotation matrix.
     *
     * @param theta The rotation angle (radians).
     * @param axis  The axis of rotation (CR_X_AXIS, CR_Y_AXIS, or CR_Z_AXIS).
     * @return A 3x3 rotation matrix.
     */
    public static double[][] singleAxisRotation(double theta, int axis) {
        double sin_t = Math.sin(theta);
        double cos_t = Math.cos(theta);
        double[][] m = new double[3][3];

        switch (axis) {
            case CR_X_AXIS:
                m[0][0] = 1.0;   m[0][1] = 0.0;    m[0][2] = 0.0;
                m[1][0] = 0.0;   m[1][1] = cos_t;  m[1][2] = sin_t;
                m[2][0] = 0.0;   m[2][1] = -sin_t; m[2][2] = cos_t;
                break;
            case CR_Y_AXIS:
                m[0][0] = cos_t;  m[0][1] = 0.0; m[0][2] = -sin_t;
                m[1][0] = 0.0;    m[1][1] = 1.0; m[1][2] = 0.0;
                m[2][0] = sin_t;  m[2][1] = 0.0; m[2][2] = cos_t;
                break;
            case CR_Z_AXIS:
                m[0][0] = cos_t;  m[0][1] = sin_t; m[0][2] = 0.0;
                m[1][0] = -sin_t; m[1][1] = cos_t; m[1][2] = 0.0;
                m[2][0] = 0.0;    m[2][1] = 0.0;   m[2][2] = 1.0;
                break;
            default:
                throw new IllegalArgumentException("[ChimeraRotation] Illegal axis constant for rotation: " + axis);
        }
        return m;
    }

    /**
     * Multiplies two 3x3 matrices: result = A * B.
     *
     * @param A The first matrix.
     * @param B The second matrix.
     * @return The product matrix.
     */
    public static double[][] matrixMultiply(double[][] A, double[][] B) {
        double[][] result = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result[i][j] = 0.0;
                for (int k = 0; k < 3; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }

    /**
     * Returns the transpose of a 3x3 matrix.
     *
     * @param A The matrix to transpose.
     * @return The transposed matrix.
     */
    public static double[][] transpose(double[][] A) {
        double[][] T = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                T[i][j] = A[j][i];
            }
        }
        return T;
    }

    /**
     * Returns a copy of a 3x3 matrix.
     *
     * @param A The matrix to copy.
     * @return A new copy of the matrix.
     */
    public static double[][] copyMatrix(double[][] A) {
        double[][] copy = new double[3][3];
        for (int i = 0; i < 3; i++) {
            System.arraycopy(A[i], 0, copy[i], 0, 3);
        }
        return copy;
    }

    /**
     * Multiplies a 3x3 matrix by a 3-element vector.
     *
     * @param A      The 3x3 matrix.
     * @param vector A 3-element vector.
     * @return The resulting 3-element vector.
     */
    public static double[] multiplyMatrixVector(double[][] A, double[] vector) {
        double[] result = new double[3];
        for (int i = 0; i < 3; i++) {
            result[i] = A[i][0] * vector[0] + A[i][1] * vector[1] + A[i][2] * vector[2];
        }
        return result;
    }

    // --------------------------
    // Debug and rotation operations
    // --------------------------

    /**
     * Dumps the rotation matrix, its inverse, and their product to standard output.
     */
    public void rotationDump() {
        System.out.println("Rotation Matrix: ");
        printMatrix(matrix, 14);
        System.out.println("Inverse Rotation Matrix: ");
        printMatrix(invMatrix, 14);
        System.out.println("Product Matrix: ");
        printMatrix(matrixMultiply(matrix, invMatrix), 14);
    }

    /**
     * Prints a 3x3 matrix with the specified field width.
     *
     * @param A          The matrix to print.
     * @param fieldWidth The field width for each element.
     */
    public static void printMatrix(double[][] A, int fieldWidth) {
        for (int i = 0; i < 3; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < 3; j++) {
                sb.append(String.format("%" + fieldWidth + ".6f ", A[i][j]));
            }
            System.out.println(sb.toString());
        }
    }

    /**
     * Applies this rotation to the given point.
     * <p>
     * If this rotation is the identity, the point is unchanged.
     * </p>
     *
     * @param v The point to rotate (modified in place).
     */
    public void cv_rotate(Point3D.Double v) {
        if (identity) {
            return;
        }
        double[] vec = new double[]{v.x, v.y, v.z};
        double[] res = multiplyMatrixVector(matrix, vec);
        v.setLocation(res[0], res[1], res[2]);
    }

    /**
     * Applies the inverse of this rotation to the given point.
     * <p>
     * If this rotation is the identity, the point is unchanged.
     * </p>
     *
     * @param v The point to rotate (modified in place).
     */
    public void inverseRotate(Point3D.Double v) {
        if (identity) {
            return;
        }
        double[] vec = new double[]{v.x, v.y, v.z};
        double[] res = multiplyMatrixVector(invMatrix, vec);
        v.setLocation(res[0], res[1], res[2]);
    }
}
