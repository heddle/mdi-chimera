package edu.cnu.mdi.chimera.grid;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Spherical surface grid aligned with GSM coordinates.
 * <p>
 * The grid is defined on the surface of a sphere centered at Earth center.
 * The polar angle {@code theta} is measured from +GSM Z, and the azimuthal
 * angle {@code phi} is measured in the GSM X-Y plane from +GSM X toward
 * +GSM Y.
 * </p>
 */
public final class SphericalGrid {

    /** Tolerance used for angular endpoint validation. */
    private static final double ANGLE_TOL = 1.0e-6;

    /** Theta grid vertices, spanning 0 to pi. */
    private final Grid1D thetaGrid;

    /** Phi grid vertices, spanning -pi to pi. */
    private final Grid1D phiGrid;

    /** Sphere radius in the owning length unit, normally Earth radii. */
    private final double radius;

    /**
     * Creates a spherical grid.
     *
     * @param thetaArray theta grid vertices in radians, spanning 0 to pi
     * @param phiArray phi grid vertices in radians, spanning -pi to pi
     * @param radius sphere radius
     */
    public SphericalGrid(double[] thetaArray, double[] phiArray, double radius) {
        if (radius <= 0.0 || !Double.isFinite(radius)) {
            throw new IllegalArgumentException("Spherical grid radius must be positive and finite.");
        }

        double[] thetaCopy = Arrays.copyOf(thetaArray, thetaArray.length);
        double[] phiCopy = Arrays.copyOf(phiArray, phiArray.length);

        Arrays.sort(thetaCopy);
        Arrays.sort(phiCopy);

        if (Math.abs(thetaCopy[0]) > ANGLE_TOL) {
            throw new IllegalArgumentException("Theta grid must start at 0.");
        }
        if (Math.abs(thetaCopy[thetaCopy.length - 1] - Math.PI) > ANGLE_TOL) {
            throw new IllegalArgumentException("Theta grid must end at pi.");
        }
        if (Math.abs(phiCopy[0] + Math.PI) > ANGLE_TOL) {
            throw new IllegalArgumentException("Phi grid must start at -pi.");
        }
        if (Math.abs(phiCopy[phiCopy.length - 1] - Math.PI) > ANGLE_TOL) {
            throw new IllegalArgumentException("Phi grid must end at pi.");
        }

        this.thetaGrid = new Grid1D(thetaCopy);
        this.phiGrid = new Grid1D(phiCopy);
        this.radius = radius;
    }

    /**
     * Copy constructor.
     *
     * @param other the grid to copy
     */
    public SphericalGrid(SphericalGrid other) {
        this(requireNonNull(other).thetaGrid.getPoints(), other.phiGrid.getPoints(), other.radius);
    }

    /**
     * Gets the theta grid.
     *
     * @return a defensive copy of the theta grid
     */
    public Grid1D getThetaGrid() {
        return new Grid1D(thetaGrid);
    }

    /**
     * Gets the phi grid.
     *
     * @return a defensive copy of the phi grid
     */
    public Grid1D getPhiGrid() {
        return new Grid1D(phiGrid);
    }

    /**
     * Gets the sphere radius.
     *
     * @return the sphere radius
     */
    public double getRadius() {
        return radius;
    }

    /**
     * Gets the number of theta cells.
     *
     * @return the number of theta cells
     */
    public int getNumThetaCells() {
        return thetaGrid.numCells();
    }

    /**
     * Gets the number of phi cells.
     *
     * @return the number of phi cells
     */
    public int getNumPhiCells() {
        return phiGrid.numCells();
    }

    /**
     * Locates angular coordinates in the spherical grid.
     *
     * @param theta polar angle in radians
     * @param phi azimuthal angle in radians
     * @param indices output array of length at least two
     * @return the supplied indices array, filled with theta and phi cell indices
     */
    public int[] getIndices(double theta, double phi, int[] indices) {
        if (indices == null || indices.length < 2) {
            throw new IllegalArgumentException("indices must have length at least 2.");
        }

        indices[0] = thetaGrid.locateInterval(theta);
        indices[1] = phiGrid.locateInterval(normalizePhi(phi));
        return indices;
    }

    /**
     * Converts Cartesian GSM coordinates to spherical coordinates.
     *
     * @param x GSM x
     * @param y GSM y
     * @param z GSM z
     * @return a two-element array containing theta and phi in radians
     */
    public static double[] cartesianToThetaPhi(double x, double y, double z) {
        double r = Math.sqrt(x * x + y * y + z * z);
        if (r == 0.0) {
            throw new IllegalArgumentException("Cannot compute theta/phi for the origin.");
        }

        double theta = Math.acos(z / r);
        double phi = normalizePhi(Math.atan2(y, x));
        return new double[] { theta, phi };
    }

    /**
     * Converts spherical coordinates on this grid's radius to Cartesian GSM.
     *
     * @param theta polar angle
     * @param phi azimuthal angle
     * @return a three-element array containing x, y, z
     */
    public double[] thetaPhiToCartesian(double theta, double phi) {
        double sinTheta = Math.sin(theta);
        double x = radius * sinTheta * Math.cos(phi);
        double y = radius * sinTheta * Math.sin(phi);
        double z = radius * Math.cos(theta);
        return new double[] { x, y, z };
    }

    /**
     * Creates a spherical grid with generated theta and phi vertices.
     *
     * @param radius sphere radius
     * @param numThetaCells number of theta cells
     * @param numPhiCells number of phi cells
     * @param thetaSpacing theta spacing option
     * @return the generated spherical grid
     */
    public static SphericalGrid generated(double radius, int numThetaCells,
            int numPhiCells, ThetaSpacing thetaSpacing) {

        if (numThetaCells < 1) {
            throw new IllegalArgumentException("numThetaCells must be at least 1.");
        }
        if (numPhiCells < 1) {
            throw new IllegalArgumentException("numPhiCells must be at least 1.");
        }
        if (thetaSpacing == null) {
            thetaSpacing = ThetaSpacing.UNIFORM_THETA;
        }

        double[] theta = new double[numThetaCells + 1];
        double[] phi = new double[numPhiCells + 1];

        switch (thetaSpacing) {
        case UNIFORM_THETA -> {
            for (int i = 0; i <= numThetaCells; i++) {
                theta[i] = Math.PI * i / numThetaCells;
            }
        }
        case UNIFORM_COS_THETA -> {
            for (int i = 0; i <= numThetaCells; i++) {
                double mu = 1.0 - 2.0 * i / numThetaCells;
                theta[i] = Math.acos(mu);
            }
        }
        default -> throw new IllegalStateException("Unexpected theta spacing: " + thetaSpacing);
        }

        for (int i = 0; i <= numPhiCells; i++) {
            phi[i] = -Math.PI + 2.0 * Math.PI * i / numPhiCells;
        }

        return new SphericalGrid(theta, phi, radius);
    }

    /**
     * Creates a short text summary.
     *
     * @return a summary string
     */
    public String summary() {
        return String.format(
                "Spherical GSM grid: radius=%g, theta cells=%d, phi cells=%d",
                radius, getNumThetaCells(), getNumPhiCells());
    }

    /**
     * Normalizes an azimuthal angle to the range [-pi, pi].
     *
     * @param phi the angle to normalize
     * @return the normalized angle
     */
    public static double normalizePhi(double phi) {
        double result = phi;
        while (result < -Math.PI) {
            result += 2.0 * Math.PI;
        }
        while (result > Math.PI) {
            result -= 2.0 * Math.PI;
        }
        return result;
    }

    /**
     * Null check helper.
     *
     * @param other object to check
     * @return the same object
     */
    private static SphericalGrid requireNonNull(SphericalGrid other) {
        if (other == null) {
            throw new IllegalArgumentException("Cannot copy a null SphericalGrid.");
        }
        return other;
    }
    
    /**
 	 * Creates a list of human-readable feedback strings for the Cartesian grid.
 	 *
 	 * @return a list of feedback strings
 	 */
     public ArrayList<String> feedbackStrings() {
	 	ArrayList<String> feedback = new ArrayList<>();
	 	feedback.add(String.format("Radius: %g", radius));
	 	feedback.add(String.format("Theta: %d cells from %g to %g",
	 			getNumThetaCells(), Math.toDegrees(thetaGrid.min()), Math.toDegrees(thetaGrid.max())));
	 	feedback.add(String.format("Phi: %d cells from %g to %g",
	 			getNumPhiCells(), Math.toDegrees(phiGrid.min()), Math.toDegrees(phiGrid.max())));
	 	return feedback;
     }
}