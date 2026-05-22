package edu.cnu.mdi.chimera.model;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;

/**
 * Immutable specification of a compositegrid.
 * <p>
 * A grid specification contains the Cartesian grid, the spherical
 * surface grid, and the physical interpretation of the coordinates. It does
 * not contain computed prepatches, theta patches, final patches, or area and
 * perimeter results. Those should be stored in a separate result object after
 * the algorithm runs.
 * </p>
 */
public final class ChimeraGridSpec {

    /** Coordinate system. */
    private final CoordinateSystem coordinateSystem;

    /** Length unit. */
    private final LengthUnit lengthUnit;

    /** Cartesian grid. */
    private final CartesianGrid cartesianGrid;

    /** Spherical grid. */
    private final SphericalGrid sphericalGrid;

    /** Optional name for display and logging. */
    private final String name;

    /**
     * Creates a grid specification.
     *
     * @param name display name
     * @param coordinateSystem coordinate system
     * @param lengthUnit length unit
     * @param cartesianGrid Cartesian grid
     * @param sphericalGrid spherical grid
     */
    public ChimeraGridSpec(String name, CoordinateSystem coordinateSystem,
            LengthUnit lengthUnit, CartesianGrid cartesianGrid,
            SphericalGrid sphericalGrid) {

        if (cartesianGrid == null) {
            throw new IllegalArgumentException("Cartesian grid must not be null.");
        }
        if (sphericalGrid == null) {
            throw new IllegalArgumentException("Spherical grid must not be null.");
        }

        this.name = (name == null || name.isBlank()) ? "Untitled Grid" : name.trim();
        this.coordinateSystem = coordinateSystem == null ? CoordinateSystem.GSM : coordinateSystem;
        this.lengthUnit = lengthUnit == null ? LengthUnit.EARTH_RADII : lengthUnit;
        this.cartesianGrid = new CartesianGrid(cartesianGrid);
        this.sphericalGrid = new SphericalGrid(sphericalGrid);

        validateSphereEnclosed();
    }

    /**
     * Gets the display name.
     *
     * @return the display name
     */
    public String getName() {
        return name;
    }

    /**
     * Gets the coordinate system.
     *
     * @return the coordinate system
     */
    public CoordinateSystem getCoordinateSystem() {
        return coordinateSystem;
    }

    /**
     * Gets the length unit.
     *
     * @return the length unit
     */
    public LengthUnit getLengthUnit() {
        return lengthUnit;
    }

    /**
     * Gets a copy of the Cartesian grid.
     *
     * @return the Cartesian grid
     */
    public CartesianGrid getCartesianGrid() {
        return new CartesianGrid(cartesianGrid);
    }

    /**
     * Gets a copy of the spherical grid.
     *
     * @return the spherical grid
     */
    public SphericalGrid getSphericalGrid() {
        return new SphericalGrid(sphericalGrid);
    }

    /**
     * Creates a human-readable summary suitable for logging.
     *
     * @return a summary string
     */
    public String summary() {
        String sep = System.lineSeparator();
        return "Grid: " + name + sep
                + "Coordinate system: " + coordinateSystem + sep
                + "Length unit: " + lengthUnit + sep
                + cartesianGrid.summary() + sep
                + sphericalGrid.summary();
    }

    /**
     * Checks that the spherical surface is enclosed by the Cartesian grid.
     * <p>
     * This is intentionally conservative and assumes the sphere is centered at
     * the global origin.
     * </p>
     */
    private void validateSphereEnclosed() {
        double r = sphericalGrid.getRadius();

        if (cartesianGrid.getXMin() > -r || cartesianGrid.getXMax() < r
                || cartesianGrid.getYMin() > -r || cartesianGrid.getYMax() < r
                || cartesianGrid.getZMin() > -r || cartesianGrid.getZMax() < r) {
            throw new IllegalArgumentException(
                    "The Cartesian grid must enclose the spherical grid of radius " + r + ".");
        }
    }
    
    /**
	 * Gets the Cartesian and spherical grid cell indices for a point on the
	 * spherical surface.
	 *
	 * @param theta the spherical theta coordinate
	 * @param phi the spherical phi coordinate
	 * @param radius the spherical radius
	 * @param indices an array of length 5 to hold the output indices; on return,
	 *        this contains [nx, ny, nz, ntheta, nphi]
	 */
    public void getPatchIndices(double theta, double phi, double radius, int[] indices) {
    	
    	double sinTheta = Math.sin(theta);
    	
    	double x = radius * sinTheta * Math.cos(phi);
    	double y = radius * sinTheta * Math.sin(phi);
    	double z = radius * Math.cos(theta);   	
    	
		indices[0] = cartesianGrid.getXGrid().cellIndex(x);
		indices[1] = cartesianGrid.getYGrid().cellIndex(y);
		indices[2] = cartesianGrid.getZGrid().cellIndex(z);
		indices[3] = sphericalGrid.getThetaGrid().cellIndex(theta);
		indices[4] = sphericalGrid.getPhiGrid().cellIndex(phi);
	}
}