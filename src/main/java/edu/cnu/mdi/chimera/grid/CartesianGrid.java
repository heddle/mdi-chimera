package edu.cnu.mdi.chimera.grid;

import java.awt.geom.Point2D;
import java.util.ArrayList;

/**
 * Cartesian grid in GSM coordinates.
 * <p>
 * Coordinates are interpreted as GSM coordinates measured in the length unit
 * specified by the owning {@code ChimeraGridSpec}. For the default 
 * application, that unit is Earth radii.
 * </p>
 */
public final class CartesianGrid {

    /** Grid vertices along GSM x. */
    private final Grid1D xGrid;

    /** Grid vertices along GSM y. */
    private final Grid1D yGrid;

    /** Grid vertices along GSM z. */
    private final Grid1D zGrid;

    /** Optional x offset of the grid origin relative to Earth center. */
    private final double xOffset;

    /** Optional y offset of the grid origin relative to Earth center. */
    private final double yOffset;

    /** Optional z offset of the grid origin relative to Earth center. */
    private final double zOffset;

    /**
     * Creates a Cartesian grid.
     *
     * @param xGrid the x grid
     * @param yGrid the y grid
     * @param zGrid the z grid
     * @param xOffset the x offset
     * @param yOffset the y offset
     * @param zOffset the z offset
     */
    public CartesianGrid(Grid1D xGrid, Grid1D yGrid, Grid1D zGrid,
            double xOffset, double yOffset, double zOffset) {

        if (xGrid == null || yGrid == null || zGrid == null) {
            throw new IllegalArgumentException("Cartesian grid axes must not be null.");
        }

        validateFinite("xOffset", xOffset);
        validateFinite("yOffset", yOffset);
        validateFinite("zOffset", zOffset);

        this.xGrid = new Grid1D(xGrid);
        this.yGrid = new Grid1D(yGrid);
        this.zGrid = new Grid1D(zGrid);
        this.xOffset = xOffset;
        this.yOffset = yOffset;
        this.zOffset = zOffset;
    }

    /**
     * Creates a Cartesian grid from coordinate arrays.
     *
     * @param xArray x grid vertices
     * @param yArray y grid vertices
     * @param zArray z grid vertices
     * @param xOffset x offset
     * @param yOffset y offset
     * @param zOffset z offset
     */
    public CartesianGrid(double[] xArray, double[] yArray, double[] zArray,
            double xOffset, double yOffset, double zOffset) {
        this(new Grid1D(xArray), new Grid1D(yArray), new Grid1D(zArray),
                xOffset, yOffset, zOffset);
    }

    /**
     * Copy constructor.
     *
     * @param other the grid to copy
     */
    public CartesianGrid(CartesianGrid other) {
        this(requireNonNull(other).xGrid, other.yGrid, other.zGrid,
                other.xOffset, other.yOffset, other.zOffset);
    }

    /**
     * Gets the x grid.
     *
     * @return the x grid
     */
    public Grid1D getXGrid() {
        return new Grid1D(xGrid);
    }

    /**
     * Gets the y grid.
     *
     * @return the y grid
     */
    public Grid1D getYGrid() {
        return new Grid1D(yGrid);
    }

    /**
     * Gets the z grid.
     *
     * @return the z grid
     */
    public Grid1D getZGrid() {
        return new Grid1D(zGrid);
    }

    /**
     * Gets the x offset.
     *
     * @return the x offset
     */
    public double getXOffset() {
        return xOffset;
    }

    /**
     * Gets the y offset.
     *
     * @return the y offset
     */
    public double getYOffset() {
        return yOffset;
    }

    /**
     * Gets the z offset.
     *
     * @return the z offset
     */
    public double getZOffset() {
        return zOffset;
    }

    /**
     * Gets the minimum global x value.
     *
     * @return the minimum global x value
     */
    public double getXMin() {
        return xGrid.min() + xOffset;
    }

    /**
     * Gets the maximum global x value.
     *
     * @return the maximum global x value
     */
    public double getXMax() {
        return xGrid.max() + xOffset;
    }

    /**
     * Gets the minimum global y value.
     *
     * @return the minimum global y value
     */
    public double getYMin() {
        return yGrid.min() + yOffset;
    }

    /**
     * Gets the maximum global y value.
     *
     * @return the maximum global y value
     */
    public double getYMax() {
        return yGrid.max() + yOffset;
    }

    /**
     * Gets the minimum global z value.
     *
     * @return the minimum global z value
     */
    public double getZMin() {
        return zGrid.min() + zOffset;
    }

    /**
     * Gets the maximum global z value.
     *
     * @return the maximum global z value
     */
    public double getZMax() {
        return zGrid.max() + zOffset;
    }

    /**
     * Gets the number of x cells.
     *
     * @return the number of x cells
     */
    public int getNumXCells() {
        return xGrid.numCells();
    }

    /**
     * Gets the number of y cells.
     *
     * @return the number of y cells
     */
    public int getNumYCells() {
        return yGrid.numCells();
    }

    /**
     * Gets the number of z cells.
     *
     * @return the number of z cells
     */
    public int getNumZCells() {
        return zGrid.numCells();
    }

    /**
     * Locates a point in the Cartesian grid.
     *
     * @param x the global x coordinate
     * @param y the global y coordinate
     * @param z the global z coordinate
     * @param indices output array of length at least three
     * @return the supplied indices array, filled with x, y, z cell indices
     */
    public int[] getIndices(double x, double y, double z, int[] indices) {
        if (indices == null || indices.length < 3) {
            throw new IllegalArgumentException("indices must have length at least 3.");
        }

        indices[0] = xGrid.locateInterval(x - xOffset);
        indices[1] = yGrid.locateInterval(y - yOffset);
        indices[2] = zGrid.locateInterval(z - zOffset);
        return indices;
    }

    /**
     * Creates a uniform Cartesian grid.
     *
     * @param xmin minimum x
     * @param xmax maximum x
     * @param nx number of x cells
     * @param ymin minimum y
     * @param ymax maximum y
     * @param ny number of y cells
     * @param zmin minimum z
     * @param zmax maximum z
     * @param nz number of z cells
     * @return the generated Cartesian grid
     */
    public static CartesianGrid uniform(double xmin, double xmax, int nx,
            double ymin, double ymax, int ny,
            double zmin, double zmax, int nz) {
        return new CartesianGrid(
                Grid1D.uniform(xmin, xmax, nx),
                Grid1D.uniform(ymin, ymax, ny),
                Grid1D.uniform(zmin, zmax, nz),
                0.0, 0.0, 0.0);
    }

    /**
     * Creates a short text summary.
     *
     * @return a summary string
     */
    public String summary() {
        return String.format(
                "Cartesian GSM grid: X[%g, %g] %d cells, Y[%g, %g] %d cells, Z[%g, %g] %d cells, offsets=(%g,%g,%g)",
                getXMin(), getXMax(), getNumXCells(),
                getYMin(), getYMax(), getNumYCells(),
                getZMin(), getZMax(), getNumZCells(),
                xOffset, yOffset, zOffset);
    }
    
    /**
	 * Creates a list of human-readable feedback strings for the Cartesian grid.
	 *
	 * @return a list of feedback strings
	 */
    public ArrayList<String> feedbackStrings() {
    	ArrayList<String> feedback = new ArrayList<>();
		feedback.add(String.format("X: %d cells from %g to %g", getNumXCells(), getXMin(), getXMax()));
		feedback.add(String.format("Y: %d cells from %g to %g", getNumYCells(), getYMin(), getYMax()));
		feedback.add(String.format("Z: %d cells from %g to %g", getNumZCells(), getZMin(), getZMax()));
		return feedback;
    }

    /**
     * Validates a finite double.
     *
     * @param name field name
     * @param value field value
     */
    private static void validateFinite(String name, double value) {
        if (!Double.isFinite(value)) {
            throw new IllegalArgumentException(name + " must be finite.");
        }
    }

    /**
     * Null check helper.
     *
     * @param other object to check
     * @return the same object
     */
    private static CartesianGrid requireNonNull(CartesianGrid other) {
        if (other == null) {
            throw new IllegalArgumentException("Cannot copy a null CartesianGrid.");
        }
        return other;
    }
}