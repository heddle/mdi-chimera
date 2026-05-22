package edu.cnu.mdi.chimera.grid;

import java.util.Arrays;

/**
 * Immutable one-dimensional grid.
 * <p>
 * A {@code Grid1D} stores an ordered set of grid vertices. Cells are the
 * intervals between adjacent vertices. Thus a grid with {@code N} points has
 * {@code N - 1} cells.
 * </p>
 */
public final class Grid1D {

    /** Grid vertex coordinates in strictly increasing order. */
    private final double[] points;

    /** Maximum spacing between adjacent grid vertices. */
    private final double maxSpacing;

    /**
     * Creates a one-dimensional grid.
     * <p>
     * The input array is copied and sorted. After sorting, the values must be
     * finite and strictly increasing.
     * </p>
     *
     * @param inputPoints the grid vertex coordinates
     * @throws IllegalArgumentException if fewer than two points are supplied, if
     *         any value is not finite, or if duplicate points are present
     */
    public Grid1D(double[] inputPoints) {
        if (inputPoints == null || inputPoints.length < 2) {
            throw new IllegalArgumentException("A Grid1D requires at least two points.");
        }

        points = Arrays.copyOf(inputPoints, inputPoints.length);
        Arrays.sort(points);

        for (int i = 0; i < points.length; i++) {
            if (!Double.isFinite(points[i])) {
                throw new IllegalArgumentException("Grid point " + i + " is not finite.");
            }
            if (i > 0 && points[i] <= points[i - 1]) {
                throw new IllegalArgumentException("Grid points must be strictly increasing.");
            }
        }

        maxSpacing = computeMaxSpacing();
    }

    /**
     * Copy constructor.
     *
     * @param other the grid to copy
     * @throws IllegalArgumentException if {@code other} is {@code null}
     */
    public Grid1D(Grid1D other) {
        this(requireNonNullPoints(other));
    }

    /**
     * Gets the number of grid vertices.
     *
     * @return the number of grid vertices
     */
    public int numPoints() {
        return points.length;
    }

    /**
     * Gets the number of cells.
     *
     * @return {@code numPoints() - 1}
     */
    public int numCells() {
        return points.length - 1;
    }

    /**
     * Gets the grid value at a vertex index.
     *
     * @param index the vertex index
     * @return the coordinate value
     * @throws IndexOutOfBoundsException if the index is invalid
     */
    public double valueAt(int index) {
        if (index < 0 || index >= points.length) {
            throw new IndexOutOfBoundsException("Grid index out of range: " + index);
        }
        return points[index];
    }

    /**
     * Gets the minimum grid coordinate.
     *
     * @return the minimum coordinate
     */
    public double min() {
        return points[0];
    }

    /**
     * Gets the maximum grid coordinate.
     *
     * @return the maximum coordinate
     */
    public double max() {
        return points[points.length - 1];
    }

    /**
     * Gets a defensive copy of the grid points.
     *
     * @return a copy of the grid vertices
     */
    public double[] getPoints() {
        return Arrays.copyOf(points, points.length);
    }
    
    /**
	 * Checks if a value is within a tolerance of any grid vertex.
	 *
	 * @param value the coordinate value to check
	 * @param tol the tolerance for matching vertices
	 * @return {@code true} if the value is within {@code tol} of any vertex, {@code false} otherwise
	 */
    public boolean valueIsAVertex(double value, double tol) {
    			if (!Double.isFinite(value)) {
			return false;
		}

		for (double point : points) {
			if (Math.abs(value - point) <= tol) {
				return true;
			}
		}

		return false;
	}

    /**
     * Locates the interval containing a value.
     * <p>
     * The return value is the index {@code i} such that
     * {@code points[i] <= value < points[i + 1]}. If {@code value} is exactly the
     * final grid point, the final interval is returned.
     * </p>
     *
     * @param value the coordinate value to locate
     * @return the interval index, or {@code -1} if the value is out of range
     */
    public int locateInterval(double value) {
        if (!Double.isFinite(value) || value < points[0] || value > points[points.length - 1]) {
            return -1;
        }

        int index = Arrays.binarySearch(points, value);

        if (index >= 0) {
            return Math.min(index, points.length - 2);
        }

        int insertionPoint = -index - 1;
        if (insertionPoint == 0 || insertionPoint >= points.length + 1) {
            return -1;
        }

        return insertionPoint - 1;
    }
    
    
    /**
     * Locates the cell containing a coordinate value.
     * <p>
     * A grid with {@code N} points has {@code N - 1} cells, indexed from
     * {@code 0} through {@code numCells() - 1}. This method returns the cell index
     * {@code i} such that:
     * </p>
     *
     * <pre>
     * points[i] <= value < points[i + 1]
     * </pre>
     *
     * <p>
     * If {@code value} is exactly equal to the final grid point, the last cell is
     * returned. This is useful for closed-domain feedback and patch indexing, where
     * the upper boundary should still belong to the final cell rather than being
     * reported as outside the grid.
     * </p>
     *
     * @param value the coordinate value
     * @return the cell index in {@code [0, numCells() - 1]}, or {@code -1} if the
     *         value is outside the grid
     */
    public int cellIndex(double value) {
        return locateInterval(value);
    }

    /**
     * Gets the index of the grid vertex closest to a value.
     *
     * @param value the coordinate value
     * @return the closest vertex index, or {@code -1} if the value is outside
     *         the grid range
     */
    public int closestIndex(double value) {
        int interval = locateInterval(value);
        if (interval < 0) {
            return -1;
        }

        if (interval == points.length - 1) {
            return interval;
        }

        double d0 = Math.abs(value - points[interval]);
        double d1 = Math.abs(value - points[interval + 1]);
        return (d1 < d0) ? interval + 1 : interval;
    }

    /**
     * Gets the average spacing between adjacent vertices.
     *
     * @return the average grid spacing
     */
    public double getAverageSpacing() {
        return (max() - min()) / numCells();
    }

    /**
     * Gets the maximum spacing between adjacent vertices.
     *
     * @return the maximum grid spacing
     */
    public double getMaxSpacing() {
        return maxSpacing;
    }

    /**
     * Gets conservative interval limits for a sphere centered on the origin.
     * <p>
     * This is a bulk filter used to avoid checking grid cells that cannot
     * possibly intersect a sphere of the given radius.
     * </p>
     *
     * @param radius the sphere radius
     * @return a two-element array containing lower and upper cell-index limits
     */
    public int[] bulkFilterLimits(double radius) {
        if (radius <= 0.0 || !Double.isFinite(radius)) {
            throw new IllegalArgumentException("Radius must be positive and finite.");
        }

        int lower = locateInterval(-radius) - 1;
        int upper = locateInterval(radius) + 1;

        lower = Math.max(0, lower);
        upper = Math.min(numCells() - 1, upper);

        return new int[] { lower, upper };
    }

    /**
     * Creates a uniformly spaced grid.
     *
     * @param min the minimum coordinate
     * @param max the maximum coordinate
     * @param numCells the number of cells
     * @return the generated grid
     */
    public static Grid1D uniform(double min, double max, int numCells) {
        if (!Double.isFinite(min) || !Double.isFinite(max) || max <= min) {
            throw new IllegalArgumentException("Uniform grid requires finite min < max.");
        }
        if (numCells < 1) {
            throw new IllegalArgumentException("Uniform grid requires at least one cell.");
        }

        double[] values = new double[numCells + 1];
        double step = (max - min) / numCells;

        for (int i = 0; i <= numCells; i++) {
            values[i] = min + i * step;
        }

        return new Grid1D(values);
    }

    /**
     * Computes the maximum adjacent spacing.
     *
     * @return the maximum spacing
     */
    private double computeMaxSpacing() {
        double max = 0.0;
        for (int i = 1; i < points.length; i++) {
            max = Math.max(max, points[i] - points[i - 1]);
        }
        return max;
    }

    /**
     * Helper for copy construction.
     *
     * @param other the grid to check
     * @return the grid points from {@code other}
     */
    private static double[] requireNonNullPoints(Grid1D other) {
        if (other == null) {
            throw new IllegalArgumentException("Cannot copy a null Grid1D.");
        }
        return other.points;
    }
}