package edu.cnu.mdi.chimera.grid;

import java.util.Arrays;

/**
 * Immutable one-dimensional grid.
 *
 * <p>A {@code Grid1D} stores an ordered sequence of grid vertices. Cells are the
 * closed-open intervals {@code [points[i], points[i+1])} between adjacent vertices,
 * except that the final interval is closed on both ends. A grid with {@code N}
 * vertices therefore has {@code N - 1} cells.</p>
 *
 * <p>All coordinate queries ({@link #locateInterval}, {@link #cellIndex},
 * {@link #closestIndex}, {@link #valueIsAVertex}) return {@code -1} when the
 * supplied value is outside the closed range {@code [min(), max()]}.</p>
 */
public final class Grid1D {

    /**
     * Absolute tolerance used by {@link #valueIsAVertex} when deciding whether a
     * floating-point value is close enough to a vertex to be treated as one.
     * The value {@code 1.0e-12} is well below any physically meaningful grid
     * spacing while still absorbing typical double-precision rounding error.
     */
    private static final double VERTEX_TOL = 1.0e-12;

    /** Grid vertex coordinates in strictly increasing order. */
    private final double[] points;

    /** Maximum spacing between adjacent grid vertices, cached at construction. */
    private final double maxSpacing;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    /**
     * Creates a one-dimensional grid from an array of vertex coordinates.
     *
     * <p>The input array is copied and sorted before validation. After sorting,
     * all values must be finite and strictly increasing (no duplicates).</p>
     *
     * @param inputPoints the grid vertex coordinates; must contain at least two
     *                    distinct finite values
     * @throws IllegalArgumentException if {@code inputPoints} is {@code null},
     *         contains fewer than two elements, contains a non-finite value, or
     *         contains duplicate coordinates
     */
    public Grid1D(double[] inputPoints) {
        if (inputPoints == null || inputPoints.length < 2) {
            throw new IllegalArgumentException("A Grid1D requires at least two points.");
        }

        points = Arrays.copyOf(inputPoints, inputPoints.length);
        Arrays.sort(points);

        for (int i = 0; i < points.length; i++) {
            if (!Double.isFinite(points[i])) {
                throw new IllegalArgumentException(
                        "Grid point " + i + " is not finite: " + points[i]);
            }
            if (i > 0 && points[i] <= points[i - 1]) {
                throw new IllegalArgumentException(
                        "Grid points must be strictly increasing; duplicate or "
                        + "inverted pair at index " + (i - 1) + " and " + i + ".");
            }
        }

        maxSpacing = computeMaxSpacing();
    }

    /**
     * Copy constructor. Creates a new {@code Grid1D} that is an independent copy
     * of {@code other}.
     *
     * @param other the grid to copy; must not be {@code null}
     * @throws IllegalArgumentException if {@code other} is {@code null}
     */
    public Grid1D(Grid1D other) {
        if (other == null) {
            throw new IllegalArgumentException("Cannot copy a null Grid1D.");
        }
        // other.points is already validated and sorted; copy directly.
        points = Arrays.copyOf(other.points, other.points.length);
        maxSpacing = other.maxSpacing;
    }

    // -------------------------------------------------------------------------
    // Dimensional queries
    // -------------------------------------------------------------------------

    /**
     * Returns the number of grid vertices.
     *
     * @return the vertex count; always &ge; 2
     */
    public int numPoints() {
        return points.length;
    }

    /**
     * Returns the number of cells.
     *
     * <p>Equivalent to {@code numPoints() - 1}.</p>
     *
     * @return the cell count; always &ge; 1
     */
    public int numCells() {
        return points.length - 1;
    }

    // -------------------------------------------------------------------------
    // Coordinate accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the coordinate of the vertex at the given index.
     *
     * @param index a zero-based vertex index in {@code [0, numPoints() - 1]}
     * @return the coordinate value at {@code index}
     * @throws IndexOutOfBoundsException if {@code index} is negative or &ge;
     *         {@link #numPoints()}
     */
    public double valueAt(int index) {
        if (index < 0 || index >= points.length) {
            throw new IndexOutOfBoundsException("Grid index out of range: " + index);
        }
        return points[index];
    }

    /**
     * Returns the minimum (leftmost) vertex coordinate.
     *
     * @return the minimum coordinate; equivalent to {@code valueAt(0)}
     */
    public double min() {
        return points[0];
    }

    /**
     * Returns the maximum (rightmost) vertex coordinate.
     *
     * @return the maximum coordinate; equivalent to {@code valueAt(numPoints() - 1)}
     */
    public double max() {
        return points[points.length - 1];
    }

    /**
     * Returns a defensive copy of the vertex coordinate array.
     *
     * @return a copy of the internal points array in strictly increasing order
     */
    public double[] getPoints() {
        return Arrays.copyOf(points, points.length);
    }

    // -------------------------------------------------------------------------
    // Vertex / interval queries
    // -------------------------------------------------------------------------

    /**
     * Determines whether a value coincides with a grid vertex, within the
     * tolerance {@value #VERTEX_TOL}.
     *
     * <p>The method first attempts an exact binary search. If that fails it
     * checks both immediate neighbors of the binary-search insertion point.
     * Only the nearest neighbor that satisfies
     * {@code |points[i] - value| <= VERTEX_TOL} is returned; if both neighbors
     * are within tolerance (which cannot occur on a valid grid because consecutive
     * vertices are strictly separated by much more than 2*VERTEX_TOL), 
     * the lower-index neighbor is preferred.</p>
     *
     * @param value the coordinate to test; non-finite values always return
     *              {@code -1}
     * @return the zero-based index of the matching vertex, or {@code -1} if no
     *         vertex is within {@value #VERTEX_TOL} of {@code value}
     */
    public int valueIsAVertex(double value) {
        if (!Double.isFinite(value)) {
            return -1;
        }

        int index = Arrays.binarySearch(points, value);
        if (index >= 0) {
            return index;
        }

        // binarySearch returns -(insertionPoint) - 1 on a miss.
        // insertionPoint is in [0, points.length].
        int ip = -index - 1;

        // Check the vertex immediately below the insertion point.
        if (ip > 0 && Math.abs(points[ip - 1] - value) <= VERTEX_TOL) {
            return ip - 1;
        }
        // Check the vertex at the insertion point itself.
        if (ip < points.length && Math.abs(points[ip] - value) <= VERTEX_TOL) {
            return ip;
        }

        return -1;
    }

    /**
     * Locates the grid interval that contains a coordinate value.
     *
     * <p>Returns the index {@code i} such that
     * {@code points[i] <= value < points[i + 1]}, except that when
     * {@code value} equals the final vertex exactly, the last interval index
     * ({@code numCells() - 1}) is returned so that the upper boundary is always
     * considered part of the grid.</p>
     *
     * @param value the coordinate value to locate
     * @return the interval index in {@code [0, numCells() - 1]}, or {@code -1}
     *         if {@code value} is non-finite or outside {@code [min(), max()]}
     */
    public int locateInterval(double value) {
        if (!Double.isFinite(value) || value < points[0] || value > points[points.length - 1]) {
            return -1;
        }

        int index = Arrays.binarySearch(points, value);
        if (index >= 0) {
            // Exact vertex hit: clamp to the last valid interval.
            return Math.min(index, points.length - 2);
        }

        // Miss: insertionPoint is in [1, points.length - 1] because the range
        // check above already excluded values below points[0] and above
        // points[points.length - 1].
        int insertionPoint = -index - 1;
        return insertionPoint - 1;
    }

    /**
     * Returns the cell index that contains a coordinate value.
     *
     * <p>A grid with {@code N} vertices has {@code N - 1} cells, indexed from
     * {@code 0} through {@code numCells() - 1}. The cell index {@code i}
     * satisfies:</p>
     * <pre>
     *   points[i] &lt;= value &lt; points[i + 1]
     * </pre>
     * <p>If {@code value} equals the final vertex exactly, the last cell is
     * returned. This ensures that closed-domain logic (e.g.&nbsp;patch selection
     * at {@code max()}) never reports an out-of-bounds cell.</p>
     *
     * @param value the coordinate value
     * @return the cell index in {@code [0, numCells() - 1]}, or {@code -1} if
     *         the value is outside the grid domain
     * @see #locateInterval(double)
     */
    public int cellIndex(double value) {
        return locateInterval(value);
    }

    /**
     * Returns the index of the grid vertex closest to a coordinate value.
     *
     * <p>When two vertices are equidistant, the lower-index vertex is returned.</p>
     *
     * @param value the coordinate value
     * @return the index of the closest vertex in {@code [0, numPoints() - 1]},
     *         or {@code -1} if the value is outside the grid domain
     */
    public int closestIndex(double value) {
        int interval = locateInterval(value);
        if (interval < 0) {
            return -1;
        }

        // locateInterval clamps to numCells()-1, so interval+1 is always valid.
        double d0 = Math.abs(value - points[interval]);
        double d1 = Math.abs(value - points[interval + 1]);
        return (d1 < d0) ? interval + 1 : interval;
    }

    // -------------------------------------------------------------------------
    // Spacing queries
    // -------------------------------------------------------------------------

    /**
     * Returns the average spacing between adjacent vertices.
     *
     * <p>Computed as {@code (max() - min()) / numCells()}.</p>
     *
     * @return the average vertex spacing; always positive
     */
    public double getAverageSpacing() {
        return (max() - min()) / numCells();
    }

    /**
     * Returns the maximum spacing between any two adjacent vertices.
     *
     * <p>This value is computed once at construction and cached.</p>
     *
     * @return the maximum adjacent spacing; always positive
     */
    public double getMaxSpacing() {
        return maxSpacing;
    }

    // -------------------------------------------------------------------------
    // Filtering utilities
    // -------------------------------------------------------------------------

    /**
     * Returns conservative cell-index bounds for a sphere centered on the origin.
     *
     * <p>This is a bulk filter intended to avoid testing cells that cannot
     * possibly intersect a sphere of the given radius. The returned bounds may
     * include cells that do not actually intersect the sphere (they are
     * conservative); they will never exclude a cell that does.</p>
     *
     * @param radius the sphere radius; must be positive and finite
     * @return a two-element array {@code {lower, upper}} where both values are
     * clamped to {@code [0, numCells() - 1]}, or {@code {-1, -1}} if
     * this grid does not overlap {@code [-radius, radius]}
     * @throws IllegalArgumentException if {@code radius} is not positive and finite
     */
    public int[] bulkFilterLimits(double radius) {
        if (!Double.isFinite(radius) || radius <= 0.0) {
            throw new IllegalArgumentException("Radius must be positive and finite; got " + radius + ".");
        }

        double lowerTarget = -radius;
        double upperTarget = radius;

        // If the entire sphere is outside the grid, return an empty/invalid range or throw
        if (upperTarget < min() || lowerTarget > max()) {
            return new int[] { -1, -1 }; // Or return new int[]{0, 0} depending on consumer contract
        }

        int rawLower = locateInterval(lowerTarget);
        int lower = (lowerTarget < min()) ? 0 : Math.max(0, rawLower - 1);

        int rawUpper = locateInterval(upperTarget);
        int upper = (upperTarget > max()) ? numCells() - 1 : Math.min(numCells() - 1, rawUpper + 1);

        return new int[] { lower, upper };
    }
    
    // -------------------------------------------------------------------------
    // Factory methods
    // -------------------------------------------------------------------------

    /**
     * Creates a uniformly spaced grid.
     *
     * <p>The {@code i}-th vertex is placed at
     * {@code min + i * (max - min) / numCells} for
     * {@code i} in {@code [0, numCells]}. The endpoints {@code min} and
     * {@code max} are set exactly, not computed by the recurrence, to avoid
     * floating-point drift at the boundaries.</p>
     *
     * @param min      the minimum coordinate; must be finite and strictly less
     *                 than {@code max}
     * @param max      the maximum coordinate; must be finite and strictly greater
     *                 than {@code min}
     * @param numCells the number of cells; must be &ge; 1
     * @return a new uniformly spaced {@code Grid1D}
     * @throws IllegalArgumentException if {@code min} or {@code max} is
     *         non-finite, if {@code min >= max}, or if {@code numCells < 1}
     */
    public static Grid1D uniform(double min, double max, int numCells) {
        if (!Double.isFinite(min) || !Double.isFinite(max) || max <= min) {
            throw new IllegalArgumentException(
                    "Uniform grid requires finite min < max; got min=" + min + ", max=" + max + ".");
        }
        if (numCells < 1) {
            throw new IllegalArgumentException(
                    "Uniform grid requires at least one cell; got " + numCells + ".");
        }

        double[] values = new double[numCells + 1];
        for (int i = 0; i <= numCells; i++) {
            // This safely guarantees values strictly increase from min to max
            values[i] = min * ((double)(numCells - i) / numCells) + max * ((double)i / numCells);
        }
        
        // Pin endpoints exactly to avoid floating-point drift.
        values[0] = min;
        values[numCells] = max;

        return new Grid1D(values);
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Computes the maximum spacing between adjacent vertices.
     *
     * <p>Called once from the constructor after {@link #points} has been
     * populated and validated.</p>
     *
     * @return the maximum adjacent spacing
     */
    private double computeMaxSpacing() {
        double max = 0.0;
        for (int i = 1; i < points.length; i++) {
            max = Math.max(max, points[i] - points[i - 1]);
        }
        return max;
    }
}