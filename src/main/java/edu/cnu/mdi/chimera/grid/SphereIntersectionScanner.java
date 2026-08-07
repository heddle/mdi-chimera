package edu.cnu.mdi.chimera.grid;

import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.cell.IntersectionType;
import edu.cnu.mdi.chimera.cell.KissGeometry;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Scans a {@link CartesianGrid} and identifies all cubic cells that are
 * intersected by — or tangentially touched by — the surface of a sphere
 * centred at the origin.
 *
 * <h2>Main intersection pass</h2>
 * <ol>
 *   <li>Use {@link Grid1D#bulkFilterLimits(double)} on each axis to obtain
 *       conservative index ranges; cells outside those ranges are skipped.</li>
 *   <li>For each candidate cell, compute the corner-inside bitmask.</li>
 *   <li>Cells with all corners inside (0xFF) are skipped (no surface).</li>
 *   <li>Cells with at least one corner inside and one outside are genuine
 *       edge-crossing intersections; a {@link Cell} is constructed for each.</li>
 *   <li>All-outside cells (0x00) are forwarded to the Kiss test.</li>
 *   <li>After construction, polar cells are flagged.</li>
 * </ol>
 *
 * <h2>Kiss detection</h2>
 * <p>A Kiss occurs when the sphere's curvature causes it to penetrate a face
 * of a cell whose corners are all outside the sphere — a spherical cap poking
 * through without touching any corner. For each axis-aligned face the test is:
 * (1) the origin's projection onto the face plane must fall strictly inside
 * the face quad, and (2) the perpendicular distance from the origin to the
 * face plane must be strictly less than {@code radius}. Condition (1) ensures
 * the penetration is genuinely interior to this face and not something a
 * neighbouring cell's corner bits would already catch. No plane construction,
 * rotation, or tolerance is required.</p>
 *
 * <h2>Polar detection</h2>
 * <p>A cell is polar if its AABB brackets the north pole {@code (0,0,+r)} or
 * the south pole {@code (0,0,-r)}.</p>
 */
public class SphereIntersectionScanner {

    /** The Cartesian grid to scan. */
    private final CartesianGrid grid;

    /** Radius of the sphere centred at the origin. */
    private final double radius;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * Creates a scanner for the given grid and sphere radius.
     *
     * @param grid   the Cartesian grid; must not be {@code null}
     * @param radius sphere radius; must be positive and finite
     */
    public SphereIntersectionScanner(CartesianGrid grid, double radius) {
        if (grid == null) {
            throw new IllegalArgumentException("CartesianGrid must not be null.");
        }
        if (!Double.isFinite(radius) || radius <= 0.0) {
            throw new IllegalArgumentException("Radius must be positive and finite.");
        }
        this.grid   = grid;
        this.radius = radius;
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Runs the scan and returns all cells whose surface is intersected by or
     * tangent to the sphere. The list is unmodifiable; order is not guaranteed.
     *
     * @return an unmodifiable list of {@link Cell} objects
     */
    public List<Cell> scan() {
        List<Cell> result = new ArrayList<>();

        int[] xLimits = grid.getXGrid().bulkFilterLimits(radius);
        int[] yLimits = grid.getYGrid().bulkFilterLimits(radius);
        int[] zLimits = grid.getZGrid().bulkFilterLimits(radius);

        int xLo = xLimits[0], xHi = xLimits[1];
        int yLo = yLimits[0], yHi = yLimits[1];
        int zLo = zLimits[0], zHi = zLimits[1];

        double r2 = radius * radius;

        for (int iz = zLo; iz <= zHi; iz++) {
            for (int iy = yLo; iy <= yHi; iy++) {
                for (int ix = xLo; ix <= xHi; ix++) {

                    int cornerBits = computeCornerBits(ix, iy, iz, r2);

                    if (cornerBits == 0xFF) {
                        // All corners inside — no surface intersection.
                        continue;
                    }

                    if (cornerBits == 0x00) {
                        // All corners outside — only a Kiss is possible.
                        Cell kiss = tryKiss(ix, iy, iz);
                        if (kiss != null) {
                            flagPolarCell(kiss, ix, iy, iz);
                            result.add(kiss);
                        }
                        continue;
                    }

                    // Mixed: genuine edge-crossing intersection.
                    Cell cell = new Cell(grid, ix, iy, iz, radius);
                    flagPolarCell(cell, ix, iy, iz);
                    result.add(cell);
                }
            }
        }

        return Collections.unmodifiableList(result);
    }

    // -----------------------------------------------------------------------
    // Kiss detection
    // -----------------------------------------------------------------------

    /**
     * Tests whether a cell whose corners are all outside the sphere is
     * nevertheless tangent to (kissed by) the sphere surface.
     *
     * <p>For each of the six axis-aligned faces of the cell the closest point
     * on the face to the origin is computed by clamping:</p>
     * <ul>
     *   <li>The coordinate perpendicular to the face is fixed at the face's
     *       constant value.</li>
     *   <li>The two in-plane coordinates are each clamped to the cell's
     *       interval for that axis.</li>
     * </ul>
     *
     * <p>If the distance from the origin to this clamped point is within
     * {@link #kissTol} of {@code radius}, the cell touches the sphere and is
     * returned as a Kiss cell. We test all six faces (not just the nearest)
     * because on a coarse grid, more than one face of an all-outside cell can
     * be tangent simultaneously (corner Kiss).</p>
     *
     * @param ix lower x cell index
     * @param iy lower y cell index
     * @param iz lower z cell index
     * @return a Kiss {@link Cell} if the sphere is tangent to any face,
     *         or {@code null} if the cell is genuinely non-intersecting
     */
    private Cell tryKiss(int ix, int iy, int iz) {
        Grid1D xGrid = grid.getXGrid();
        Grid1D yGrid = grid.getYGrid();
        Grid1D zGrid = grid.getZGrid();

        double x0 = xGrid.valueAt(ix),  x1 = xGrid.valueAt(ix + 1);
        double y0 = yGrid.valueAt(iy),  y1 = yGrid.valueAt(iy + 1);
        double z0 = zGrid.valueAt(iz),  z1 = zGrid.valueAt(iz + 1);

        // For each face: (fixedAxis, fixedValue, clampAxis1Min, clampAxis1Max,
        //                  clampAxis2Min, clampAxis2Max)
        // Face numbering matches GridSupport:
        //   0: z=z0,  1: z=z1,  2: y=y0,  3: y=y1,  4: x=x0,  5: x=x1
        if (faceIsKiss(z0, x0, x1, y0, y1)   // face 0: z=z0, clamp x and y
         || faceIsKiss(z1, x0, x1, y0, y1)   // face 1: z=z1
         || faceIsKiss(y0, x0, x1, z0, z1)   // face 2: y=y0, clamp x and z
         || faceIsKiss(y1, x0, x1, z0, z1)   // face 3: y=y1
         || faceIsKiss(x0, y0, y1, z0, z1)   // face 4: x=x0, clamp y and z
         || faceIsKiss(x1, y0, y1, z0, z1))  // face 5: x=x1
        {
            Cell cell = new Cell(grid, ix, iy, iz, radius);
            // Override the type set by Cell's constructor (which saw no edge
            // crossings and left it UNKNOWN) with the correct Kiss type.
            
            cell.setIntersectionType(IntersectionType.KISS);
            
            
            KissGeometry geometry = new KissGeometry(cell, grid, radius);
            cell.setKissGeometry(geometry);
            return cell;
        }

        return null;
    }

    /**
     * Returns {@code true} if the sphere penetrates this axis-aligned face
     * without any cell corner being inside the sphere.
     *
     * <p>The face has one coordinate fixed at {@code fixedVal}; the other two
     * span {@code [u0, u1]} and {@code [v0, v1]}. The closest point on the
     * face to the origin is found by clamping the origin's coordinates to the
     * face bounds:</p>
     * <pre>
     *   closest = (clamp(0, u0, u1),  clamp(0, v0, v1),  fixedVal)
     * </pre>
     * <p>If the squared distance from the origin to that point is strictly less
     * than {@code radius²}, the sphere penetrates the face.</p>
     *
     * <p>Note: there is no requirement that the projection fall in the face
     * interior. The sphere can penetrate a face even when the closest point is
     * at a face edge or corner, because in an all-outside cell no neighbouring
     * corner bit catches that region.</p>
     *
     * @param fixedVal the fixed coordinate of the face (e.g. z0 for face 0)
     * @param u0       lower bound of the first in-plane coordinate
     * @param u1       upper bound of the first in-plane coordinate
     * @param v0       lower bound of the second in-plane coordinate
     * @param v1       upper bound of the second in-plane coordinate
     * @return {@code true} if the sphere penetrates the face
     */
    private boolean faceIsKiss(double fixedVal,
                                double u0, double u1,
                                double v0, double v1) {
        // The face plane must actually intersect the sphere: if |fixedVal| >= R
        // there is no intersection circle and no arc can exist, regardless of
        // how close the clamped point is.
        if (fixedVal * fixedVal >= radius * radius) return false;

        double u = clamp(0.0, u0, u1);
        double v = clamp(0.0, v0, v1);
        double dist2 = fixedVal * fixedVal + u * u + v * v;
        return dist2 < radius * radius;
    }

    /**
     * Clamps {@code value} to {@code [lo, hi]}.
     */
    private static double clamp(double value, double lo, double hi) {
        return Math.max(lo, Math.min(hi, value));
    }

    // -----------------------------------------------------------------------
    // Corner bitmask
    // -----------------------------------------------------------------------

    /**
     * Computes the corner-inside bitmask for the cell at {@code (ix, iy, iz)}.
     * Bit {@code k} is set when corner {@code k} lies strictly inside the sphere.
     */
    private int computeCornerBits(int ix, int iy, int iz, double r2) {
        Grid1D xGrid = grid.getXGrid();
        Grid1D yGrid = grid.getYGrid();
        Grid1D zGrid = grid.getZGrid();

        double x0 = xGrid.valueAt(ix),  x1 = xGrid.valueAt(ix + 1);
        double y0 = yGrid.valueAt(iy),  y1 = yGrid.valueAt(iy + 1);
        double z0 = zGrid.valueAt(iz),  z1 = zGrid.valueAt(iz + 1);

        // Canonical corner layout (matches GridSupport):
        //   corner 0: (x0,y0,z0)  corner 1: (x1,y0,z0)
        //   corner 2: (x0,y1,z0)  corner 3: (x1,y1,z0)
        //   corner 4: (x0,y0,z1)  corner 5: (x1,y0,z1)
        //   corner 6: (x0,y1,z1)  corner 7: (x1,y1,z1)
        double[] xs = { x0, x1, x0, x1, x0, x1, x0, x1 };
        double[] ys = { y0, y0, y1, y1, y0, y0, y1, y1 };
        double[] zs = { z0, z0, z0, z0, z1, z1, z1, z1 };

        int bits = 0;
        for (int k = 0; k < 8; k++) {
            double x = xs[k], y = ys[k], z = zs[k];
            if (x * x + y * y + z * z < r2) {
                bits |= GridSupport.CORNERBITS[k];
            }
        }
        return bits;
    }
    
    
    /**
     * Returns a list of all cells of the given {@link IntersectionType}.
     *
     * @param cells the full cell list; must not be {@code null}
     * @param type  the intersection type to filter by; must not be {@code null}
     * @return a new list containing only cells of the requested type,
     *         in the same order as the input list
     */
    public static List<Cell> getCellsOfType(List<Cell> cells, IntersectionType type) {
        if (cells == null) {
            throw new IllegalArgumentException("Cell list must not be null.");
        }
        if (type == null) {
            throw new IllegalArgumentException("IntersectionType must not be null.");
        }
        List<Cell> result = new ArrayList<>();
        for (Cell cell : cells) {
            if (cell.getIntersectionType() == type) {
                result.add(cell);
            }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------------

    /**
     * Prints a diagnostic summary of a list of {@link Cell} objects to
     * standard output, showing the count for each {@link IntersectionType}
     * and the total number of cells.
     *
     * <p>Only types with a non-zero count are printed. Example output:</p>
     * <pre>
     * Cell intersection summary:
     *   Corner In        :    142
     *   Face Cut         :   2881
     *   Kiss             :     16
     *   Skew Cut         :    563
     *   ----------------------------
     *   Total            :   3618
     * </pre>
     *
     * @param cells the list of cells to summarise; must not be {@code null}
     */
    public static void printDiagnostics(List<Cell> cells) {
        if (cells == null) {
            throw new IllegalArgumentException("Cell list must not be null.");
        }

        java.util.EnumMap<IntersectionType, Integer> counts =
                new java.util.EnumMap<>(IntersectionType.class);
        for (IntersectionType type : IntersectionType.values()) {
            counts.put(type, 0);
        }
        for (Cell cell : cells) {
            IntersectionType t = cell.getIntersectionType();
            counts.put(t, counts.get(t) + 1);
        }

        // Width of the longest type label, for column alignment.
        int maxLen = 0;
        for (IntersectionType type : IntersectionType.values()) {
            maxLen = Math.max(maxLen, type.toString().length());
        }
        String fmt  = "  %-" + maxLen + "s : %6d%n";
        String rule = "  " + "-".repeat(maxLen + 9);

        System.out.println("Cell intersection summary:");
        for (IntersectionType type : IntersectionType.values()) {
            int count = counts.get(type);
            if (count > 0) {
                System.out.printf(fmt, type.toString(), count);
            }
        }
        System.out.println(rule);
        System.out.printf(fmt, "Total", cells.size());
    }

    // -----------------------------------------------------------------------
    // Polar detection
    // -----------------------------------------------------------------------

    /**
     * Sets the {@code polar} flag on {@code cell} if the cell brackets the
     * north pole {@code (0, 0, +radius)} or the south pole {@code (0, 0, -radius)}.
     */
    private void flagPolarCell(Cell cell, int ix, int iy, int iz) {
        Grid1D xGrid = grid.getXGrid();
        Grid1D yGrid = grid.getYGrid();
        Grid1D zGrid = grid.getZGrid();

        double x0 = xGrid.valueAt(ix),  x1 = xGrid.valueAt(ix + 1);
        double y0 = yGrid.valueAt(iy),  y1 = yGrid.valueAt(iy + 1);
        double z0 = zGrid.valueAt(iz),  z1 = zGrid.valueAt(iz + 1);

        boolean bracketsXY = (x0 <= 0.0 && 0.0 < x1) && (y0 <= 0.0 && 0.0 < y1);
        if (bracketsXY) {
            boolean northPole = (z0 <= radius  && radius  <= z1);
            boolean southPole = (z0 <= -radius && -radius <= z1);
            if (northPole || southPole) {
                cell.setPolar(true);
            }
        }
    }
}
