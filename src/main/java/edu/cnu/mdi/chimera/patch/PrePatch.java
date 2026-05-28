package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.List;

import edu.cnu.mdi.chimera.app.ChimeraApp;
import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.cell.IntersectionType;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.Crossing;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.edge.Edge;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.GridSupport;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A prepatch is the region on the sphere surface enclosed by the intersection
 * of a single Cartesian cell with the sphere (paper §5.2). It is identified
 * by the Cartesian 3-tuple {@code (nx, ny, nz)} and carries no spherical-grid
 * index yet — those are assigned during the theta and phi splice steps.
 *
 * <p>Its boundary is a closed loop of {@link GeneralCurve}s, one per
 * intersecting face of the cell. Each curve runs from one edge–sphere
 * intersection point to the next, staying on the sphere surface and on
 * exactly one cell face.</p>
 *
 * <h2>Factory method</h2>
 * <p>Use {@link #from(Cell, CartesianGrid, SphericalGrid)} rather than
 * constructing directly. That method:</p>
 * <ol>
 *   <li>Retrieves the ordered edge list from the cell.</li>
 *   <li>For each consecutive edge pair, finds their common face.</li>
 *   <li>Fetches three corners of that face to define the plane.</li>
 *   <li>Constructs a {@link GeneralCurve} from the first intersection
 *       point to the second.</li>
 *   <li>Passes the resulting curve list to the {@code PrePatch} constructor,
 *       which validates the closed loop.</li>
 * </ol>
 *
 * <h2>Kiss cells</h2>
 * <p>Kiss cells have no edge intersections and are handled separately; this
 * class throws {@link IllegalArgumentException} if passed one.</p>
 */
public class PrePatch extends BasePatch {

    // -----------------------------------------------------------------------
    // Construction — private; use the factory method
    // -----------------------------------------------------------------------

    /**
     * Constructs a PrePatch from an already-built curve list.
     *
     * @param curves   ordered closed loop of {@link GeneralCurve}s
     * @param nx       Cartesian x cell index
     * @param ny       Cartesian y cell index
     * @param nz       Cartesian z cell index
     */
    private PrePatch(List<BaseCurve> curves,
			int nx, int ny, int nz) {
		super(curves, nx, ny, nz);

		// sanity check, all curves should be general curves
		for (BaseCurve curve : curves) {
			if (!(curve instanceof GeneralCurve)) {
				throw new IllegalArgumentException("PrePatch constructor: all curves must be GeneralCurves.");
			}
		}
	}

    // -----------------------------------------------------------------------
    // Factory method
    // -----------------------------------------------------------------------

    /**
     * Builds a {@link PrePatch} from a classified intersection {@link Cell}.
     *
     * <p>The cell must be a genuine edge-crossing intersection cell — not a
     * Kiss cell and not fully inside or outside the sphere.</p>
     *
     * @param cell     the intersection cell (must not be Kiss or UNKNOWN)
      * @return the constructed PrePatch
     * @throws IllegalArgumentException if the cell is a Kiss, has no ordered
     *                                  edges, has an invalid common-face
     *                                  result, or if the curves do not form
     *                                  a closed loop
     */
    public static PrePatch from(Cell cell) {
    	
        if (cell.getIntersectionType() == IntersectionType.KISS) {
            throw new IllegalArgumentException(
                "PrePatch.from: Kiss cells must be handled separately.");
        }

        Edge[] edges = cell.getOrderedEdges();
        if (edges == null || edges.length == 0) {
            throw new IllegalArgumentException(
                "PrePatch.from: cell (" + cell.nx + "," + cell.ny + "," + cell.nz +
                ") has no ordered edges.");
        }
        
        CartesianGrid cartGrid = ChimeraApp.getInstance().getCartesianGrid();

        // Pre-fetch all eight cell corners once — used for every face lookup.
        double[][] cellCorners = GridSupport.getCellCorners(
                cartGrid, cell.nx, cell.ny, cell.nz);

        int n = edges.length;
        List<BaseCurve> curves = new ArrayList<>(n);

        for (int i = 0; i < n; i++) {
            Edge current = edges[i];
            Edge next    = edges[(i + 1) % n];

            // The curve on this face runs from the current edge's intersection
            // point to the next edge's intersection point.
            Point3D.Double p0 = current.getIntersection();
            Point3D.Double p1 = next.getIntersection();

            // Find the face shared by the two edges — the face this curve lies on.
            int faceIndex = current.getCommonFace(next);
            if (faceIndex < 0) {
                throw new IllegalArgumentException(String.format(
                    "PrePatch.from: edges %d and %d of cell (%d,%d,%d) share no face.",
                    i, (i + 1) % n, cell.nx, cell.ny, cell.nz));
            }

            // Get three corners of that face to define its plane.
            Point3D.Double[] faceCorners = getFaceCorners(cellCorners, faceIndex);

            GeneralCurve curve = new GeneralCurve(
					p0, p1, ChimeraApp.getInstance().getRadius(),
					faceCorners[0], faceCorners[1], faceCorners[2]);
            curve.getThetaCrossings(); // precompute theta crossings for diagnostics
            curves.add(curve);
        }

        return new PrePatch(curves, cell.nx, cell.ny, cell.nz);
    }

    // -----------------------------------------------------------------------
    // BasePatch abstract implementation
    // -----------------------------------------------------------------------

    /**
     * Tests whether a sphere-surface point lies within this prepatch.
     *
     * <p>The point is located in the Cartesian grid (which Cartesian cell does
     * it belong to?) and in the spherical grid (which theta/phi cell?). It is
     * inside this prepatch if and only if its Cartesian cell indices match
     * {@code (nx, ny, nz)}.</p>
     *
     * @param x Cartesian x coordinate (on the sphere surface)
     * @param y Cartesian y coordinate
     * @param z Cartesian z coordinate
     * @return {@code true} if the point belongs to this prepatch's Cartesian cell
     */
    @Override
    public boolean containsPoint(double x, double y, double z) {
    	CartesianGrid cartesianGrid = ChimeraApp.getInstance().getCartesianGrid();
        int[] idx = cartesianGrid.getIndices(x, y, z, new int[3]);
        return idx[0] == nx && idx[1] == ny && idx[2] == nz;
    }

    // -----------------------------------------------------------------------
    // Static helpers
    // -----------------------------------------------------------------------

    /**
     * Extracts three {@link Point3D.Double} corners of the given face from the
     * full 8-corner cell array.
     *
     * <p>Only three corners are needed to define the face plane; we take the
     * first three of the four face corners in canonical order.</p>
     *
     * @param cellCorners 8-element array of (x,y,z) cell corners
     * @param faceIndex   face index [0, 5]
     * @return array of three {@link Point3D.Double} corners
     */
    private static Point3D.Double[] getFaceCorners(double[][] cellCorners,
                                                   int faceIndex) {
        int[] ci = GridSupport.getFaceCornerIndices(faceIndex);
        // ci has 4 elements; we only need 3 for the plane.
        return new Point3D.Double[]{
            new Point3D.Double(cellCorners[ci[0]]),
            new Point3D.Double(cellCorners[ci[1]]),
            new Point3D.Double(cellCorners[ci[2]])
        };
    }

    // -----------------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------------

    /**
     * Prints a human-readable summary of this prepatch.
     */
    public void printSummary() {
        System.out.printf("PrePatch [%d,%d,%d]  curves=%d  polar=N%s S%s%n",
                nx, ny, nz, curves.size(),
                enclosesNorthPole() ? "✓" : "✗",
                enclosesSouthPole() ? "✓" : "✗");
        for (int i = 0; i < curves.size(); i++) {
            System.out.printf("  curve[%d] %s%n", i, curves.get(i));
        }
        System.out.printf("  area (n=5):  %.6e%n", areaEstimate(5));
        System.out.printf("  perimeter:   %.6f%n", perimeter());
    }

    /**
     * Builds a PrePatch list from all non-Kiss intersecting cells and prints
     * a one-line diagnostic summary.
     *
     * @param cells    the full list of intersecting cells from the scanner
     * @param cartGrid the Cartesian grid
     * @param sphGrid  the spherical grid
     * @return list of PrePatches, one per non-Kiss cell
     */
    public static List<PrePatch> buildAll(List<Cell> cells) {
        List<PrePatch> patches = new ArrayList<>();
        int failed = 0;

        for (Cell cell : cells) {
            if (cell.getIntersectionType() == IntersectionType.KISS) {
                continue;
            }
            try {
                patches.add(PrePatch.from(cell));
            } catch (IllegalArgumentException ex) {
                failed++;
                System.err.printf("[PrePatch] Failed for cell (%d,%d,%d): %s%n",
                        cell.nx, cell.ny, cell.nz, ex.getMessage());
            }
        }

        System.out.printf("PrePatch build: %d succeeded, %d failed%n",
                patches.size(), failed);
        return patches;
    }
    
    /**
	 * Counts the total number of theta crossings across all curves in this
	 * prepatch. This is a diagnostic method.
	 *
	 * @return total number of theta crossings in this prepatch
	 */
    public int getTotalNumberThetaCrossings() {
		int count = 0;
		for (BaseCurve curve : curves) {
			if (curve instanceof GeneralCurve) {
				GeneralCurve gc = (GeneralCurve) curve;
				count += gc.thetaCrossings.size();
			}
		}
		return count;
	}
    
    /**
	 * Retrieves a list of all theta crossings across all curves in this
	 * prepatch. 
	 * @return list of all theta crossings in this prepatch
	 */
	public List<Crossing> getAllThetaCrossings() {
		List<Crossing> crossings = new ArrayList<>();
		for (BaseCurve curve : curves) {
			if (curve instanceof GeneralCurve) {
				GeneralCurve gc = (GeneralCurve) curve;
				crossings.addAll(gc.thetaCrossings);
			}
		}
		return Crossing.removeDuplicates(crossings);
	}

}