package edu.cnu.mdi.chimera.cell;

import edu.cnu.mdi.chimera.edge.Edge;
import edu.cnu.mdi.chimera.edge.EdgeOrdering;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.GridSupport;
import edu.cnu.mdi.util.Bits;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * One cubic cell of a {@link CartesianGrid}, annotated with the result of
 * intersecting it against a sphere centred at the origin.
 *
 * <p>
 * The cell is identified by its lower-corner grid indices {@code (nx, ny, nz)}.
 * During classification the eight corners are tested against the sphere radius;
 * the resulting inside/outside bitmask drives edge-intersection finding,
 * polygon construction, and {@link IntersectionType} assignment.
 * </p>
 */
public class Cell {

	// -----------------------------------------------------------------------
	// Grid position
	// -----------------------------------------------------------------------

	/** X cell index (lower corner). */
	public final int nx;

	/** Y cell index (lower corner). */
	public final int ny;

	/** Z cell index (lower corner). */
	public final int nz;

	// Only non-null for kisses
	private KissGeometry kissGeometry = null;

	// -----------------------------------------------------------------------
	// Sphere-intersection geometry
	// -----------------------------------------------------------------------

	/**
	 * Bitmask of corners that lie strictly inside the sphere. Bit {@code k} is set
	 * when corner {@code k} satisfies {@code r < radius}. Uses
	 * {@link GridSupport#CORNERBITS}.
	 */
	private int cornerBits;

	/**
	 * The eight corners of the cell in canonical order (see {@link GridSupport}).
	 * Each element is a {@code double[3]} holding (x, y, z).
	 */
	private double[][] corners;

	/**
	 * Edges that straddle the sphere surface, reordered into a closed loop by
	 * {@link EdgeOrdering}. {@code null} until computed.
	 */
	private Edge[] orderedEdges;

	/**
	 * The polygon formed by the sphere-surface intersection points, in the same
	 * cyclic order as {@code orderedEdges}. {@code null} until computed.
	 */
	private Point3D.Double[] intersectionPolygon;

	// -----------------------------------------------------------------------
	// Classification
	// -----------------------------------------------------------------------

	/** How the sphere intersects this cell. */
	private IntersectionType intersectionType = IntersectionType.UNKNOWN;

	/**
	 * Whether this cell contains a pole of the spherical grid (θ = 0 or θ = π).
	 * Polar cells require special treatment when mapping intersection polygons onto
	 * the spherical surface grid.
	 */
	private boolean polar = false;

	// -----------------------------------------------------------------------
	// Construction
	// -----------------------------------------------------------------------

	/**
	 * Creates a Cell and immediately classifies its intersection with the sphere.
	 *
	 * @param grid   the Cartesian grid that owns this cell
	 * @param nx     x cell index
	 * @param ny     y cell index
	 * @param nz     z cell index
	 * @param radius sphere radius (centred at the origin)
	 */
	public Cell(CartesianGrid grid, int nx, int ny, int nz, double radius) {
		this.nx = nx;
		this.ny = ny;
		this.nz = nz;

		corners = GridSupport.getCellCorners(grid, nx, ny, nz);
		cornerBits = computeCornerBits(corners, radius);

		int numInside = Bits.countBits(cornerBits);

		// Cells entirely inside or outside have no surface intersection.
		if (numInside == 0 || numInside == 8) {
			intersectionType = IntersectionType.UNKNOWN;
			return;
		}

		buildEdgesAndPolygon(grid, radius);
		intersectionType = classify(cornerBits, numInside);
	}

	// -----------------------------------------------------------------------
	// Private helpers
	// -----------------------------------------------------------------------

	/**
	 * Computes the corner bitmask: bit {@code k} is set if corner {@code k} is
	 * strictly inside the sphere.
	 */
	private static int computeCornerBits(double[][] corners, double radius) {
		int bits = 0;
		double r2 = radius * radius;
		for (int k = 0; k < 8; k++) {
			double x = corners[k][0], y = corners[k][1], z = corners[k][2];
			if (x * x + y * y + z * z < r2) {
				bits |= GridSupport.CORNERBITS[k];
			}
		}
		return bits;
	}

	/**
	 * Finds the intersecting edges, constructs {@link Edge} objects, reorders them
	 * into a closed loop, and extracts the intersection polygon.
	 */
	private void buildEdgesAndPolygon(CartesianGrid grid, double radius) {
		int[] edgeIndices = GridSupport.findIntersectingEdges(cornerBits);
		if (edgeIndices.length == 0) {
			return;
		}

		Edge[] rawEdges = new Edge[edgeIndices.length];
		for (int i = 0; i < edgeIndices.length; i++) {
			int[] edgeCorners = GridSupport.getCornersOfEdges(edgeIndices[i]);
			rawEdges[i] = new Edge(grid, edgeCorners[0], edgeCorners[1], nx, ny, nz, radius);
		}

		orderedEdges = EdgeOrdering.reorderEdges(rawEdges);

		intersectionPolygon = new Point3D.Double[orderedEdges.length];
		for (int i = 0; i < orderedEdges.length; i++) {
			intersectionPolygon[i] = orderedEdges[i].getIntersection();
		}
	}

	/**
	 * Assigns an {@link IntersectionType} from the corner bitmask.
	 *
	 * <p>
	 * The mapping follows the standard marching-cubes topology for a sphere
	 * intersecting a cube:
	 * </p>
	 * <ul>
	 * <li>1 or 7 corners inside → {@code CORNERIN} / {@code CORNEROUT}</li>
	 * <li>2 or 6 corners inside, adjacent on one face → {@code DOUBLECORNERIN} /
	 * {@code DOUBLECORNEROUT}</li>
	 * <li>2 or 6 corners inside, sharing one edge → face cut variants</li>
	 * <li>3 or 5 corners inside → {@code SKEWCUT} or {@code CORNERPULL} /
	 * {@code CORNERPUSH}</li>
	 * <li>4 corners inside, forming a face → {@code FACECUT}</li>
	 * <li>Degenerate tangent cases → {@code KISS}</li>
	 * </ul>
	 *
	 * <p>
	 * By convention the type is always described from the minority side (≤ 4 inside
	 * corners); the complement is handled symmetrically.
	 * </p>
	 */
	private static IntersectionType classify(int cornerBits, int numInside) {
		// Normalise to the minority side for lookup symmetry.
		int bits = (numInside <= 4) ? cornerBits : (~cornerBits & 0xFF);
		int n = (numInside <= 4) ? numInside : (8 - numInside);

		return switch (n) {
		case 1 -> (numInside == 1) ? IntersectionType.CORNERIN : IntersectionType.CORNEROUT;
		case 2 -> classifyTwo(bits, numInside);
		case 3 -> classifyThree(bits, numInside);
		case 4 -> classifyFour(bits);
		default -> IntersectionType.UNKNOWN;
		};
	}

	/** Two corners inside: either a double-corner clip or a face-edge cut. */
	private static IntersectionType classifyTwo(int bits, int numInside) {
		// Two corners share a face if they appear together in any face definition.
		for (int face = 0; face < 6; face++) {
			int[] fc = GridSupport.getFaceCornerIndices(face);
			int count = 0;
			for (int c : fc) {
				if ((bits & GridSupport.CORNERBITS[c]) != 0)
					count++;
			}
			if (count == 2) {
				return (numInside == 2) ? IntersectionType.DOUBLECORNERIN : IntersectionType.DOUBLECORNEROUT;
			}
		}
		// Both inside corners share only an edge, not a face — skewed slice.
		return IntersectionType.SKEWCUT;
	}

	/** Three corners inside: corner-pull, corner-push, or skew. */
	private static IntersectionType classifyThree(int bits, int numInside) {
		// If all three lie on one face it's a pull/push; otherwise skew.
		for (int face = 0; face < 6; face++) {
			int[] fc = GridSupport.getFaceCornerIndices(face);
			int count = 0;
			for (int c : fc) {
				if ((bits & GridSupport.CORNERBITS[c]) != 0)
					count++;
			}
			if (count == 3) {
				return (numInside == 3) ? IntersectionType.CORNERPULL : IntersectionType.CORNERPUSH;
			}
		}
		return IntersectionType.SKEWCUT;
	}

	/** Four corners inside: clean face cut if they form a face, otherwise skew. */
	private static IntersectionType classifyFour(int bits) {
		for (int face = 0; face < 6; face++) {
			int[] fc = GridSupport.getFaceCornerIndices(face);
			int count = 0;
			for (int c : fc) {
				if ((bits & GridSupport.CORNERBITS[c]) != 0)
					count++;
			}
			if (count == 4)
				return IntersectionType.FACECUT;
		}
		return IntersectionType.SKEWCUT;
	}

	/**
	 * Sets the intersection type for this cell.
	 * 
	 * @param type the IntersectionType to assign to this cell
	 */
	public void setIntersectionType(IntersectionType type) {
		this.intersectionType = type;
	}

	// -----------------------------------------------------------------------
	// Accessors
	// -----------------------------------------------------------------------

	/** @return the inside-corner bitmask */
	public int getCornerBits() {
		return cornerBits;
	}

	/** @return the eight cell corners in canonical order */
	public double[][] getCorners() {
		return corners;
	}

	/** @return the intersection type */
	public IntersectionType getIntersectionType() {
		return intersectionType;
	}

	/** @return true if this cell contains a spherical-grid pole */
	public boolean isPolar() {
		return polar;
	}

	/** @param polar true if this cell straddles a pole */
	public void setPolar(boolean polar) {
		this.polar = polar;
	}

	/**
	 * @return the intersecting edges ordered into a closed loop, or {@code null} if
	 *         the cell has no surface intersection
	 */
	public Edge[] getOrderedEdges() {
		return orderedEdges;
	}

	/**
	 * @return the sphere-surface polygon vertices in cyclic order, or {@code null}
	 *         if the cell has no surface intersection
	 */
	public Point3D.Double[] getIntersectionPolygon() {
		return intersectionPolygon;
	}

	/**
	 * @return true if the sphere surface intersects this cell
	 */
	public boolean hasIntersection() {
		return orderedEdges != null && orderedEdges.length > 0;
	}
	
	/**
	 * @return true if this cell is a degenerate tangent case (a "kiss" with the
	 *         sphere)
	 */
	public boolean isKiss() {
		return intersectionType == IntersectionType.KISS;
	}
	
	/**
	 * Gets the KissGeometry for this cell, which contains information about the
	 * @return the KissGeometry for this cell, or null if this cell is not a kiss case
	 */
	public KissGeometry getKissGeometry() {
		return kissGeometry;
	}
	
	/**
	 * Sets the KissGeometry for this cell.
	 * @param kissGeometry the KissGeometry to set for this cell
	 */
	public void setKissGeometry(KissGeometry kissGeometry) {
		if (this.intersectionType != IntersectionType.KISS) {
			throw new IllegalStateException("Can only set KissGeometry for cells classified as KISS");
		}
		this.kissGeometry = kissGeometry;
	}

	@Override
	public String toString() {
		return String.format("Cell[%d,%d,%d] type=%s polar=%b corners=0x%02X", nx, ny, nz, intersectionType, polar,
				cornerBits);
	}
}