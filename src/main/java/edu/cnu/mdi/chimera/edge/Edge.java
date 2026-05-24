package edu.cnu.mdi.chimera.edge;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.GridSupport;
import edu.cnu.mdi.chimera.util.Point3D;

public class Edge {
	
	private static final double TOL = 1.0e-14;

	private Point3D.Double startPoint;
	private Point3D.Double endPoint;
	private Point3D.Double intersection;

	//each corner is on three faces
	int[][] cornerfaces;

	//each edge is on two faces
	private int[] edgeFaces;

	/**
	 * Constructor for an Edge.
	 * @param grid The Cartesian grid.
	 * @param corner0 The canonical index of the first corner.
	 * @param corner1 The canonical index of the second corner.
	 * @param nx The x index of the cell in the grid.
	 * @param ny The y index of the cell in the grid.
	 * @param nz The z index of the cell in the grid.
	 * @param radius The radius of the sphere.
     */
	public Edge(CartesianGrid grid, int corner0, int corner1, int nx, int ny, int nz, double radius) {

		startPoint = new Point3D.Double(GridSupport.getCellCorner(grid, nx, ny, nz, corner0));
		endPoint = new Point3D.Double(GridSupport.getCellCorner(grid, nx, ny, nz, corner1));

		cornerfaces = new int[2][];
		cornerfaces[0] = GridSupport.getCornerFaces(corner0);
		cornerfaces[1] = GridSupport.getCornerFaces(corner1);

		//each edge is on two faces
		int gridIndex = GridSupport.getEdgeIndex(corner0, corner1);
		if (gridIndex < 0) {
			throw new IllegalArgumentException("Invalid edge indices: " + corner0 + ", " + corner1);
		}
		edgeFaces = GridSupport.getEdgeFaces(gridIndex);

		//find the intersection point
        intersection = findSphereIntersection(startPoint, endPoint, radius);
	}

	/**
	 * Get the faces that the edge is on.
	 *
	 * @return The faces.
	 */
	public int[] getEdgeFaces() {
		return edgeFaces;
	}

	/**
	 * Get the starting point of the edge.
	 *
	 * @return The starting point.
	 */
	public Point3D.Double getStartPoint() {
		return startPoint;
	}

	/**
	 * Get the ending point of the edge.
	 * @return The ending point.
	 */
	public Point3D.Double getEndPoint() {
		return endPoint;
	}


	/**
	 * Get a common face with this and another edge
	 * @param other The other edge.
	 * @return A common face or -1 if no common face
	 */
	public int getCommonFace(Edge other) {
		for (int element : edgeFaces) {
			for (int element2 : other.edgeFaces) {
				if (element == element2) {
					return element;
				}
			}
		}
		return -1;
	}

	/**
	 * Get the intersection point of the edge with the sphere.
	 *
	 * @return The intersection point.
	 */
	public Point3D.Double getIntersection() {
		return intersection;
	}

	/**
     * Finds the intersection point on the segment p0->p1 with a sphere centered at the origin.
     * One point must be inside the sphere and one must be outside; otherwise, an exception is thrown.
     *
     * @param p0     One endpoint of the segment.
     * @param p1     The other endpoint of the segment.
     * @param radius The radius of the sphere centered at the origin.
     * @return The intersection point on the segment.
     * @throws IllegalArgumentException if both points are inside or both are outside the sphere,
     *                                  or if the intersection does not lie on the segment.
     */
    public static Point3D.Double findSphereIntersection(Point3D.Double p0, Point3D.Double p1, double radius) {
        // Compute the distances of the points from the origin.
        double d0 = p0.length();
        double d1 = p1.length();

        // Check that one point is inside (distance < radius) and one is outside (distance > radius).
        if ((d0 < radius && d1 < radius) || (d0 > radius && d1 > radius)) {
            throw new IllegalArgumentException(
                "One point must be inside and one point must be outside the sphere."
            );
        }

        // Compute the direction vector from p0 to p1.
        double dx = p1.x - p0.x;
        double dy = p1.y - p0.y;
        double dz = p1.z - p0.z;

        // The parametric equation of the line is:
        //   p(t) = p0 + t * (p1 - p0),   for t in [0,1].
        // We need to find t such that ||p(t)||^2 = radius^2.
        //
        // Write p(t) = (p0.x + t*dx, p0.y + t*dy, p0.z + t*dz) and square the norm:
        //   (p0.x + t*dx)^2 + (p0.y + t*dy)^2 + (p0.z + t*dz)^2 = radius^2.
        // This expands to:
        //   A*t^2 + B*t + C = 0,
        // where:
        double A = dx * dx + dy * dy + dz * dz;
        double B = 2 * (p0.x * dx + p0.y * dy + p0.z * dz);
        double C = p0.x * p0.x + p0.y * p0.y + p0.z * p0.z - radius * radius;

        // Compute the discriminant.
        double discriminant = B * B - 4 * A * C;
        if (discriminant < 0) {
            throw new IllegalArgumentException("No intersection exists between the line and the sphere.");
        }

        double sqrtD = Math.sqrt(discriminant);
        double t;

        // Choose the correct root based on which point is inside.
        // If p0 is inside (d0 < radius), we must go from p0 outwards to exit the sphere.
        // Thus, we choose the larger t (i.e., t = (-B + sqrtD) / (2*A)).
        // Otherwise, if p0 is outside (and p1 inside), we choose the smaller t.
        if (d0 < radius) {
            t = (-B + sqrtD) / (2 * A);
        } else {
            t = (-B - sqrtD) / (2 * A);
        }
        
        if (t < -TOL || t > 1 + TOL) {
			throw new IllegalArgumentException("The intersection does not lie on the segment. t = " + t);
		}
        t = Math.max(0, Math.min(1, t)); // Clamp t to [0, 1] to handle numerical issues.

        // Compute the intersection point.
        double ix = p0.x + t * dx;
        double iy = p0.y + t * dy;
        double iz = p0.z + t * dz;

        return new Point3D.Double(ix, iy, iz);
    }

}
