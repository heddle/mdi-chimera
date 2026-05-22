package edu.cnu.mdi.chimera.grid;

import java.util.ArrayList;
import java.util.List;

import edu.cnu.mdi.util.Bits;

public class GridSupport {

	//used for bitwise operations
	public static final int CORNER0BIT = 01;
	public static final int CORNER1BIT = 02;
	public static final int CORNER2BIT = 04;
	public static final int CORNER3BIT = 010;
	public static final int CORNER4BIT = 020;
	public static final int CORNER5BIT = 040;
	public static final int CORNER6BIT = 0100;
	public static final int CORNER7BIT = 0200;

	public static final int[] CORNERBITS = { CORNER0BIT, CORNER1BIT, CORNER2BIT, CORNER3BIT, CORNER4BIT, CORNER5BIT,
			CORNER6BIT, CORNER7BIT };

	/**
	 * Get the indices of the corners of the given face.
	 *
	 * @param face [0, 5] for the faces of a cube. This assumes the canonical
	 *             corner numbering scheme:
	 *             <ol>
	 *             <li value="0">Corner 0 at (x0, y0, z0)</li>
	 *             <li value="1">Corner 1 at (x1, y0, z0)</li>
	 *             <li value="2">Corner 2 at (x0, y1, z0)</li>
	 *             <li value="3">Corner 3 at (x1, y1, z0)</li>
	 *             <li value="4">Corner 4 at (x0, y0, z1)</li>
	 *             <li value="5">Corner 5 at (x1, y0, z1)</li>
	 *             <li value="6">Corner 6 at (x0, y1, z1)</li>
	 *             <li value="7">Corner 7 at (x1, y1, z1)</li>
	 *             </ol>
	 *
	 *             The faces are numbered as follows:
	 *             <ol>
	 *             <li value="0">Face 0 (xy plane at z=z0)</li>
	 *             <li value="1">Face 1 (xy plane at z=z1)</li>
	 *             <li value="2">Face 2 (xz plane at y=y0)</li>
	 *             <li value="3">Face 3 (xz plane at y=y1)</li>
	 *             <li value="4">Face 4 (yz plane at x=x0)</li>
	 *             <li value="5">Face 5 (yz plane at x=x1)</li>
	 * @return
	 */
    public static int[] getFaceCornerIndices(int face) {
        return switch (face) {
            case 0 -> new int[]{0, 1, 3, 2};
            case 1 -> new int[]{4, 5, 7, 6};
            case 2 -> new int[]{0, 1, 5, 4};
            case 3 -> new int[]{2, 3, 7, 6};
            case 4 -> new int[]{0, 2, 6, 4};
            case 5 -> new int[]{1, 3, 7, 5};
            default -> throw new IllegalArgumentException("Invalid face index.");
        };
    }

    /**
     * Returns the indices of the three faces on which the given corner lies.
     * The canonical numbering for the corners is assumed:
     * <pre>
     * Corner 0: (x0, y0, z0)
     * Corner 1: (x1, y0, z0)
     * Corner 2: (x0, y1, z0)
     * Corner 3: (x1, y1, z0)
     * Corner 4: (x0, y0, z1)
     * Corner 5: (x1, y0, z1)
     * Corner 6: (x0, y1, z1)
     * Corner 7: (x1, y1, z1)
     * </pre>
     * The faces are numbered as:
     * <pre>
     * Face 0: {0, 1, 3, 2} (xy plane at z=z0)
     * Face 1: {4, 5, 7, 6} (xy plane at z=z1)
     * Face 2: {0, 1, 5, 4} (xz plane at y=y0)
     * Face 3: {2, 3, 7, 6} (xz plane at y=y1)
     * Face 4: {0, 2, 6, 4} (yz plane at x=x0)
     * Face 5: {1, 3, 7, 5} (yz plane at x=x1)
     * </pre>
     *
     * @param cornerIndex the index of the corner (0 to 7)
     * @return an int array of length 3 containing the face indices for the corner.
     * @throws IllegalArgumentException if the cornerIndex is not in [0, 7].
     */
    public static int[] getCornerFaces(int cornerIndex) {
        if (cornerIndex < 0 || cornerIndex > 7) {
            throw new IllegalArgumentException("Corner index must be between 0 and 7.");
        }

        // Determine the face for the z-coordinate.
        // If the third bit (value 4) is set, then z = z1, and the face is 1.
        // Otherwise, z = z0 and the face is 0.
        int faceZ = ((cornerIndex & 4) == 0) ? 0 : 1;

        // Determine the face for the y-coordinate.
        // If the second bit (value 2) is set, then y = y1, and the face is 3.
        // Otherwise, y = y0 and the face is 2.
        int faceY = ((cornerIndex & 2) == 0) ? 2 : 3;

        // Determine the face for the x-coordinate.
        // If the least-significant bit (value 1) is set, then x = x1, and the face is 5.
        // Otherwise, x = x0 and the face is 4.
        int faceX = ((cornerIndex & 1) == 0) ? 4 : 5;

        return new int[] { faceZ, faceY, faceX };
    }

    /**
     * Returns the two face indices (as an int[2]) that the given edge is on.
     * The edge is specified by its canonical index (0-11) as defined by
     * getCornersOfEdges.
     *
     * @param edgeIndex the canonical edge index (0 to 11)
     * @return an int array of length 2 containing the indices of the two faces that the edge lies on.
     * @throws IllegalArgumentException if the edge index is invalid or if the edge does not lie on exactly two faces.
     */
    public static int[] getEdgeFaces(int edgeIndex) {
        // Get the two corner indices that define the edge.
        int[] edgeCorners = getCornersOfEdges(edgeIndex);
        int cornerA = edgeCorners[0];
        int cornerB = edgeCorners[1];

        // Get the faces associated with each corner.
        int[] facesA = getCornerFaces(cornerA);  // returns three face indices for cornerA
        int[] facesB = getCornerFaces(cornerB);  // returns three face indices for cornerB

        // Find the common faces between the two corners.
        List<Integer> commonFaces = new ArrayList<>(2);
        for (int face : facesA) {
            for (int faceB : facesB) {
                if (face == faceB) {
                    commonFaces.add(face);
                }
            }
        }

        // There should be exactly two faces in common.
        if (commonFaces.size() != 2) {
            throw new IllegalArgumentException("Edge " + edgeIndex + " does not lie on exactly two faces.");
        }

        // Convert the result to an int[] and return.
        int[] result = new int[2];
        result[0] = commonFaces.get(0);
        result[1] = commonFaces.get(1);
        return result;
    }


    /**
     * Get the canonical edge index for the given pair of corners.
     * @param edge the pair of corners that define the edge. The corners are
     *            numbered as follows:
     *            <ol>
     *            <li value="0">Corner 0 at (x0, y0, z0)</li>
     *            <li value="1">Corner 1 at (x1, y0, z0)</li>
     *            <li value="2">Corner 2 at (x0, y1, z0)</li>
     *            <li value="3">Corner 3 at (x1, y1, z0)</li>
     *            <li value="4">Corner 4 at (x0, y0, z1)</li>
     *            <li value="5">Corner 5 at (x1, y0, z1)</li>
     *            <li value="6">Corner 6 at (x0, y1, z1)</li>
     *            <li value="7">Corner 7 at (x1, y1, z1)</li>
     *            </ol>
     *
     * @return the canonical edge index for the given pair of corners.
     */
    public static int[] getCornersOfEdges(int edge) {
		return switch (edge) {
		case 0 -> new int[] { 0, 1 };
		case 1 -> new int[] { 0, 2 };
		case 2 -> new int[] { 0, 4 };
		case 3 -> new int[] { 1, 3 };
		case 4 -> new int[] { 1, 5 };
		case 5 -> new int[] { 2, 3 };
		case 6 -> new int[] { 2, 6 };
		case 7 -> new int[] { 3, 7 };
		case 8 -> new int[] { 4, 5 };
		case 9 -> new int[] { 4, 6 };
		case 10 -> new int[] { 5, 7 };
		case 11 -> new int[] { 6, 7 };


		default -> throw new IllegalArgumentException("Invalid edge index.");
		};
    }


    /**
     * Determines which edges of the rectangular prism (cell) intersect the surface
     * of the sphere. An edge is intersecting if one of its corners is inside the sphere
     * (its bit is set in cornerBits) and the other is outside.
     *
     * @param cornerBits an int with the inside corners set (using bit constants)
     * @return an array containing the indices (0-11) of the intersecting edges.
     */
    public static int[] findIntersectingEdges(int cornerBits) {

    	int numInside = Bits.countBits(cornerBits);
    	if (numInside == 0) {
			return new int[0];
		}
    	if (numInside > 7) {
    		System.err.println("findIntersectingEdges: too many inside corners");
    		System.exit(-1);
    	}


        List<Integer> intersectingEdges = new ArrayList<>();

        // Loop over all 12 edges.
        for (int edge = 0; edge < 12; edge++) {
            int[] edgeCorners = getCornersOfEdges(edge);

            // Get the bit constants for the two corners of the edge.
            int bitA = CORNERBITS[edgeCorners[0]];
            int bitB = CORNERBITS[edgeCorners[1]];

            // Use the Bits utility class to check if the bit is set.
            boolean inA = Bits.check(cornerBits, bitA);
            boolean inB = Bits.check(cornerBits, bitB);

            // If one is inside and the other is outside, then the edge must intersect.
            if (inA != inB) {
                intersectingEdges.add(edge);
            }
        }

        // Convert the List<Integer> to int[].
        int[] result = new int[intersectingEdges.size()];
        for (int i = 0; i < intersectingEdges.size(); i++) {
            result[i] = intersectingEdges.get(i);
        }
        return result;
    }

	/**
	 * Get the canonical edge index for the given pair of corners.
	 *
	 * @param corner1 the first corner
	 * @param corner2 the second corner
	 * @return the canonical edge index for the given pair of corners.
	 * If the corners don't make a valid edge, -1 is returned.
	 */
    public static int getEdgeIndex(int corner1, int corner2) {
		if (corner1 < 0 || corner1 > 7 || corner2 < 0 || corner2 > 7) {
			return -1;
		}

		if (corner1 == corner2) {
			return -1;
		}

		// swap if necessary
		if (corner1 > corner2) {
			int temp = corner1;
			corner1 = corner2;
			corner2 = temp;
		}

		return switch (corner1) {
		case 0 -> switch (corner2) {
		case 1 -> 0;
		case 2 -> 1;
		case 4 -> 2;
		default -> -1;

		};
		case 1 -> switch (corner2) {
		case 3 -> 3;
		case 5 -> 4;
		default -> -1;
		};

		case 2 -> switch (corner2) {
		case 3 -> 5;
		case 6 -> 6;
		default -> -1;
		};

		case 3 -> switch (corner2) {
		case 7 -> 7;
		default -> -1;
		};

		case 4 -> switch (corner2) {
		case 5 -> 8;
		case 6 -> 9;
		default -> -1;
		};

		case 5 -> switch (corner2) {
		case 7 -> 10;
		default -> -1;
		};

		case 6 -> switch (corner2) {
		case 7 -> 11;
		default -> -1;
		};

		default -> -1;
		};
	}


	/**
	 * Get the indices of the corners of the given face.
	 *
	 * @param corners the 8 cell corners stored as an array of 8 double arrays,
	 *               each with 3 elements (x, y, z).
	 * @param face    [0, 5] for the faces of a cube. This assumes the canonical
	 *                corner numbering scheme:
	 *                <ol>
	 *                <li value="0">Corner 0 at (x0, y0, z0)</li>
	 *                <li value="1">Corner 1 at (x1, y0, z0)</li>
	 *                <li value="2">Corner 2 at (x0, y1, z0)</li>
	 *                <li value="3">Corner 3 at (x1, y1, z0)</li>
	 *                <li value="4">Corner 4 at (x0, y0, z1)</li>
	 *                <li value="5">Corner 5 at (x1, y0, z1)</li>
	 *                <li value="6">Corner 6 at (x0, y1, z1)</li>
	 *                <li value="7">Corner 7 at (x1, y1, z1)</li>
	 *                </ol>
	 *
	 *                The faces are numbered as follows:
	 *                <ol>
	 *                <li value="0">Face 0 (xy plane at z=z0)</li>
	 *                <li value="1">Face 1 (xy plane at z=z1)</li>
	 *                <li value="2">Face 2 (xz plane at y=y0)</li>
	 *                <li value="3">Face 3 (xz plane at y=y1)</li>
	 *                <li value="4">Face 4 (yz plane at x=x0)</li>
	 *                <li value="5">Face 5 (yz plane at x=x1)</li>
	 * @return the 4 corners of the face from the 8 cell corners. Will be an array
	 *        of 4 double arrays, each with 3 elements (x, y, z).
	 */

	public static double[][] getFaceCorners(double[][] corners, int face) {
		int[] cornerIndices = getFaceCornerIndices(face);
		double[][] faceCorners = new double[4][3];
		for (int i = 0; i < 4; i++) {
			faceCorners[i] = corners[cornerIndices[i]];
		}
		return faceCorners;
	}


    /**
     * Compute the centroid of the face.
 	 * @param corners the 8 cell corners stored as an array of 8 double arrays,
	 *              each with 3 elements (x, y, z).
	 * @param face [0, 5] for the faces of a cube. This assumes the canonical numbering.
    * @return the centroid
     */
    public static double[] computeCentroid(double[][] corners, int face) {
    	double[][] faceCorners = getFaceCorners(corners, face);
        double[] centroid = {0, 0, 0};
        for (int i = 0; i < 4; i++) {
            centroid[0] += faceCorners[i][0];
            centroid[1] += faceCorners[i][1];
            centroid[2] += faceCorners[i][2];
        }
        centroid[0] /= 4.0;
        centroid[1] /= 4.0;
        centroid[2] /= 4.0;
        return centroid;
    }


	/**
	 * Get the indices of the corners of the given face.
	 * @param corners the 8 cell corners stored as an array of 8 double arrays,
	 *              each with 3 elements (x, y, z).
	 * @param face [0, 5] for the faces of a cube. This assumes the canonical numbering.
	 * @return the average distance squared of the face corners from the origin,
	 * which is the center of the sphere.
	 */
	public static double faceAverageDistanceSquare(double[][] corners, int face) {
		double[][] faceCorners = getFaceCorners(corners, face);
		double rsqsum = 0;

		for (double[] corner : faceCorners) {
			double x = corner[0], y = corner[1], z = corner[2];
			double distsq = x * x + y * y + z * z;
			rsqsum += distsq;
		}
		return rsqsum / 4;
	}

	/**
	 * Get the face that is closest to the origin based on its centroid
	 * @param corners the corners of the cell
	 * @return the face index [0, 5] that is closest
	 */
	public static int getClosestFaceToOrigin(double[][] corners) {
    	int closestFace = -1;
    	double minDist = Double.MAX_VALUE;
		for (int face = 0; face < 6; face++) {

			double[] centroid = GridSupport.computeCentroid(corners, face);
			double distsq = centroid[0]*centroid[0] + centroid[1]*centroid[1] + centroid[2]*centroid[2];


//			double distsq = GridSupport.faceAverageDistanceSquare(corners, face);
			if (distsq < minDist) {
				minDist = distsq;
				closestFace = face;
			}
		}
        return closestFace;
	}




    /**
     * Returns the Cartesian coordinates of the corners of the cell at the given indices.
     * @param ix the index of the smaller x coordinate
     * @param iy the index of the smaller y coordinate
     * @param iz the index of the smaller z coordinate
     * @return the Cartesian coordinates of the eight cell corners. The corners
     * are ordered in the canonical numbering as follows:
     * <ol>
     * <li value="0">Corner 0 at (x0, y0, z0)</li>
     * <li value="1">Corner 1 at (x1, y0, z0)</li>
     * <li value="2">Corner 2 at (x0, y1, z0)</li>
     * <li value="3">Corner 3 at (x1, y1, z0)</li>
     * <li value="4">Corner 4 at (x0, y0, z1)</li>
     * <li value="5">Corner 5 at (x1, y0, z1)</li>
     * <li value="6">Corner 6 at (x0, y1, z1)</li>
     * <li value="7">Corner 7 at (x1, y1, z1)</li>
     * </ol>
     * The corners are returned as an array of 8 double arrays, each with 3 elements
     * (x, y, z).
     * @return the Cartesian coordinates of the corners of the cell at the given
     * indices. Will be an array of 8 double arrays, each with 3 elements (x, y, z).
     */
    public static double[][] getCellCorners(CartesianGrid cartesianGrid, int ix, int iy, int iz) {
        double[][] corners = new double[8][3];

        Grid1D xGrid = cartesianGrid.getXGrid();
        Grid1D yGrid = cartesianGrid.getYGrid();
        Grid1D zGrid = cartesianGrid.getZGrid();

        double x0 = xGrid.valueAt(ix);
        double y0 = yGrid.valueAt(iy);
        double z0 = zGrid.valueAt(iz);

        double x1 = xGrid.valueAt(ix + 1);
        double y1 = yGrid.valueAt(iy + 1);
        double z1 = zGrid.valueAt(iz + 1);

        //note the canonical order of the corners
        corners[0] = new double[]{x0, y0, z0};
        corners[1] = new double[]{x1, y0, z0};
        corners[2] = new double[]{x0, y1, z0};
        corners[3] = new double[]{x1, y1, z0};
        corners[4] = new double[]{x0, y0, z1};
        corners[5] = new double[]{x1, y0, z1};
        corners[6] = new double[]{x0, y1, z1};
        corners[7] = new double[]{x1, y1, z1};

        return corners;
    }

    /**
     * Get the Cartesian coordinates of the corner of the cell at the given indices.
     * @param cartesianGrid the Cartesian grid
     * @param ix the index of the smaller x coordinate
     * @param iy the index of the smaller y coordinate
     * @param iz the index of the smaller z coordinate
     * @param corner the corner index [0, 7] for the corners of a cube. This assumes the canonical numbering
     * @return the Cartesian coordinates of the corner of the cell at the given indices.
     */
	public static double[] getCellCorner(CartesianGrid cartesianGrid, int ix, int iy, int iz, int corner) {
		double[][] corners = getCellCorners(cartesianGrid, ix, iy, iz);
		return corners[corner];
	}

}
