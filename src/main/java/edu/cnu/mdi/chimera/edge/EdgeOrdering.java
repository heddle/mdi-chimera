package edu.cnu.mdi.chimera.edge;

import java.util.ArrayList;
import java.util.List;

/**
 * Utility class that reorders an array of Edge objects so that
 * each consecutive pair of edges shares at least one face, and
 * the last edge shares a face with the first edge.
 */
public class EdgeOrdering {

    /**
     * Reorders the given array of edges so that each adjacent pair shares a common face,
     * and the last edge shares a common face with the first. If no such ordering exists,
     * an IllegalArgumentException is thrown.
     *
     * @param edges an array of Edge objects that intersect a sphere.
     * @return a new array of Edge objects ordered to form a continuous loop.
     * @throws IllegalArgumentException if the edges cannot be ordered to form a loop.
     */
    public static Edge[] reorderEdges(Edge[] edges) {
        if (edges == null || edges.length == 0) {
            return edges;
        }

        // Try each edge as a possible starting edge.
        for (int start = 0; start < edges.length; start++) {
            boolean[] used = new boolean[edges.length];
            List<Edge> ordering = new ArrayList<>();

            // start with edges[start]
            used[start] = true;
            ordering.add(edges[start]);

            if (search(ordering, used, edges)) {
                // Found a valid ordering.
                return ordering.toArray(new Edge[ordering.size()]);
            }
        }

        // No ordering was found.
        throw new IllegalArgumentException("Cannot order edges to form a closed loop with common faces.");
    }

    /**
     * Recursively searches for an ordering of edges (stored in the 'ordering' list) such that
     * each adjacent pair shares at least one face, and (when complete) the last edge shares a face
     * with the first.
     *
     * @param ordering the current partial ordering.
     * @param used a boolean array marking which edges have already been used.
     * @param edges the complete array of edges.
     * @return true if a complete valid ordering has been found; false otherwise.
     */
    private static boolean search(List<Edge> ordering, boolean[] used, Edge[] edges) {
        // If all edges have been used, check that the last edge connects back to the first.
        if (ordering.size() == edges.length) {
            if (hasCommonFace(ordering.get(ordering.size() - 1), ordering.get(0))) {
                return true;
            }
            return false;
        }

        // Get the last edge in the current ordering.
        Edge lastEdge = ordering.get(ordering.size() - 1);

        // Try all unused edges that share a common face with the last edge.
        for (int i = 0; i < edges.length; i++) {
            if (!used[i] && hasCommonFace(lastEdge, edges[i])) {
                used[i] = true;
                ordering.add(edges[i]);

                if (search(ordering, used, edges)) {
                    return true;
                }

                // Backtrack.
                ordering.remove(ordering.size() - 1);
                used[i] = false;
            }
        }
        return false;
    }

    /**
     * Returns true if the two edges share at least one common face.
     * Each edge is on exactly two faces (accessible via getEdgeFaces()).
     *
     * @param e1 the first edge.
     * @param e2 the second edge.
     * @return true if the two edges share a face; false otherwise.
     */
    private static boolean hasCommonFace(Edge e1, Edge e2) {
        int[] faces1 = e1.getEdgeFaces();
        int[] faces2 = e2.getEdgeFaces();

        // Check for any common face.
        for (int face1 : faces1) {
            for (int face2 : faces2) {
                if (face1 == face2) {
                    return true;
                }
            }
        }
        return false;
    }
}
