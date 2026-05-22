package edu.cnu.mdi.chimera.area;

import java.util.List;

import edu.cnu.mdi.chimera.util.Point3D;

/**
 * Computes the area of a spherical polygon from sampled boundary points.
 *
 * <h2>Algorithm</h2>
 * <p>Fan triangulation from vertex 0 using the Oosterom-Strackee formula
 * for each spherical triangle:</p>
 * <pre>
 *   area = 2 · atan2( det(a, b, c),  1 + a·b + b·c + c·a )
 * </pre>
 * <p>where {@code det(a,b,c) = a · (b × c)} and {@code a, b, c} are unit
 * vectors. The signed area is accumulated across the fan and then wrapped
 * into {@code (-2π, 2π]} before taking the absolute value.</p>
 *
 * <h2>Polar patches</h2>
 * <p>No special pole-enclosure flag is required. When the fan from vertex 0
 * covers the complement of the enclosed region (the large cap), the signed
 * area exceeds 2π and the wrapping step reduces it to the correct small value.
 * When the fan covers the enclosed region directly, the signed area is already
 * small. In both cases {@code abs(wrapped)} gives the correct enclosed area.</p>
 *
 * <h2>Input</h2>
 * <p>Vertices are {@link Point3D.Double} Cartesian points on the sphere
 * surface. They are normalized to unit vectors internally, so the returned
 * area is on the unit sphere (steradians). Divide by {@code 4π} to obtain
 * the fraction of the full sphere.</p>
 */
public class SphericalPolygonArea {

    private static final double TOL = 1.0e-14;
    private static final int    MIN_VERTICES = 3;

    private SphericalPolygonArea() {}

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Returns the area of the spherical polygon as a fraction of the full
     * sphere area {@code 4πR²}, in {@code [0, 1]}.
     *
     * @param vertices ordered boundary vertices on the sphere; at least 3 required
     * @return area fraction in {@code [0, 1]}
     */
    public static double computeAreaFraction(List<Point3D.Double> vertices) {
        return signedAreaUnitSphere(vertices) / (4.0 * Math.PI);
    }

    /**
     * Returns the area of the polygon in steradians (on the unit sphere),
     * in {@code [0, 4π]}.
     *
     * @param vertices ordered boundary vertices; at least 3 required
     * @return unsigned area in steradians
     */
    public static double computeSteradians(List<Point3D.Double> vertices) {
        return signedAreaUnitSphere(vertices);
    }

    // -----------------------------------------------------------------------
    // Core algorithm
    // -----------------------------------------------------------------------

    /**
     * Computes the unsigned area of the spherical polygon on the unit sphere.
     *
     * <p>The signed fan area is wrapped into {@code (-2π, 2π]} and then
     * {@code abs()} is taken. This correctly handles both polar and non-polar
     * patches regardless of boundary winding direction.</p>
     */
    private static double signedAreaUnitSphere(List<Point3D.Double> vertices) {
        if (vertices == null || vertices.size() < MIN_VERTICES) return 0.0;

        // Normalize to unit vectors, skipping degenerate points.
        double[][] uv = new double[vertices.size()][3];
        int m = 0;
        for (Point3D.Double p : vertices) {
            double len = Math.sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            if (!Double.isFinite(len) || len < TOL) continue;
            uv[m][0] = p.x/len;  uv[m][1] = p.y/len;  uv[m][2] = p.z/len;
            m++;
        }
        if (m < MIN_VERTICES) return 0.0;

        // Fan triangulation from uv[0].
        double[] a = uv[0];
        double area = 0.0;
        for (int i = 1; i < m - 1; i++) {
            double[] b = uv[i], c = uv[i + 1];
            // det = a · (b × c)
            double det  = a[0]*(b[1]*c[2]-b[2]*c[1])
                        + a[1]*(b[2]*c[0]-b[0]*c[2])
                        + a[2]*(b[0]*c[1]-b[1]*c[0]);
            double denom = 1.0
                    + (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
                    + (b[0]*c[0] + b[1]*c[1] + b[2]*c[2])
                    + (a[0]*c[0] + a[1]*c[1] + a[2]*c[2]);
            double tri = 2.0 * Math.atan2(det, denom);
            if (Double.isFinite(tri)) area += tri;
        }

        // Wrap into (-2π, 2π] to handle polar-patch fan complement.
        double full = 4.0 * Math.PI;
        while (area >  2.0 * Math.PI) area -= full;
        while (area < -2.0 * Math.PI) area += full;

        return Math.abs(area);
    }
}
