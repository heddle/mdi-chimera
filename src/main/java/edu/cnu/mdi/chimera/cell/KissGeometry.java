package edu.cnu.mdi.chimera.cell;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.util.Point3D;
import edu.cnu.mdi.chimera.util.SphericalVector;

import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Computes the intersection boundary and normalized surface area for a
 * {@link IntersectionType#KISS Kiss} cell.
 *
 * <h2>Geometry</h2>
 * <p>A Kiss cell has all eight corners outside the sphere, but the sphere
 * surface bulges through one or more faces. On each axis-aligned face the
 * intersection is a circular arc: the face plane is at perpendicular distance
 * {@code |d|} from the origin, so the circle of intersection has radius
 * {@code rho = sqrt(R² - d²)}, centred at the projection of the origin onto
 * the face plane. The arc is the portion of that circle that stays within
 * the face's rectangular bounds.</p>
 *
 * <h2>Boundary representation</h2>
 * <p>Each arc segment is represented by a {@link KissArc}, which exposes
 * {@link KissArc#theta(double)} and {@link KissArc#phi(double)} as
 * {@link UnivariateFunction}s over a natural parameter {@code t} (the arc
 * angle in radians, running from {@code tStart} to {@code tEnd}).
 * The full boundary is a {@link List} of {@code KissArc}s ordered so that
 * the end of arc {@code k} connects to the start of arc {@code k+1}.</p>
 *
 * <h2>Area</h2>
 * <p>The normalized area {@code A / (4πR²)} is computed from a dense sample
 * of boundary points using the spherical-excess shoelace formula: for each
 * consecutive pair of unit-sphere boundary points the signed area of the
 * spherical triangle they form with an arbitrary reference point is summed.
 * This is exact in the limit of dense sampling and requires no special
 * treatment of the coordinate singularities at the poles.</p>
 */
public class KissGeometry {

    // -----------------------------------------------------------------------
    // Inner class: one arc segment on one face
    // -----------------------------------------------------------------------

    /**
     * One circular arc of the Kiss boundary, lying on an axis-aligned face of
     * the Cartesian cell.
     *
     * <p>The natural parameter {@code t} is the arc angle in radians,
     * running from {@link #tStart} to {@link #tEnd} (both in
     * {@code [-π, π]}). The 3D Cartesian point at parameter {@code t} is:</p>
     * <pre>
     *   inPlane1 = centerU + rho * cos(t)   (clamped to [u0, u1])
     *   inPlane2 = centerV + rho * sin(t)   (clamped to [v0, v1])
     *   fixed    = fixedVal
     * </pre>
     * <p>with the three values assigned to (x, y, z) according to
     * {@link #faceAxis}.</p>
     */
    public static class KissArc {

        /**
         * Identifies which Cartesian axis is perpendicular to the face.
         * 0 = x-face, 1 = y-face, 2 = z-face.
         */
        public final int faceAxis;

        /** The fixed coordinate value of the face plane. */
        public final double fixedVal;

        /** In-plane circle centre, first in-plane axis. */
        public final double centerU;

        /** In-plane circle centre, second in-plane axis. */
        public final double centerV;

        /** Radius of the intersection circle on the face plane. */
        public final double rho;

        /** Start of the arc parameter range (radians). */
        public final double tStart;

        /** End of the arc parameter range (radians). */
        public final double tEnd;

        /** Sphere radius. */
        private final double R;

        /**
         * Constructs a KissArc.
         *
         * @param faceAxis  axis perpendicular to the face (0=x, 1=y, 2=z)
         * @param fixedVal  the face's fixed coordinate
         * @param centerU   circle centre, first in-plane axis
         * @param centerV   circle centre, second in-plane axis
         * @param rho       in-plane circle radius
         * @param tStart    arc start angle (radians)
         * @param tEnd      arc end angle (radians)
         * @param R         sphere radius
         */
        KissArc(int faceAxis, double fixedVal,
                double centerU, double centerV,
                double rho, double tStart, double tEnd, double R) {
            this.faceAxis = faceAxis;
            this.fixedVal = fixedVal;
            this.centerU  = centerU;
            this.centerV  = centerV;
            this.rho      = rho;
            this.tStart   = tStart;
            this.tEnd     = tEnd;
            this.R        = R;
        }

        /**
         * Returns the Cartesian point on the arc at parameter {@code t}.
         *
         * @param t arc angle in radians, in {@code [tStart, tEnd]}
         * @return the 3D point on the sphere surface
         */
        public Point3D.Double cartesian(double t) {
            double u = centerU + rho * Math.cos(t);
            double v = centerV + rho * Math.sin(t);
            return switch (faceAxis) {
                case 0 -> new Point3D.Double(fixedVal, u, v); // x-face: u=y, v=z
                case 1 -> new Point3D.Double(u, fixedVal, v); // y-face: u=x, v=z
                case 2 -> new Point3D.Double(u, v, fixedVal); // z-face: u=x, v=y
                default -> throw new IllegalStateException("Bad faceAxis: " + faceAxis);
            };
        }

        /**
         * Returns the polar angle θ on the sphere as a {@link UnivariateFunction}
         * of the arc parameter {@code t}.
         *
         * @return θ(t) in radians, domain {@code [tStart, tEnd]}
         */
        public UnivariateFunction theta() {
            return t -> {
                Point3D.Double p = cartesian(t);
                return Math.acos(p.z / R);
            };
        }

        /**
         * Returns the azimuthal angle φ on the sphere as a {@link UnivariateFunction}
         * of the arc parameter {@code t}.
         *
         * @return φ(t) in radians in {@code [-π, π]}, domain {@code [tStart, tEnd]}
         */
        public UnivariateFunction phi() {
            return t -> {
                Point3D.Double p = cartesian(t);
                return Math.atan2(p.y, p.x);
            };
        }

        /**
         * Returns the arc length (in radians of arc angle, i.e. {@code tEnd - tStart}).
         *
         * @return arc parameter span
         */
        public double arcSpan() {
            return tEnd - tStart;
        }

        @Override
        public String toString() {
            return String.format(
                "KissArc[axis=%d fixedVal=%.4f rho=%.4f t=[%.4f,%.4f] span=%.4f rad]",
                faceAxis, fixedVal, rho, tStart, tEnd, arcSpan());
        }
    }

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    private final Cell cell;
    private final double R;
    private final double x0, x1, y0, y1, z0, z1;

    /** The ordered list of arc segments forming the closed boundary. */
    private final List<KissArc> arcs;

    /**
     * Number of sample points per arc used for area integration.
     * More points give higher accuracy; 200 per arc is more than sufficient
     * for typical cell sizes relative to the sphere radius.
     */
    private static final int SAMPLES_PER_ARC = 200;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * Constructs a KissGeometry for the given Kiss cell.
     *
     * @param cell the Kiss cell; must have
     *             {@link IntersectionType#KISS} type
     * @param grid the Cartesian grid that owns the cell
     * @param R    the sphere radius
     * @throws IllegalArgumentException if the cell is not a Kiss cell
     */
    public KissGeometry(Cell cell, CartesianGrid grid, double R) {
        if (cell.getIntersectionType() != IntersectionType.KISS) {
            throw new IllegalArgumentException(
                "KissGeometry requires a KISS cell, got: " + cell.getIntersectionType());
        }

        this.cell = cell;
        this.R    = R;

        Grid1D xGrid = grid.getXGrid();
        Grid1D yGrid = grid.getYGrid();
        Grid1D zGrid = grid.getZGrid();

        x0 = xGrid.valueAt(cell.nx);     x1 = xGrid.valueAt(cell.nx + 1);
        y0 = yGrid.valueAt(cell.ny);     y1 = yGrid.valueAt(cell.ny + 1);
        z0 = zGrid.valueAt(cell.nz);     z1 = zGrid.valueAt(cell.nz + 1);

        arcs = buildArcs();
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Returns the ordered list of {@link KissArc}s forming the closed
     * intersection boundary. Each arc is a circular arc on one face of the
     * cell; consecutive arcs share an endpoint.
     *
     * @return unmodifiable list of arcs
     */
    public List<KissArc> getArcs() {
        return java.util.Collections.unmodifiableList(arcs);
    }

    /**
     * Samples the boundary at uniform arc-parameter spacing and returns the
     * points as {@link SphericalVector}s at radius {@code R}.
     *
     * @param samplesPerArc number of sample points per arc segment
     * @return list of boundary points in order around the boundary
     */
    public List<SphericalVector> sampleBoundary(int samplesPerArc) {
        List<SphericalVector> pts = new ArrayList<>();
        for (KissArc arc : arcs) {
            double dt = arc.arcSpan() / samplesPerArc;
            for (int i = 0; i < samplesPerArc; i++) {
                double t = arc.tStart + i * dt;
                Point3D.Double p = arc.cartesian(t);
                pts.add(new SphericalVector(p));
            }
        }
        return pts;
    }

    /**
     * Computes the normalized surface area of the Kiss region:
     * {@code A / (4πR²)}, i.e. as a fraction of the full sphere area.
     *
     * <p>The boundary is sampled at {@link #SAMPLES_PER_ARC} points per arc
     * and the spherical-excess shoelace formula is applied. For a polygon
     * with vertices {@code p_i} on the unit sphere, the signed area is:</p>
     * <pre>
     *   area = |Σ atan2( p_i · (p_{i+1} × p_{i+2}),
     *                    1 + p_i·p_{i+1} + p_{i+1}·p_{i+2} + p_i·p_{i+2} )|
     * </pre>
     * <p>This is the spherical analog of the 2-D shoelace formula and is
     * exact for a spherical polygon with straight (great-circle) edges.</p>
     *
     * @return normalized area in {@code [0, 1]}
     */
    public double normalizedArea() {
        // Collect unit-sphere sample points.
        List<Point3D.Double> pts = new ArrayList<>();
        for (KissArc arc : arcs) {
            double dt = arc.arcSpan() / SAMPLES_PER_ARC;
            for (int i = 0; i < SAMPLES_PER_ARC; i++) {
                double t = arc.tStart + i * dt;
                Point3D.Double p = arc.cartesian(t);
                // Project onto unit sphere.
                double len = p.length();
                pts.add(new Point3D.Double(p.x / len, p.y / len, p.z / len));
            }
        }

        int n = pts.size();
        if (n < 3) return 0.0;

        // Spherical excess shoelace: sum signed triangle areas.
        // Each term uses three consecutive unit-sphere points.
        double signedArea = 0.0;
        for (int i = 0; i < n; i++) {
            Point3D.Double a = pts.get(i);
            Point3D.Double b = pts.get((i + 1) % n);
            Point3D.Double c = pts.get((i + 2) % n);

            // Cross product b × c.
            Point3D.Double cross = new Point3D.Double();
            Point3D.Double.crossProduct(b, c, cross);

            double numerator   = a.dot(cross);                         // a · (b × c)
            double denominator = 1.0 + a.dot(b) + b.dot(c) + a.dot(c);
            signedArea += Math.atan2(numerator, denominator);
        }

        // The shoelace sum gives twice the signed spherical area in steradians.
        // Divide by 2 for the area, then by 4π for normalization.
        double area = Math.abs(signedArea) / 2.0;
        return area / (4.0 * Math.PI);
    }

    // -----------------------------------------------------------------------
    // Arc construction
    // -----------------------------------------------------------------------

    /**
     * Finds all faces that the sphere penetrates and builds the corresponding
     * {@link KissArc}s, then stitches them into a single ordered closed loop.
     *
     * <p>Face layout (matches {@link edu.cnu.mdi.chimera.grid.GridSupport}):</p>
     * <pre>
     *   face 0: z = z0  (faceAxis=2, fixedVal=z0, u=x, v=y)
     *   face 1: z = z1  (faceAxis=2, fixedVal=z1, u=x, v=y)
     *   face 2: y = y0  (faceAxis=1, fixedVal=y0, u=x, v=z)
     *   face 3: y = y1  (faceAxis=1, fixedVal=y1, u=x, v=z)
     *   face 4: x = x0  (faceAxis=0, fixedVal=x0, u=y, v=z)
     *   face 5: x = x1  (faceAxis=0, fixedVal=x1, u=y, v=z)
     * </pre>
     */
    private List<KissArc> buildArcs() {
        List<KissArc> raw = new ArrayList<>();

        // z-faces (faceAxis=2, u=x, v=y)
        tryAddArc(raw, 2, z0, x0, x1, y0, y1);
        tryAddArc(raw, 2, z1, x0, x1, y0, y1);
        // y-faces (faceAxis=1, u=x, v=z)
        tryAddArc(raw, 1, y0, x0, x1, z0, z1);
        tryAddArc(raw, 1, y1, x0, x1, z0, z1);
        // x-faces (faceAxis=0, u=y, v=z)
        tryAddArc(raw, 0, x0, y0, y1, z0, z1);
        tryAddArc(raw, 0, x1, y0, y1, z0, z1);

        if (raw.isEmpty()) {
            throw new IllegalStateException(
                "Kiss cell has no penetrating faces — cell classification error.");
        }

        return raw.size() == 1 ? raw : stitchArcs(raw);
    }

    /**
     * Attempts to build a {@link KissArc} for one face and, if the sphere
     * penetrates that face, adds it to {@code arcs}.
     *
     * <p>The intersection circle on the face has:</p>
     * <ul>
     *   <li>centre at {@code (clamp(0,u0,u1), clamp(0,v0,v1))} — the
     *       projection of the origin onto the face, clamped to the face
     *       quad</li>
     *   <li>radius {@code rho = sqrt(R² - d²)} where {@code d = fixedVal}</li>
     * </ul>
     * <p>The arc is the subset of the circle that lies within {@code [u0,u1]}
     * × {@code [v0,v1]}. The angular bounds are found by intersecting the
     * circle with each face edge and taking the angular range that lies inside
     * the quad.</p>
     *
     * @param arcs     accumulator list
     * @param faceAxis 0=x-face, 1=y-face, 2=z-face
     * @param fixedVal the face's constant coordinate
     * @param u0       lower bound, first in-plane axis
     * @param u1       upper bound, first in-plane axis
     * @param v0       lower bound, second in-plane axis
     * @param v1       upper bound, second in-plane axis
     */
    private void tryAddArc(List<KissArc> arcs, int faceAxis, double fixedVal,
                            double u0, double u1, double v0, double v1) {
        double d2 = fixedVal * fixedVal;
        double R2 = R * R;
        if (d2 >= R2) return; // face plane doesn't intersect sphere

        double rho = Math.sqrt(R2 - d2);

        // The intersection circle on this face plane is always centred at the
        // projection of the origin onto the plane, which in the two in-plane
        // coordinates is simply (0, 0). The clamped-point trick is only used
        // by the scanner to test proximity; here we need the true circle.
        double cu = 0.0;
        double cv = 0.0;

        // Find the arc angular interval within the face bounds [u0,u1]×[v0,v1].
        // The circle is  u = rho*cos(t),  v = rho*sin(t).
        // We need the connected arc where u ∈ [u0,u1] and v ∈ [v0,v1].
        double[] tRange = arcRangeInBox(cu, cv, rho, u0, u1, v0, v1);
        if (tRange == null) return; // circle entirely outside face

        arcs.add(new KissArc(faceAxis, fixedVal, cu, cv, rho,
                             tRange[0], tRange[1], R));
    }

    /**
     * Finds the angular range {@code [tStart, tEnd]} (in radians) of the
     * circle {@code (cu + rho·cos(t), cv + rho·sin(t))} that lies within
     * the box {@code [u0,u1] × [v0,v1]}.
     *
     * <p>Strategy: the circle is clipped by each of the four box edges in
     * turn. Each edge gives at most two exclusion intervals in {@code t};
     * the surviving interval is the intersection of the inclusions. We
     * represent the surviving arc as a start angle and a sweep, using the
     * canonical half-open interval {@code [tStart, tStart + span)}.
     *
     * <p>For the common case where the circle centre is inside the box and
     * {@code rho} is small (full circle inside), we return
     * {@code [-π, π]}.</p>
     *
     * @return {@code {tStart, tEnd}} or {@code null} if no arc survives
     */
    private double[] arcRangeInBox(double cu, double cv, double rho,
                                   double u0, double u1, double v0, double v1) {

        // Build a list of exclusion arcs from each box edge, then subtract
        // them from [−π, π] to find the surviving arc.
        // An exclusion arc for edge u=U (constant) is the range of t where
        // the circle is on the wrong side of that edge.

        // Clip to u >= u0: exclude t where cu + rho*cos(t) < u0
        //   => cos(t) < (u0 - cu)/rho
        // Clip to u <= u1: exclude t where cu + rho*cos(t) > u1
        //   => cos(t) > (u1 - cu)/rho
        // Clip to v >= v0: exclude t where cv + rho*sin(t) < v0
        //   => sin(t) < (v0 - cv)/rho
        // Clip to v <= v1: exclude t where cv + rho*sin(t) > v1
        //   => sin(t) > (v1 - cv)/rho

        // We work with a sorted list of (start, end) surviving intervals
        // on [−π, π], initially the whole circle.
        List<double[]> surviving = new ArrayList<>();
        surviving.add(new double[]{-Math.PI, Math.PI});

        // Clip by each of the four edges.
        surviving = clipByCosine(surviving, rho, u0 - cu, true);  // u >= u0
        surviving = clipByCosine(surviving, rho, u1 - cu, false); // u <= u1
        surviving = clipBySine(surviving,   rho, v0 - cv, true);  // v >= v0
        surviving = clipBySine(surviving,   rho, v1 - cv, false); // v <= v1

        if (surviving.isEmpty()) return null;

        // Take the largest surviving arc (there should normally be just one
        // for a well-formed Kiss cell).
        double[] best = surviving.get(0);
        for (double[] seg : surviving) {
            if (seg[1] - seg[0] > best[1] - best[0]) best = seg;
        }
        if (best[1] - best[0] < 1e-10) return null;
        return best;
    }

    /**
     * Clips a list of angular intervals by the constraint
     * {@code rho*cos(t) >= threshold} (if {@code lowerBound}) or
     * {@code rho*cos(t) <= threshold} (if not {@code lowerBound}).
     */
    private static List<double[]> clipByCosine(List<double[]> intervals,
                                               double rho, double threshold,
                                               boolean lowerBound) {
        double ratio = threshold / rho;
        if (ratio <= -1.0) {
            // Constraint always satisfied (lower) or never (upper).
            return lowerBound ? intervals : new ArrayList<>();
        }
        if (ratio >= 1.0) {
            return lowerBound ? new ArrayList<>() : intervals;
        }
        double alpha = Math.acos(ratio); // in [0, π]
        // cos(t) >= ratio  iff  t ∈ [−alpha, alpha]
        // cos(t) <= ratio  iff  t ∈ [−π, −alpha] ∪ [alpha, π]
        double[] keepLow, keepHigh;
        if (lowerBound) {
            keepLow  = new double[]{-alpha, alpha};
            keepHigh = null;
        } else {
            keepLow  = new double[]{-Math.PI, -alpha};
            keepHigh = new double[]{alpha, Math.PI};
        }
        return intersectIntervals(intervals, keepLow, keepHigh);
    }

    /**
     * Clips a list of angular intervals by the constraint
     * {@code rho*sin(t) >= threshold} (if {@code lowerBound}) or
     * {@code rho*sin(t) <= threshold} (if not {@code lowerBound}).
     */
    private static List<double[]> clipBySine(List<double[]> intervals,
                                             double rho, double threshold,
                                             boolean lowerBound) {
        double ratio = threshold / rho;
        if (ratio <= -1.0) {
            return lowerBound ? intervals : new ArrayList<>();
        }
        if (ratio >= 1.0) {
            return lowerBound ? new ArrayList<>() : intervals;
        }
        double alpha = Math.asin(ratio); // in [−π/2, π/2]
        // sin(t) >= ratio  iff  t ∈ [alpha, π − alpha]
        // sin(t) <= ratio  iff  t ∈ [−π, alpha] ∪ [π − alpha, π]
        double[] keepLow, keepHigh;
        if (lowerBound) {
            keepLow  = new double[]{alpha, Math.PI - alpha};
            keepHigh = null;
        } else {
            keepLow  = new double[]{-Math.PI, alpha};
            keepHigh = new double[]{Math.PI - alpha, Math.PI};
        }
        return intersectIntervals(intervals, keepLow, keepHigh);
    }

    /**
     * Intersects each interval in {@code existing} with the union of
     * {@code keep1} (and optionally {@code keep2}), returning the
     * surviving sub-intervals. {@code keep2} may be {@code null}.
     */
    private static List<double[]> intersectIntervals(List<double[]> existing,
                                                     double[] keep1,
                                                     double[] keep2) {
        List<double[]> result = new ArrayList<>();
        for (double[] seg : existing) {
            double[] i1 = intersect(seg, keep1);
            if (i1 != null) result.add(i1);
            if (keep2 != null) {
                double[] i2 = intersect(seg, keep2);
                if (i2 != null) result.add(i2);
            }
        }
        return result;
    }

    /** Returns the intersection of two intervals, or null if empty. */
    private static double[] intersect(double[] a, double[] b) {
        double lo = Math.max(a[0], b[0]);
        double hi = Math.min(a[1], b[1]);
        return lo < hi ? new double[]{lo, hi} : null;
    }

    // -----------------------------------------------------------------------
    // Arc stitching
    // -----------------------------------------------------------------------

    /**
     * Orders multiple arcs into a single closed loop by matching endpoints.
     * Each arc's end point (at {@code tEnd}) should be within a small
     * tolerance of the next arc's start point (at {@code tStart}).
     *
     * <p>Uses a greedy nearest-neighbour search on 3D endpoint proximity.</p>
     *
     * @param raw unordered list of arcs
     * @return ordered list forming a closed loop
     */
    private List<KissArc> stitchArcs(List<KissArc> raw) {
        List<KissArc> ordered = new ArrayList<>();
        List<KissArc> remaining = new ArrayList<>(raw);

        ordered.add(remaining.remove(0));

        while (!remaining.isEmpty()) {
            KissArc last = ordered.get(ordered.size() - 1);
            Point3D.Double tail = last.cartesian(last.tEnd);

            // Find the arc whose start is closest to the current tail.
            int bestIdx = -1;
            boolean bestReversed = false;
            double bestDist = Double.MAX_VALUE;

            for (int i = 0; i < remaining.size(); i++) {
                KissArc candidate = remaining.get(i);
                double dStart = Point3D.Double.distanceSq(tail, candidate.cartesian(candidate.tStart));
                double dEnd   = Point3D.Double.distanceSq(tail, candidate.cartesian(candidate.tEnd));

                if (dStart < bestDist) { bestDist = dStart; bestIdx = i; bestReversed = false; }
                if (dEnd   < bestDist) { bestDist = dEnd;   bestIdx = i; bestReversed = true;  }
            }

            KissArc next = remaining.remove(bestIdx);
            if (bestReversed) {
                // Reverse the arc by swapping tStart and tEnd and negating the
                // cos/sin argument direction. We do this by wrapping in a reversed arc.
                next = reverseArc(next);
            }
            ordered.add(next);
        }
        return ordered;
    }

    /**
     * Returns a new {@link KissArc} that traverses the same circle in the
     * opposite direction (swaps {@code tStart} and {@code tEnd}, adjusts
     * centre and rho to maintain the same point set).
     */
    private static KissArc reverseArc(KissArc arc) {
        // Reversing: new tStart = -arc.tEnd (modulo), new tEnd = -arc.tStart.
        // Equivalently, remap t' = arc.tStart + arc.tEnd - t, but the simplest
        // correct approach is to negate the parameterisation by replacing t with
        // (tStart + tEnd - t). Since KissArc.cartesian uses cos(t) and sin(t)
        // directly we create a subclass-free wrapper via a fresh KissArc with
        // swapped bounds and a negated centre offset, which is equivalent to
        // reflecting the circle. Instead we just swap tStart/tEnd and let
        // sampleBoundary iterate from tEnd down to tStart.
        //
        // The cleanest option: return a KissArc with tStart = -arc.tEnd and
        // tEnd = -arc.tStart, cu negated in the cos direction, which is
        // achieved by shifting the angle base by π and swapping sign.
        // Simpler: reparameterise as t -> (arc.tStart + arc.tEnd) - t.
        // We capture this by storing a flag — but KissArc is immutable.
        // Since all callers only use tStart/tEnd for iteration bounds and
        // cartesian(t) for point evaluation, and since the arc is symmetric
        // under t -> (S+E-t), just swap the bounds: cartesian will still
        // produce the same circle, just iterated in the other direction.
        return new KissArc(arc.faceAxis, arc.fixedVal, arc.centerU, arc.centerV,
                           arc.rho, arc.tEnd, arc.tStart, arc.R) {
            @Override
            public Point3D.Double cartesian(double t) {
                // Remap so iteration from tStart(=old tEnd) to tEnd(=old tStart)
                // traverses the arc in reverse.
                double tReversed = arc.tStart + arc.tEnd - t;
                return super.cartesian(tReversed);
            }
        };
    }

    // -----------------------------------------------------------------------
    // Utility
    // -----------------------------------------------------------------------

    private static double clamp(double value, double lo, double hi) {
        return Math.max(lo, Math.min(hi, value));
    }


    /**
     * Returns the geographic center of the Kiss region as
     * {@code [latitude, longitude]} in radians.
     *
     * <p>The center is the point on the sphere surface that is closest to the
     * cell — i.e. the sphere surface point in the direction of the closest
     * point on any penetrating face to the origin. For each face whose plane
     * intersects the sphere, the closest point on the face to the origin is
     * found (origin projected onto the face plane, clamped to the face bounds),
     * and the overall closest such point across all penetrating faces is
     * projected outward onto the sphere. That 3D point is then converted to
     * latitude and longitude.</p>
     *
     * <ul>
     *   <li>latitude  = π/2 − θ  (positive toward +Z, i.e. north)</li>
     *   <li>longitude = φ = atan2(y, x), normalised to [−π, π]</li>
     * </ul>
     *
     * @return {@code double[]{latitude, longitude}} in radians
     */
    public double[] getCenter() {
        // Find the closest point on any penetrating face to the origin.
        // For each face: fix one coord, clamp the other two to face bounds.
        double bestDist2 = Double.MAX_VALUE;
        double bestX = 0, bestY = 0, bestZ = 0;

        // z-faces (faceAxis=2): fixed=z, u=x, v=y
        for (double fz : new double[]{z0, z1}) {
            if (fz * fz < R * R) {
                double cx = clamp(0.0, x0, x1);
                double cy = clamp(0.0, y0, y1);
                double d2 = fz*fz + cx*cx + cy*cy;
                if (d2 < bestDist2) { bestDist2 = d2; bestX = cx; bestY = cy; bestZ = fz; }
            }
        }
        // y-faces (faceAxis=1): fixed=y, u=x, v=z
        for (double fy : new double[]{y0, y1}) {
            if (fy * fy < R * R) {
                double cx = clamp(0.0, x0, x1);
                double cz = clamp(0.0, z0, z1);
                double d2 = fy*fy + cx*cx + cz*cz;
                if (d2 < bestDist2) { bestDist2 = d2; bestX = cx; bestY = fy; bestZ = cz; }
            }
        }
        // x-faces (faceAxis=0): fixed=x, u=y, v=z
        for (double fx : new double[]{x0, x1}) {
            if (fx * fx < R * R) {
                double cy = clamp(0.0, y0, y1);
                double cz = clamp(0.0, z0, z1);
                double d2 = fx*fx + cy*cy + cz*cz;
                if (d2 < bestDist2) { bestDist2 = d2; bestX = fx; bestY = cy; bestZ = cz; }
            }
        }

        // Project the closest point outward onto the sphere surface.
        double len = Math.sqrt(bestDist2);
        double sx = bestX / len * R;
        double sy = bestY / len * R;
        double sz = bestZ / len * R;

        // Convert to lat/lon: latitude = pi/2 - theta = asin(z/R), lon = atan2(y,x).
        double lat = Math.asin(sz / R);
        double lon = Math.atan2(sy, sx);

        return new double[]{lat, lon};
    }

    // -----------------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------------

    /**
     * Prints a human-readable summary of the arcs and the normalized area.
     */
    public void printSummary() {
        System.out.printf("KissGeometry for cell [%d,%d,%d]:%n",
                cell.nx, cell.ny, cell.nz);
        System.out.printf("  Bounds: x=[%.4f,%.4f] y=[%.4f,%.4f] z=[%.4f,%.4f]%n",
                x0, x1, y0, y1, z0, z1);
        System.out.printf("  Arcs: %d%n", arcs.size());
        for (int i = 0; i < arcs.size(); i++) {
            System.out.printf("    [%d] %s%n", i, arcs.get(i));
        }
        System.out.printf("  Normalized area A/(4πR²): %.6e%n", normalizedArea());
    }
}