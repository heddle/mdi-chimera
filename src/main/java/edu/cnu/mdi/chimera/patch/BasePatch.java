package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.List;

import edu.cnu.mdi.chimera.app.ChimeraApp;
import edu.cnu.mdi.chimera.area.SphericalPolygonArea;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.Crossing;
import edu.cnu.mdi.chimera.curve.BaseCurve.PoleStatus;
import edu.cnu.mdi.chimera.curve.CompositeCurve;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;
import edu.cnu.mdi.chimera.util.ThetaPhi;

/**
 * Abstract base class for all Mosaic patches.
 *
 * <p>A patch is a region on the sphere surface enclosed by a closed loop of
 * {@link BaseCurve}s, uniquely identified by a 5-tuple {@code (nx, ny, nz)}
 * from the Cartesian grid and {@code (nTheta, nPhi)} from the spherical grid.
 * Prepatches carry only the Cartesian 3-tuple until the splice steps assign
 * the spherical indices.</p>
 *
 * <h2>Area</h2>
 * <p>Estimated via {@link SphericalPolygonArea} using fan triangulation with
 * the Oosterom-Strackee formula. Vertices are sampled directly as
 * {@link Point3D.Double} Cartesian points from {@link BaseCurve#getPoint},
 * bypassing {@link edu.cnu.mdi.chimera.util.ThetaPhi} entirely to avoid any
 * theta/latitude or degree/radian conversion errors.</p>
 *
 * <h2>Pole enclosure</h2>
 * <p>Computed once and cached using the winding-number algorithm from
 * paper Appendix A, including the inter-curve phi jumps at each junction.</p>
 */
public abstract class BasePatch {

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    protected static final double LOOP_TOL = 1.0e-6;
    protected static final double TOL      = 1.0e-8;

    /** Samples per curve for area estimation. The paper uses n=5. */
    protected static final int DEFAULT_AREA_SAMPLES = 50;

    private static final double WINDING_TOL = 0.01; // accumulated rounding over 100 steps/curve

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    public final List<BaseCurve> curves;
    public final double          radius;
    public final int             nx, ny, nz;
    public final int             nTheta, nPhi;

    private Boolean _enclosesNorthPole = null;
    private Boolean _enclosesSouthPole = null;
    
    // A single curve made from the set of BaseCurves
    public final CompositeCurve compositeCurve;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    public BasePatch(List<BaseCurve> curves,
                     int nx, int ny, int nz,
                     int nTheta, int nPhi) {

        if (curves == null || curves.size() < 2)
            throw new IllegalArgumentException("A patch requires at least two curves.");

        this.curves    = curves;
        this.radius    = ChimeraApp.getInstance().getRadius();
        this.nx = nx;  this.ny = ny;  this.nz = nz;
        this.nTheta = nTheta;  this.nPhi = nPhi;

        if (!validateLoop(curves)) {
            throw new IllegalArgumentException(
                "Curves do not form a closed loop for patch (" +
                nx+","+ny+","+nz+","+nTheta+","+nPhi+").");
        }
        this.compositeCurve = new CompositeCurve(curves);
    }

    /** Convenience constructor for prepatches (nTheta = nPhi = -1). */
    public BasePatch(List<BaseCurve> curves, int nx, int ny, int nz) {
        this(curves, nx, ny, nz, -1, -1);
    }

    // -----------------------------------------------------------------------
    // Abstract interface
    // -----------------------------------------------------------------------

    public abstract boolean containsPoint(double x, double y, double z);

    // -----------------------------------------------------------------------
    // Loop validation
    // -----------------------------------------------------------------------

    public static boolean validateLoop(List<BaseCurve> curves) {
        int n = curves.size();
        for (int i = 0; i < n; i++) {
            Point3D.Double end   = curves.get(i).p1;
            Point3D.Double start = curves.get((i + 1) % n).p0;
            double dist = Point3D.Double.distance(end, start);
            if (dist > LOOP_TOL) {
                System.err.printf("[BasePatch] Loop gap at curve %d→%d: distance=%.3e%n",
                    i, (i + 1) % n, dist);
                return false;
            }
        }
        return true;
    }

    public static boolean pointsAreClose(Point3D.Double p1, Point3D.Double p2) {
        return Point3D.Double.distance(p1, p2) < TOL;
    }

    // -----------------------------------------------------------------------
    // Pole enclosure  (paper Appendix A)
    // -----------------------------------------------------------------------

    public boolean enclosesNorthPole() {
        if (_enclosesNorthPole == null) computePoleEnclosure();
        return _enclosesNorthPole;
    }

    public boolean enclosesSouthPole() {
        if (_enclosesSouthPole == null) computePoleEnclosure();
        return _enclosesSouthPole;
    }

    public boolean polar() {
        return enclosesNorthPole() || enclosesSouthPole();
    }

    private void computePoleEnclosure() {
        double totalWinding = 0.0;
        double thetaSum     = 0.0;
        int    thetaCount   = 0;

        int n = curves.size();
        for (int i = 0; i < n; i++) {
            BaseCurve curve = curves.get(i);

            PoleStatus onCurve = curve.poleOnCurve();
            if (onCurve == PoleStatus.NORTH_ON_CURVE) {
                _enclosesNorthPole = true;  _enclosesSouthPole = false;  return;
            }
            if (onCurve == PoleStatus.SOUTH_ON_CURVE) {
                _enclosesNorthPole = false;  _enclosesSouthPole = true;  return;
            }

            totalWinding += curve.windingContribution();

            // Inter-curve phi jump at each junction.
            BaseCurve next   = curves.get((i + 1) % n);
            double phiEnd    = MathUtil.normalizeAngle(curve.sv1.phi);
            double phiStart  = MathUtil.normalizeAngle(next.sv0.phi);
            totalWinding    += MathUtil.normalizeAngle(phiStart - phiEnd);

            thetaSum += curve.averageTheta();
            thetaCount++;
        }

        double avgTheta = (thetaCount > 0) ? thetaSum / thetaCount : Math.PI / 2;

        if (Math.abs(Math.abs(totalWinding) - 2.0 * Math.PI) < WINDING_TOL) {
            _enclosesNorthPole = avgTheta < Math.PI / 2;
            _enclosesSouthPole = !_enclosesNorthPole;
        } else {
            _enclosesNorthPole = false;
            _enclosesSouthPole = false;
        }
    }

    // -----------------------------------------------------------------------
    // Geometry
    // -----------------------------------------------------------------------


    /**
     * Returns a sampled list of boundary vertices as {@link ThetaPhi} points,
     * suitable for drawing.
     *
     * <p>Delegates to {@link #getBoundaryPoints(int)} and converts each
     * Cartesian point to spherical coordinates.</p>
     *
     * @param n samples per curve (≥ 1)
     * @return ordered list of {@link ThetaPhi} boundary points
     */
    public List<edu.cnu.mdi.chimera.util.ThetaPhi> getSphericalVertices(int n) {
        List<Point3D.Double> pts = getBoundaryPoints(n);
        List<edu.cnu.mdi.chimera.util.ThetaPhi> result = new ArrayList<>(pts.size());
        for (Point3D.Double p : pts) {
            edu.cnu.mdi.chimera.util.SphericalVector sv = new edu.cnu.mdi.chimera.util.SphericalVector(p);
            result.add(new edu.cnu.mdi.chimera.util.ThetaPhi(radius, sv.theta, sv.phi));
        }
        return result;
    }

    /**
     * Returns boundary vertices using the default sample count
     * ({@value #DEFAULT_AREA_SAMPLES} samples per curve).
     *
     * @return ordered list of {@link edu.cnu.mdi.chimera.util.ThetaPhi} boundary points
     */
    public List<ThetaPhi> getSphericalVertices() {
        return getSphericalVertices(DEFAULT_AREA_SAMPLES);
    }

    /**
     * Samples the boundary as Cartesian {@link Point3D.Double} points on the
     * sphere surface — one point per step per curve, at
     * {@code t = step, 2·step, ..., 1} (i.e. i = 1..n inclusive).
     *
     * <p>Starting at {@code t = step} rather than {@code t = 0} avoids
     * duplicating the junction point shared with the previous curve's end.
     * Each junction appears exactly once, as the {@code t = 1} sample of
     * the curve that ends there.</p>
     *
     * @param n samples per curve (≥ 1)
     * @return ordered list of boundary points
     */
    public List<Point3D.Double> getBoundaryPoints(int n) {
        if (n < 1) throw new IllegalArgumentException("n must be at least 1.");
        List<Point3D.Double> pts = new ArrayList<>(curves.size() * n);
        double step = 1.0 / n;
        for (BaseCurve curve : curves) {
            for (int i = 1; i <= n; i++) {
                pts.add(curve.getPoint(i * step));
            }
        }
        return pts;
    }

    /**
     * Estimates the normalized area of this patch as a fraction of {@code 4πR²}.
     *
     * @param n samples per curve
     * @return normalized area in {@code [0, 1]}
     */
    public double areaEstimate(int n) {
        List<Point3D.Double> pts = getBoundaryPoints(n);
        return SphericalPolygonArea.computeAreaFraction(pts);
    }

    /** Area estimate using {@value #DEFAULT_AREA_SAMPLES} samples per curve. */
    public double areaEstimate() {
        return areaEstimate(DEFAULT_AREA_SAMPLES);
    }

    public double perimeter() {
        double total = 0.0;
        for (BaseCurve curve : curves) total += curve.arcLength();
        return total;
    }

    // -----------------------------------------------------------------------
    // Index accessors
    // -----------------------------------------------------------------------

    public boolean isFullyIndexed()  { return nTheta >= 0 && nPhi >= 0; }
    public String  cartesianIndex()  { return String.format("(%d,%d,%d)", nx, ny, nz); }

    public String fullIndex() {
        return cartesianIndex() +
               (isFullyIndexed() ? String.format("(%d,%d)", nTheta, nPhi) : "(?,?)");
    }

    @Override
    public String toString() {
        return String.format("%s[%s curves=%d polar=N%s S%s]",
            getClass().getSimpleName(), fullIndex(), curves.size(),
            enclosesNorthPole() ? "✓" : "✗",
            enclosesSouthPole() ? "✓" : "✗");
    }
    

}
