package edu.cnu.mdi.chimera.curve;

import org.apache.commons.math3.analysis.UnivariateFunction;

import edu.cnu.mdi.chimera.util.ChimeraPlane;
import edu.cnu.mdi.chimera.util.ChimeraRotation;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A GENERAL curve: the intersection of the sphere surface with one face of a
 * rectangular grid cell (paper §5).
 *
 * <p>A GENERAL curve lies simultaneously on the sphere and exactly one
 * cell face. Its endpoints are the two points where cell edges pierce the
 * sphere surface. The parameterization is constructed via the paper's
 * rotation algorithm:</p>
 * <ol>
 *   <li>Build the face plane P from three face corners and obtain its
 *       normal n̂ and rotation matrix R that maps ẑ → n̂.</li>
 *   <li>Rotate both endpoints into the primed frame. In that frame they
 *       share a common z' = const, hence a common θ* = θ'.</li>
 *   <li>Read off φ'₀ and φ'₁; normalize Δφ' to (−π, π] to take the
 *       shorter arc (branch-cut fix).</li>
 *   <li>Parameterize in the primed frame: θ'(t) = θ*, φ'(t) = φ'₀ + t·Δφ'.</li>
 *   <li>Rotate each point back with R⁻¹ to recover (R, θ(t), φ(t)) in the
 *       original frame.</li>
 * </ol>
 *
 * <p>Steps 1–3 are performed once at construction time; the returned
 * {@link UnivariateFunction}s for θ and φ apply step 5 on every evaluation,
 * which is O(1) per sample (one matrix–vector multiply and two trig
 * calls).</p>
 */
public class GeneralCurve extends BaseCurve {

    /**
     * The inverse rotation matrix (3×3) that maps the primed frame back to
     * the original frame. Stored as a flat array for efficiency.
     */
    private final double[][] invMatrix;

    /**
     * The constant polar angle θ* in the primed (rotated) frame.
     * Both endpoints share this value by construction.
     */
    private final double thetaStar;

    /** φ'₀: azimuthal angle of the start point in the primed frame. */
    private final double primePhi0;

    /** Δφ': signed azimuthal sweep in the primed frame, normalised to (−π, π]. */
    private final double deltaPrimePhi;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * Constructs a GENERAL curve using the paper's rotation algorithm.
     *
     * @param p0      start point on the sphere (edge–sphere intersection)
     * @param p1      end point on the sphere (edge–sphere intersection)
     * @param r       sphere radius
     * @param corner0 any corner of the cell face containing this curve
     * @param corner1 a second corner of the same face
     * @param corner2 a third corner of the same face
     * @throws IllegalArgumentException if the three corners are collinear
     */
    public GeneralCurve(Point3D.Double p0, Point3D.Double p1, double r,
                        Point3D.Double corner0,
                        Point3D.Double corner1,
                        Point3D.Double corner2) {
        super(p0, p1, r);

        // ------------------------------------------------------------------
        // Steps 3–5 (paper): build plane and rotation matrix R.
        // ChimeraPlane computes the face normal and the rotation that maps
        // ẑ onto that normal (equivalently, maps the normal to ẑ in the
        // inverse direction).
        // ------------------------------------------------------------------
        ChimeraPlane plane = new ChimeraPlane(corner0, corner1, corner2);

        // rmat.matrix maps original → primed (ẑ aligns with face normal).
        // rmat.invMatrix maps primed → original.
        double[][] rotMatrix = plane.rmat.matrix;
        invMatrix            = plane.rmat.invMatrix;

        // ------------------------------------------------------------------
        // Steps 6–7 (paper): rotate both endpoints into the primed frame
        // and read off θ* and φ'₀, φ'₁.
        // ------------------------------------------------------------------
        double[] v0prime = ChimeraRotation.multiplyMatrixVector(rotMatrix,
                new double[]{p0.x, p0.y, p0.z});
        double[] v1prime = ChimeraRotation.multiplyMatrixVector(rotMatrix,
                new double[]{p1.x, p1.y, p1.z});

        // Both rotated points are on the sphere; their z' values should
        // agree. Use their average as θ* for numerical stability.
        double zPrime0 = v0prime[2];
        double zPrime1 = v1prime[2];

        if (Math.abs(zPrime0 - zPrime1) > TOL * r) {
            System.err.printf(
                "[GeneralCurve] z' values differ after rotation: %.6f vs %.6f%n",
                zPrime0, zPrime1);
        }

        double zPrime = 0.5 * (zPrime0 + zPrime1);
        // Clamp to [-R, R] before acos to guard against tiny numerical overruns.
        thetaStar = Math.acos(Math.max(-1.0, Math.min(1.0, zPrime / r)));

        // ------------------------------------------------------------------
        // Step 8 (paper, with branch-cut fix): compute Δφ' and normalise.
        // ------------------------------------------------------------------
        primePhi0    = MathUtil.normalizeAngle(Math.atan2(v0prime[1], v0prime[0]));
        double phi1p = MathUtil.normalizeAngle(Math.atan2(v1prime[1], v1prime[0]));
        deltaPrimePhi = MathUtil.normalizeAngle(phi1p - primePhi0);
    }

    /**
     * Private constructor used by {@link #reverse()}.
     */
    private GeneralCurve(Point3D.Double p0, Point3D.Double p1, double r,
                         double[][] invMatrix, double thetaStar,
                         double primePhi0, double deltaPrimePhi) {
        super(p0, p1, r);
        this.invMatrix      = invMatrix;
        this.thetaStar      = thetaStar;
        this.primePhi0      = primePhi0;
        this.deltaPrimePhi  = deltaPrimePhi;
    }

    // -----------------------------------------------------------------------
    // BaseCurve implementation
    // -----------------------------------------------------------------------

    /**
     * Returns θ(t) by evaluating the curve in the primed frame and rotating
     * back (paper step 9).
     *
     * <p>In the primed frame: θ'(t) = θ*, φ'(t) = φ'₀ + t·Δφ'. The
     * corresponding Cartesian primed point is rotated back by R⁻¹ and
     * the spherical coordinates of the result give θ(t).</p>
     */
    @Override
    public UnivariateFunction getThetaFunction() {
        return t -> {
            double[] xyz = primeToOriginal(t);
            double r2 = Math.sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
            return Math.acos(Math.max(-1.0, Math.min(1.0, xyz[2] / r2)));
        };
    }

    /**
     * Returns φ(t) by evaluating the curve in the primed frame and rotating
     * back (paper step 9).
     */
    @Override
    public UnivariateFunction getPhiFunction() {
        return t -> {
            double[] xyz = primeToOriginal(t);
            return MathUtil.normalizeAngle(Math.atan2(xyz[1], xyz[0]));
        };
    }

    /**
     * Overrides {@link BaseCurve#getPoint} for efficiency: rather than
     * computing (θ, φ) and re-converting to Cartesian, we rotate directly
     * from the primed frame.
     */
    @Override
    public Point3D.Double getPoint(double t) {
        double[] xyz = primeToOriginal(t);
        return new Point3D.Double(xyz[0], xyz[1], xyz[2]);
    }


    /**
     * Returns a new {@link GeneralCurve} covering the sub-interval
     * {@code [t0, t1]} of this curve's parameter domain {@code [0, 1]}.
     *
     * <p>The subrange curve shares the same face plane (same rotation
     * matrices and θ*) but has its φ' range rescaled to
     * {@code [φ'(t0), φ'(t1)]}. The endpoints are computed exactly via
     * {@link #getPoint(double)}.</p>
     *
     * @param t0 start of the sub-interval, in {@code [0, 1]}
     * @param t1 end   of the sub-interval, in {@code [0, 1]}, {@code > t0}
     * @return the subrange curve
     * @throws IllegalArgumentException if {@code t0 >= t1} or either is
     *                                  outside {@code [0, 1]}
     */
    public GeneralCurve subrange(double t0, double t1) {
        if (t0 < 0.0 || t1 > 1.0 || t0 >= t1) {
            throw new IllegalArgumentException(
                String.format("Invalid subrange [%.4f, %.4f]", t0, t1));
        }
        Point3D.Double newP0    = getPoint(t0);
        Point3D.Double newP1    = getPoint(t1);
        double newPrimePhi0     = MathUtil.normalizeAngle(primePhi0 + t0 * deltaPrimePhi);
        double newDeltaPrimePhi = (t1 - t0) * deltaPrimePhi;
        return new GeneralCurve(newP0, newP1, radius,
                                invMatrix, thetaStar,
                                newPrimePhi0, newDeltaPrimePhi);
    }

    @Override
    public GeneralCurve reverse() {
        double newPhi0     = MathUtil.normalizeAngle(primePhi0 + deltaPrimePhi);
        double newDeltaPhi = -deltaPrimePhi;
        return new GeneralCurve(p1, p0, radius,
                                invMatrix, thetaStar, newPhi0, newDeltaPhi);
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /**
     * Evaluates the curve at {@code t} in the primed frame and rotates back
     * to the original frame, returning a Cartesian xyz triple.
     *
     * <p>In the primed frame the point is the small-circle arc
     * (R, θ*, φ'(t)), which in Cartesian is:</p>
     * <pre>
     *   x' = R·sin(θ*)·cos(φ'(t))
     *   y' = R·sin(θ*)·sin(φ'(t))
     *   z' = R·cos(θ*)
     * </pre>
     */
    private double[] primeToOriginal(double t) {
        double phiPrime = MathUtil.normalizeAngle(primePhi0 + t * deltaPrimePhi);
        double sinTheta = Math.sin(thetaStar);
        double xp = radius * sinTheta * Math.cos(phiPrime);
        double yp = radius * sinTheta * Math.sin(phiPrime);
        double zp = radius * Math.cos(thetaStar);
        return ChimeraRotation.multiplyMatrixVector(invMatrix, new double[]{xp, yp, zp});
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /** @return θ* — the constant polar angle in the primed (rotated) frame */
    public double getThetaStar()     { return thetaStar; }

    /** @return φ'₀ — start azimuthal angle in the primed frame */
    public double getPrimePhi0()     { return primePhi0; }

    /** @return Δφ' — signed azimuthal sweep in the primed frame */
    public double getDeltaPrimePhi() { return deltaPrimePhi; }

    @Override
    public String toString() {
        return String.format(
            "GeneralCurve[thetaStar=%.4f primePhi0=%.4f deltaPhi=%.4f]",
            thetaStar, primePhi0, deltaPrimePhi);
    }
}
