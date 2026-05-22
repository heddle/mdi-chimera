package edu.cnu.mdi.chimera.curve;

import org.apache.commons.math3.analysis.UnivariateFunction;

import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A curve of constant azimuthal angle φ on the sphere surface.
 *
 * <p>A PHI curve lies along a meridian (great-circle arc of constant
 * longitude). Its parameterization is:</p>
 * <pre>
 *   θ(t) = θ₀ + t·Δθ,    t ∈ [0, 1]
 *   φ(t) = φ*             (constant)
 * </pre>
 * <p>where Δθ = θ₁ − θ₀. No branch-cut handling is needed for θ since it
 * is defined on [0, π] without ambiguity.</p>
 *
 * <p>The arc-length formula reduces to the closed form R·|Δθ|.</p>
 */
public class PhiCurve extends BaseCurve {

    /** The constant azimuthal angle. */
    private final double phiStar;

    /** Starting polar angle. */
    private final double theta0;

    /** Signed polar-angle sweep Δθ = θ₁ − θ₀. */
    private final double deltaTheta;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * Constructs a constant-φ curve between two sphere-surface points.
     *
     * <p>Both endpoints should have (approximately) the same azimuthal angle.
     * If they differ by more than {@link #TOL} a warning is printed and the
     * average φ is used.</p>
     *
     * @param p0 start point on the sphere
     * @param p1 end point on the sphere
     * @param r  sphere radius
     */
    public PhiCurve(Point3D.Double p0, Point3D.Double p1, double r) {
        super(p0, p1, r);

        double ph0 = sv0.phi;
        double ph1 = sv1.phi;

        if (Math.abs(MathUtil.normalizeAngle(ph0 - ph1)) > TOL) {
            System.err.printf(
                "[PhiCurve] Endpoints have differing phi values: %.6f vs %.6f%n",
                ph0, ph1);
        }

        phiStar    = MathUtil.normalizeAngle(0.5 * (ph0 + ph1));
        theta0     = sv0.theta;
        deltaTheta = sv1.theta - theta0;
    }

    /**
     * Private constructor used by {@link #reverse()}.
     */
    private PhiCurve(Point3D.Double p0, Point3D.Double p1, double r,
                     double phiStar, double theta0, double deltaTheta) {
        super(p0, p1, r);
        this.phiStar    = phiStar;
        this.theta0     = theta0;
        this.deltaTheta = deltaTheta;
    }

    // -----------------------------------------------------------------------
    // BaseCurve implementation
    // -----------------------------------------------------------------------

    @Override
    public UnivariateFunction getThetaFunction() {
        return t -> theta0 + t * deltaTheta;
    }

    @Override
    public UnivariateFunction getPhiFunction() {
        return t -> phiStar;
    }

    @Override
    public PhiCurve reverse() {
        return new PhiCurve(p1, p0, radius,
                            phiStar,
                            theta0 + deltaTheta,
                            -deltaTheta);
    }

    /**
     * Exact closed-form arc length: R·|Δθ|.
     *
     * @return arc length in the same units as {@code radius}
     */
    @Override
    public double arcLength() {
        return radius * Math.abs(deltaTheta);
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /** @return the constant azimuthal angle φ* in radians */
    public double getPhiStar()    { return phiStar; }

    /** @return the signed polar-angle sweep Δθ in radians */
    public double getDeltaTheta() { return deltaTheta; }

    @Override
    public String toString() {
        return String.format("PhiCurve[phi=%.4f theta0=%.4f theta1=%.4f]",
                phiStar, theta0, theta0 + deltaTheta);
    }
}
