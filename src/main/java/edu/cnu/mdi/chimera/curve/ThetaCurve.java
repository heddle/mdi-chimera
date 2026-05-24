package edu.cnu.mdi.chimera.curve;

import org.apache.commons.math3.analysis.UnivariateFunction;

import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A curve of constant polar angle θ on the sphere surface.
 *
 * <p>A THETA curve lies along a small circle of constant colatitude. Its
 * parameterization is:</p>
 * <pre>
 *   θ(t) = θ*             (constant)
 *   φ(t) = φ₀ + t·Δφ,    t ∈ [0, 1]
 * </pre>
 * <p>where Δφ = φ₁ − φ₀ is normalized to (−π, π] to take the short arc.</p>
 *
 * <p>The arc-length formula reduces to the closed form R·sin(θ*)·|Δφ|,
 * which overrides the numerical integration in {@link BaseCurve}.</p>
 */
public class ThetaCurve extends BaseCurve {

    /** The constant polar angle. */
    private final double thetaStar;

    /** The signed azimuthal sweep Δφ = φ₁ − φ₀, normalised to (−π, π]. */
    private final double deltaPhi;

    /** Starting azimuthal angle. */
    private final double phi0;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * Constructs a constant-θ curve between two sphere-surface points.
     *
     * <p>Both endpoints must have (approximately) the same polar angle; if
     * they differ by more than {@link #TOL} a warning is logged but no
     * exception is thrown, and the average θ is used.</p>
     *
     * @param p0 start point on the sphere
     * @param p1 end point on the sphere
     * @param r  sphere radius
     */
    public ThetaCurve(Point3D.Double p0, Point3D.Double p1, double r) {
        super(p0, p1, r);

        double th0 = sv0.theta;
        double th1 = sv1.theta;

        if (Math.abs(th0 - th1) > TOL) {
            System.err.printf(
                "[ThetaCurve] Endpoints have differing theta values: %.6f vs %.6f%n",
                th0, th1);
        }

        thetaStar = 0.5 * (th0 + th1);
        phi0      = sv0.phi;
        deltaPhi  = MathUtil.normalizeAngle(sv1.phi - phi0);
    }

    /**
     * Private constructor used by {@link #reverse()}.
     */
    private ThetaCurve(Point3D.Double p0, Point3D.Double p1, double r,
                       double thetaStar, double phi0, double deltaPhi) {
        super(p0, p1, r);
        this.thetaStar = thetaStar;
        this.phi0      = phi0;
        this.deltaPhi  = deltaPhi;
    }

    // -----------------------------------------------------------------------
    // BaseCurve implementation
    // -----------------------------------------------------------------------

    @Override
    public UnivariateFunction getThetaFunction() {
        return t -> thetaStar;
    }

    @Override
    public UnivariateFunction getPhiFunction() {
        return t -> MathUtil.normalizeAngle(phi0 + t * deltaPhi);
    }

    @Override
    public ThetaCurve reverse() {
        return new ThetaCurve(p1, p0, radius,
                              thetaStar,
                              MathUtil.normalizeAngle(phi0 + deltaPhi),
                              -deltaPhi);
    }
    
	@Override
	public boolean isConstantTheta() {
		return true;
	}

    /**
     * Exact closed-form arc length: R·sin(θ*)·|Δφ|.
     *
     * @return arc length in the same units as {@code radius}
     */
    @Override
    public double arcLength() {
        return radius * Math.sin(thetaStar) * Math.abs(deltaPhi);
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /** @return the constant polar angle θ* in radians */
    public double getThetaStar() { return thetaStar; }

    /** @return the signed azimuthal sweep Δφ in radians */
    public double getDeltaPhi()  { return deltaPhi; }

    @Override
    public String toString() {
        return String.format("ThetaCurve[theta=%.4f phi0=%.4f phi1=%.4f]",
                thetaStar, phi0, MathUtil.normalizeAngle(phi0 + deltaPhi));
    }
}
