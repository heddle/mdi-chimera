package edu.cnu.mdi.chimera.curve;

import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;

import edu.cnu.mdi.chimera.app.ChimeraApp;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;
import edu.cnu.mdi.chimera.util.SphericalVector;

/**
 * Abstract base class for all Mosaic boundary curves.
 *
 * <p>A curve is a parametrized path (R, θ(t), φ(t)), t ∈ [0, 1] on the
 * surface of a sphere of radius {@code R}. Concrete subclasses implement
 * {@link #getThetaFunction()} and {@link #getPhiFunction()}; everything else
 * is derived from those two functions.</p>
 *
 * <p>Three curve types exist:</p>
 * <ul>
 *   <li>{@link ThetaCurve} — constant polar angle θ</li>
 *   <li>{@link PhiCurve}   — constant azimuthal angle φ</li>
 *   <li>{@link GeneralCurve} — lies simultaneously on the sphere and one
 *       cell face; θ and φ both vary</li>
 * </ul>
 */
public abstract class BaseCurve {

	private static int nextCurveId = 0; // For debugging and logging purposes, assign a unique ID to each curve instance

    /** Tolerance for pole detection and other angular comparisons. */
    protected static final double TOL = 1.0e-8;

    /** Maximum number of root-finding intervals for t parameter computations. */
	protected static final int MAX_ROOT_INTERVALS = 100;

    /** Number of samples used for pole detection and winding-number scans. */
    private static final int POLE_SCAN_STEPS = 100;

    /** Number of uniform sub-intervals used for initial bracketing of θ and φ crossings. */
    protected static final int BRACKET_INTERVALS = 100;

    /** Number of sub-intervals used by the Romberg arc-length integrator. */
    private static final int ARC_LENGTH_MAX_ITER = 32;

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    /** Start point in Cartesian coordinates (on the sphere surface). */
    public final Point3D.Double p0;

    /** End point in Cartesian coordinates (on the sphere surface). */
    public final Point3D.Double p1;

    /** Start point in spherical coordinates. */
    public final SphericalVector sv0;

    /** End point in spherical coordinates. */
    public final SphericalVector sv1;

    /** Sphere radius. */
    public final double radius;

    /** Unique ID for this curve instance, useful for debugging and logging. */
    public final int curveId; // Unique ID for this curve instance, useful for debugging and logging


    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * Constructs a curve between two points on the sphere surface.
     *
     * @param p0 start point (must lie on the sphere within numerical tolerance)
     * @param p1 end point   (must lie on the sphere within numerical tolerance)
     * @param r  sphere radius; must be positive and finite
     * @throws IllegalArgumentException if {@code r} is not positive and finite
     */
    public BaseCurve(Point3D.Double p0, Point3D.Double p1, double r) {
        if (!Double.isFinite(r) || r <= 0.0) {
            throw new IllegalArgumentException("Sphere radius must be positive and finite.");
        }
        if (p0 == null || p1 == null) {
            throw new IllegalArgumentException("Curve endpoints must not be null.");
        }
        this.p0     = p0;
        this.p1     = p1;
        this.radius = r;
        this.sv0    = new SphericalVector(p0);
        this.sv1    = new SphericalVector(p1);
        this.curveId = nextCurveId++; // Assign a unique ID to this curve instance
    }

    // -----------------------------------------------------------------------
    // Abstract interface
    // -----------------------------------------------------------------------

    /**
     * Returns the parametrized polar-angle function θ(t), t ∈ [0, 1].
     *
     * @return θ as a {@link UnivariateFunction}
     */
    public abstract UnivariateFunction getThetaFunction();

    /**
     * Returns the parametrized azimuthal-angle function φ(t), t ∈ [0, 1].
     *
     * @return φ as a {@link UnivariateFunction}
     */
    public abstract UnivariateFunction getPhiFunction();

    /**
     * Returns a copy of this curve with the endpoints reversed.
     * The reversed curve traverses the same path in the opposite direction.
     *
     * @return a new curve from {@code p1} to {@code p0}
     */
    public abstract BaseCurve reverse();

    // -----------------------------------------------------------------------
    // Concrete evaluation methods
    // -----------------------------------------------------------------------

    /**
     * Returns θ(t) in radians.
     *
     * @param t parameter in [0, 1]
     * @return polar angle in radians
     */
    public double theta(double t) {
        return getThetaFunction().value(t);
    }

    /**
     * Returns φ(t) in radians, normalised to [−π, π].
     *
     * @param t parameter in [0, 1]
     * @return azimuthal angle in radians
     */
    public double phi(double t) {
        return getPhiFunction().value(t);
    }

    /**
     * Returns the 3-D Cartesian point on the sphere surface at parameter {@code t}.
     *
     * @param t parameter in [0, 1]
     * @return point on the sphere
     */
    public Point3D.Double getPoint(double t) {
        double th  = theta(t);
        double ph  = phi(t);
        double sin = Math.sin(th);
        return new Point3D.Double(
                radius * sin * Math.cos(ph),
                radius * sin * Math.sin(ph),
                radius * Math.cos(th));
    }

    /**
     * Returns true if this curve has constant θ (i.e., is a {@link ThetaCurve} or
     * was constructed as a GeneralCurve but turned out to have
     * constant theta).
     * @return true if this curve has constant θ, false otherwise
     */
    public abstract boolean isConstantTheta();

    /**
     * Returns the spherical coordinates at parameter {@code t}.
     *
     * @param t parameter in [0, 1]
     * @return a {@link SphericalVector} at radius {@code R}
     */
    public SphericalVector getSphericalVector(double t) {
        return new SphericalVector(theta(t), phi(t), radius);
    }

    // -----------------------------------------------------------------------
    // Pole detection  (paper §5.2.1 and Appendix A)
    // -----------------------------------------------------------------------

    /**
     * Pole enclosure / pole-on-curve return values, matching the paper's
     * Algorithm 1 return codes.
     */
    public enum PoleStatus {
        /** No pole enclosed or on the curve. */
        NONE(0),
        /** North pole (θ = 0) enclosed by the closed loop. */
        NORTH_ENCLOSED(1),
        /** South pole (θ = π) enclosed by the closed loop. */
        SOUTH_ENCLOSED(2),
        /** North pole lies on this curve (within {@link #TOL}). */
        NORTH_ON_CURVE(-1),
        /** South pole lies on this curve (within {@link #TOL}). */
        SOUTH_ON_CURVE(-2);

        public final int code;
        PoleStatus(int code) { this.code = code; }
    }

    /**
     * Checks whether a pole lies on this curve (paper Algorithm 2).
     *
     * <p>Samples θ(t) at {@link #POLE_SCAN_STEPS} uniform steps and returns
     * {@link PoleStatus#NORTH_ON_CURVE} if any sample has θ ≈ 0, or
     * {@link PoleStatus#SOUTH_ON_CURVE} if any sample has θ ≈ π; otherwise
     * {@link PoleStatus#NONE}.</p>
     *
     * @return pole status for this curve in isolation
     */
    public PoleStatus poleOnCurve() {
        UnivariateFunction thetaFn = getThetaFunction();
        double step = 1.0 / POLE_SCAN_STEPS;
        for (int i = 0; i <= POLE_SCAN_STEPS; i++) {
            double t  = i * step;
            double th = thetaFn.value(t);
            if (Math.abs(th) < TOL) {
				return PoleStatus.NORTH_ON_CURVE;
			}
            if (Math.abs(th - Math.PI) < TOL) {
				return PoleStatus.SOUTH_ON_CURVE;
			}
        }
        return PoleStatus.NONE;
    }
    
    /**
	 * Returns a list of crossings where this curve intersects the φ grid lines.
	 *
	 * <p>The default implementation returns {@code null}; only {@link ThetaCurve}
	 * and {@link GeneralCurve} override this method to compute actual crossings.</p>
	 *
	 * @return list of φ crossings, or {@code null} if not applicable
	 */
	public List<Crossing> getPhiCrossings() {
		// default implementation returns null; only ThetaCurve and GeneralCurve override this method to compute actual crossings
		return null;
	}
	
	/**
	 * Returns a list of sub-curves resulting from splitting this curve at φ crossings.
	 *
	 * <p>The default implementation returns {@code null}; only {@link ThetaCurve}
	 * and {@link GeneralCurve} override this method to compute actual splits at φ crossings.</p>
	 *
	 * @return list of sub-curves split at φ crossings, or {@code null} if not applicable
	 */
	public List<? extends BaseCurve> splitAtPhiCrossings() {
		// default implementation returns null; only ThetaCurve and GeneralCurve override this method to compute actual splits at φ crossings
		return null;
	}
	
	
	/**
	 * Returns true if a pole lies on this curve.
	 *
	 *@param point the pole to check (should be either the north pole at θ=0 or the south pole at θ=π)
	 *@param tolerance tolerance 
	 * @return true if a pole lies on this curve, false otherwise
	 */
	public boolean pointOnCurve(SphericalVector point, double tolerance) {
		if (point == null || !Double.isFinite(tolerance) || tolerance < 0.0) {
			return false;
		}

		double targetTheta;
		if (Math.abs(point.theta) <= tolerance) {
			targetTheta = 0.0;
		} else if (Math.abs(point.theta - Math.PI) <= tolerance) {
			targetTheta = Math.PI;
		} else {
			return false;
		}

		double step = 1.0 / POLE_SCAN_STEPS;
		for (int i = 0; i <= POLE_SCAN_STEPS; i++) {
			if (Math.abs(theta(i * step) - targetTheta) <= tolerance) {
				return true;
			}
		}
		return false;
	}
	
	
    /**
     * Computes the winding number contribution of this curve (paper Algorithm 3).
     *
     * <p>Accumulates the signed azimuthal displacement Δφ along the curve,
     * unwrapping each step to (−π, π] to handle the branch cut.</p>
     *
     * @return total signed Δφ in radians for this curve
     */
    public double windingContribution() {
        UnivariateFunction phiFn = getPhiFunction();
        double total   = 0.0;
        double prevPhi = phiFn.value(0.0);
        double step    = 1.0 / POLE_SCAN_STEPS;
        for (int i = 1; i <= POLE_SCAN_STEPS; i++) {
            double currPhi = phiFn.value(i * step);
            total   += MathUtil.normalizeAngle(currPhi - prevPhi);
            prevPhi  = currPhi;
        }
        return total;
    }

    /**
     * Computes the average θ along this curve (paper Algorithm 4).
     *
     * @return mean polar angle in radians
     */
    public double averageTheta() {
        UnivariateFunction thetaFn = getThetaFunction();
        double sum  = 0.0;
        double step = 1.0 / POLE_SCAN_STEPS;
        for (int i = 0; i <= POLE_SCAN_STEPS; i++) {
            sum += thetaFn.value(i * step);
        }
        return sum / (POLE_SCAN_STEPS + 1);
    }

    // -----------------------------------------------------------------------
    // Arc-length  (paper §6.2)
    // -----------------------------------------------------------------------

    /**
     * Computes the arc length of this curve (paper Eq. 4):
     * <pre>
     *   L = R ∫₀¹ √( (dθ/dt)² + sin²θ(t)·(dφ/dt)² ) dt
     * </pre>
     * using the Romberg integrator from Apache Commons Math.
     *
     * @return arc length in the same units as {@code radius}
     */
    public double arcLength() {
        UnivariateFunction integrand = t -> {
            double th    = theta(t);
            double dThdt = dTheta_dt(t);
            double dPhdt = dPhi_dt(t);
            double sinTh = Math.sin(th);
            return radius * Math.sqrt(dThdt * dThdt + sinTh * sinTh * dPhdt * dPhdt);
        };
        RombergIntegrator integrator = new RombergIntegrator();
        return integrator.integrate(
                (int) Math.pow(2, ARC_LENGTH_MAX_ITER), integrand, 0.0, 1.0);
    }

    // -----------------------------------------------------------------------
    // TValue support
    // -----------------------------------------------------------------------

    /**
     * Creates a {@link TValue} for the θ function at the given target.
     *
     * @param targetTheta target polar angle in radians
     * @param tolerance   root-finding tolerance
     * @return a TValue encapsulating the found parameter
     */
    public TValue thetaTValue(double targetTheta, double tolerance) {
        return new TValue(this, getThetaFunction(), targetTheta, 0.0, 1.0, tolerance);
    }

    /**
     * Creates a {@link TValue} for the φ function at the given target.
     *
     * @param targetPhi target azimuthal angle in radians
     * @param tolerance root-finding tolerance
     * @return a TValue encapsulating the found parameter
     */
    public TValue phiTValue(double targetPhi, double tolerance) {
        return new TValue(this, getPhiFunction(), targetPhi, 0.0, 1.0, tolerance);
    }

    // -----------------------------------------------------------------------
    // Private numerical-derivative helpers for arc length
    // -----------------------------------------------------------------------

    /** Centred-difference derivative of θ(t). */
    private double dTheta_dt(double t) {
        double h = 1e-6;
        double lo = Math.max(0.0, t - h);
        double hi = Math.min(1.0, t + h);
        return (theta(hi) - theta(lo)) / (hi - lo);
    }

    /** Centred-difference derivative of φ(t), with angle unwrapping. */
    private double dPhi_dt(double t) {
        double h  = 1e-6;
        double lo = Math.max(0.0, t - h);
        double hi = Math.min(1.0, t + h);
        return MathUtil.normalizeAngle(phi(hi) - phi(lo)) / (hi - lo);
    }

    /**
     * Returns the 0-based index of the θ grid cell containing the midpoint of this curve.
     * @return the 0-based index of the θ grid cell containing the midpoint of this curve
     */
    public int getMidpointThetaIndex() {
    	double midTheta = theta(0.5);
    	SphericalGrid grid = ChimeraApp.getInstance().getSphericalGrid();
    	return grid.getThetaGrid().cellIndex(midTheta);
    }

    public String shortString() {
		return String.format("Curve %d: %s %s", curveId, sv0.toStringDegrees("sv0"), sv1.toStringDegrees("sv1"));
	}



}
