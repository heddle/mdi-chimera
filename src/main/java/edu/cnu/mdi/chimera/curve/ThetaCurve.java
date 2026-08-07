package edu.cnu.mdi.chimera.curve;

import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.ArrayList;
import java.util.List;

import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridContext;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;
import edu.cnu.mdi.chimera.util.SphericalVector;

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
    
    /** Cached list of crossings where the curve crosses phi grid lines. */
    private List<Crossing> phiCrossings;

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
    ThetaCurve(Point3D.Double p0, Point3D.Double p1, double r,
                       double thetaStar, double phi0, double deltaPhi) {
        super(p0, p1, r);
        this.thetaStar = thetaStar;
        this.phi0      = phi0;
        this.deltaPhi  = deltaPhi;
    }

    /**
     * Reconstructs an explicitly parameterized directed constant-theta arc.
     *
     * <p>This factory is intended for lossless interchange-file import. Unlike
     * the ordinary constructor it does not infer or normalize the angular sweep,
     * so complementary long arcs retain their exported direction and extent.</p>
     *
     * @param p0 start point on the sphere
     * @param p1 end point on the sphere
     * @param r sphere radius
     * @param thetaStar constant colatitude in radians
     * @param phi0 starting azimuth in radians
     * @param deltaPhi signed azimuthal sweep in radians
     * @return the reconstructed directed arc
     */
    public static ThetaCurve parameterized(Point3D.Double p0, Point3D.Double p1,
            double r, double thetaStar, double phi0, double deltaPhi) {
        return new ThetaCurve(p0, p1, r, thetaStar, phi0, deltaPhi);
    }
    
    /**
     * Constructs a constant-theta curve between two sphere-surface points, allowing
     * the caller to choose either the normal short azimuthal arc or the complementary
     * long arc.
     *
     * <p>
     * The ordinary {@link #ThetaCurve(Point3D.Double, Point3D.Double, double)}
     * constructor always takes the normalized short arc in phi. That is correct for
     * most theta connectors, but it is ambiguous near poles where phi is singular.
     * In those polar cases, the long arc may be the geometrically correct boundary
     * of the small polar patch.
     * </p>
     *
     * @param p0      start point
     * @param p1      end point
     * @param r       sphere radius
     * @param longArc if true, use the complementary long phi arc
     * @return a theta curve from {@code p0} to {@code p1}
     */
    public static ThetaCurve between(Point3D.Double p0,
                                     Point3D.Double p1,
                                     double r,
                                     boolean longArc) {

        SphericalVector sv0 = new SphericalVector(p0);
        SphericalVector sv1 = new SphericalVector(p1);

        double thetaStar = 0.5 * (sv0.theta + sv1.theta);
        double phi0 = sv0.phi;

        double delta = MathUtil.normalizeAngle(sv1.phi - phi0);

        if (longArc) {
            /*
             * Replace the normalized short arc with its complementary arc.
             * Example:
             *
             *   +30 degrees  -> -330 degrees
             *   -30 degrees  -> +330 degrees
             *
             * Reversing a short curve only changes traversal direction; it does not
             * create this complementary geometry.
             */
            if (delta >= 0.0) {
                delta -= 2.0 * Math.PI;
            } else {
                delta += 2.0 * Math.PI;
            }
        }

        return new ThetaCurve(p0, p1, r, thetaStar, phi0, delta);
    }
    
    /**
     * Computes candidate crossings where this constant-theta curve meets phi grid
     * lines.
     *
     * <p>
     * Unlike {@link #getPhiFunction()}, this method uses the unwrapped
     * parameterization {@code phi0 + t * deltaPhi}. That is essential because
     * {@code deltaPhi} may represent either the short arc or the complementary long
     * arc for polar theta connectors.
     * </p>
     *
     * <p>
     * Endpoint hits are included as candidate crossings. The patch-level phi splicer
     * can decide later whether an endpoint hit is a true splice crossing or merely a
     * junction touch.
     * </p>
     *
     * @return list of candidate crossings with phi grid lines
     */
    @Override
    public List<Crossing> getPhiCrossings() {
        if (phiCrossings != null) {
            return phiCrossings;
        }

        phiCrossings = new ArrayList<>();

        /*
         * At the poles, phi is singular. A constant-theta curve at theta=0 or pi
         * should not report crossings with every meridian.
         */
        if (Math.sin(thetaStar) < TOL) {
            return phiCrossings;
        }

        if (Math.abs(deltaPhi) < TOL) {
            return phiCrossings;
        }

        SphericalGrid grid = ChimeraGridContext.sphericalGrid();
        Grid1D phiGrid = grid.getPhiGrid();

        double start = phi0;
        double end = phi0 + deltaPhi;

        double lo = Math.min(start, end);
        double hi = Math.max(start, end);

        double twoPi = 2.0 * Math.PI;
        double phiTol = Math.max(1.0e-12, 1.0e-10 * Math.abs(deltaPhi));

        for (int i = 0; i < phiGrid.numPoints(); i++) {
            double targetPhi = phiGrid.valueAt(i);

            int kMin = (int) Math.floor((lo - targetPhi) / twoPi) - 1;
            int kMax = (int) Math.ceil((hi - targetPhi) / twoPi) + 1;

            for (int k = kMin; k <= kMax; k++) {
                double phi = targetPhi + k * twoPi;

                if (phi < lo - phiTol || phi > hi + phiTol) {
                    continue;
                }

                double t = (phi - start) / deltaPhi;

                if (t < -1.0e-10 || t > 1.0 + 1.0e-10) {
                    continue;
                }

                t = Math.max(0.0, Math.min(1.0, t));

                /*
                 * Keep endpoint policy explicit. Endpoint candidate crossings are
                 * added below, so skip them here.
                 */
                if (isEndpointT(t)) {
                    continue;
                }

                phiCrossings.add(new Crossing(this, t,
                        MathUtil.normalizeAngle(targetPhi), i));
            }
        }

        addEndpointPhiCrossing(phiCrossings, phiGrid, 0.0);
        addEndpointPhiCrossing(phiCrossings, phiGrid, 1.0);

        phiCrossings = Crossing.removeDuplicates(phiCrossings);
        return phiCrossings;
    }
    
	@Override
	public List<ThetaCurve> splitAtPhiCrossings() {
		List<Crossing> crossings = getPhiCrossings();
		if (crossings.isEmpty()) {
			return List.of(this);
		}
		
	    List<Double> splitTs = new ArrayList<>();
	    for (Crossing crossing : crossings) {
	        if (crossing == null) {
	            continue;
	        }

	        double t = crossing.t();

	        if (!Double.isFinite(t)) {
	            continue;
	        }

	        if (isEndpointT(t)) {
	            continue;
	        }

	        t = Math.max(0.0, Math.min(1.0, t));	        
	        boolean duplicate = false;
	        for (double existing : splitTs) {
	            if (Math.abs(existing - t) < 1.0e-10) {
	                duplicate = true;
	                break;
	            }
	        }

	        if (!duplicate) {
	            splitTs.add(t);
	        }
	        
	    } //crossing loop
	    
	    if (splitTs.isEmpty()) {
	        return List.of(this);
	    }

	    splitTs.sort(Double::compare);

	    List<ThetaCurve> pieces = new ArrayList<>();
	    
	    double t0 = 0.0;
	    for (double t1 : splitTs) {
	        if (t1 - t0 > TOL) {
	            pieces.add(subCurve(t0, t1));
	        }
	        t0 = t1;
	    }

	    if (1.0 - t0 > TOL) {
	        pieces.add(subCurve(t0, 1.0));
	    }

	    if (pieces.isEmpty()) {
	        return List.of(this);
	    }

	    return pieces;
	}
	
	private ThetaCurve subCurve(double t0, double t1) {
	    double phiStart = phi0 + t0 * deltaPhi;
	    double phiEnd = phi0 + t1 * deltaPhi;
	    double delta = MathUtil.normalizeAngle(phiEnd - phiStart);
	    return new ThetaCurve(getPointAt(t0), getPointAt(t1), radius, thetaStar, phiStart, delta);
	}
	
	private Point3D.Double getPointAt(double t) {
	    double phi = phi0 + t * deltaPhi;
	    SphericalVector sv = new SphericalVector(thetaStar, phi, radius);
	    return sv.toCartesian();
	}

    /**
     * Adds a candidate endpoint phi crossing if the endpoint lies on a phi grid
     * line and the curve immediately leaves that line.
     */
    private void addEndpointPhiCrossing(List<Crossing> crossings,
                                        Grid1D phiGrid,
                                        double tEndpoint) {

        if (Math.sin(theta(tEndpoint)) < TOL) {
            return;
        }

        double ph = MathUtil.normalizeAngle(phi(tEndpoint));

        for (int i = 0; i < phiGrid.numPoints(); i++) {
            double targetPhi = MathUtil.normalizeAngle(phiGrid.valueAt(i));

            if (Math.abs(MathUtil.normalizeAngle(ph - targetPhi)) > TOL) {
                continue;
            }

            double probeT = (tEndpoint <= 0.5) ? 1.0e-6 : 1.0 - 1.0e-6;
            probeT = Math.max(0.0, Math.min(1.0, probeT));

            double probePhi = MathUtil.normalizeAngle(phi(probeT));

            if (Math.abs(MathUtil.normalizeAngle(probePhi - targetPhi)) < TOL) {
                continue;
            }

            crossings.add(new Crossing(this, tEndpoint, targetPhi, i));
        }
    }

    /**
     * Tests whether {@code t} is effectively an endpoint of the curve.
     */
    private static boolean isEndpointT(double t) {
        return t <= 1.0e-8 || t >= 1.0 - 1.0e-8;
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

    /** @return the starting azimuth φ₀ in radians */
    public double getPhi0() { return phi0; }

    /** @return the signed azimuthal sweep Δφ in radians */
    public double getDeltaPhi()  { return deltaPhi; }

    @Override
    public String toString() {
        return String.format("ThetaCurve[theta=%.4f phi0=%.4f phi1=%.4f]",
                thetaStar, phi0, MathUtil.normalizeAngle(phi0 + deltaPhi));
    }
}
