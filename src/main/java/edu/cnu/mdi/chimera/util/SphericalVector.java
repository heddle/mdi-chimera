package edu.cnu.mdi.chimera.util;

public class SphericalVector {

	/** the polar angle in radians */
	public double theta;

	/** the azimuthal angle in radians */
	public double phi;

	/** the radius */
	public double r;

	/**
	 * Create a spherical vector. Azimtuthal angles will be normalized to the range
	 * [-pi, pi].
	 *
	 * @param theta the polar angle in radians
	 * @param phi   the azimuthal angle in radians
	 * @param r     the radius
	 */
	public SphericalVector(double theta, double phi, double r) {
		this.theta = theta;
		this.phi = MathUtil.normalizeAngle(phi);
		this.r = r;
	}
	
	public String toString() {
		double theta = Math.toDegrees(this.theta);
		double phi = Math.toDegrees(this.phi);
		return String.format("theta: %.3f phi %.3f", theta, phi);
	}

	/**
     * Create a spherical vector from a cartesian point.
     *
     * @param p the cartesian point
     */
	public SphericalVector(Point3D.Double p) {
		r = Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
		theta = Math.acos(p.z / r);
		phi = MathUtil.normalizeAngle(Math.atan2(p.y, p.x));
	}

	/**
	 * Convert to a cartesian point.
	 *
	 * @return the cartesian point
	 */
	public void toCartesian(Point3D.Double p) {
		p.x = r * Math.sin(theta) * Math.cos(phi);
		p.y = r * Math.sin(theta) * Math.sin(phi);
		p.z = r * Math.cos(theta);
	}

	/**
	 * Convert to a cartesian point.
	 *
	 * @return the cartesian point
	 */
	public Point3D.Double toCartesian() {
		Point3D.Double p = new Point3D.Double();
		toCartesian(p);
		return p;
	}

}
