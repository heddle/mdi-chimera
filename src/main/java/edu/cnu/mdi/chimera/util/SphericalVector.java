package edu.cnu.mdi.chimera.util;

public class SphericalVector {

	/** the polar angle in radians */
	public double theta;

	/** the azimuthal angle in radians */
	public double phi;

	/** the radius */
	public double r;

	/**
	 * Create a spherical vector. Azimuthal angles will be normalized to the range
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
		return String.format("θ: %.3f φ: %.3f", theta, phi);
	}
	
	public String toStringDegrees(String name) {
		double theta = Math.toDegrees(this.theta);
		double phi = Math.toDegrees(this.phi);
		return String.format("%s θ: %.3f φ %.3f", name, theta, phi);
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
	
	/**
	 * Calculate the great-circle distance between this point and another point on the sphere.
	 *
	 * @param other the other spherical vector
	 * @return the great-circle distance in the same units as {@code r}
	 */
	public double distanceTo(SphericalVector other) {
		return gcDistance(this.theta, this.phi, other.theta, other.phi) * r;
	}

	/**
	 * Calculates the normalized (R = 1) great-circle distance between two points on
	 * a sphere using spherical coordinates.
	 * 
	 * @param theta1 Polar angle (colatitude) of point 1 in radians [0, pi]
	 * @param phi1   Azimuthal angle (longitude) of point 1 in radians [0, 2*pi] or [-pi, pi]
	 * @param theta2 Polar angle (colatitude) of point 2 in radians [0, pi]
	 * @param phi2   Azimuthal angle (longitude) of point 2 in radians [0, 2*pi] or [-pi, pi]
	 * @return The angular distance between the two points in radians [0, pi]
	 */
	public static double gcDistance(double theta1, double phi1, double theta2, double phi2) {
		// Delta phi (difference in longitude)
		double dPhi = phi2 - phi1;

		// Haversine-based calculation adapted for polar angles (theta)
		// sin²(Δθ/2) + sin(theta1) * sin(theta2) * sin²(Δɸ/2)
		double sinDTheta = Math.sin((theta2 - theta1) / 2.0);
		double sinDPhi = Math.sin(dPhi / 2.0);

		double a = (sinDTheta * sinDTheta) + Math.sin(theta1) * Math.sin(theta2) * (sinDPhi * sinDPhi);

		// Clamp 'a' to handle eventual floating-point inaccuracies yielding values
		// slightly > 1.0
		// which would cause Math.sqrt(a) or Math.asin to return NaN.
		a = Math.min(1.0, Math.max(a, 0.0));

		// Normalized distance (angular separation in radians)
		return 2.0 * Math.asin(Math.sqrt(a));
	}

}
