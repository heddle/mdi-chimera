package edu.cnu.mdi.chimera.map;

import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Shape;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

import edu.cnu.mdi.container.IContainer;
import edu.cnu.mdi.mapping.graphics.MapGraphics;
import edu.cnu.mdi.mapping.projection.EProjection;
import edu.cnu.mdi.mapping.projection.IMapProjection;
import edu.cnu.mdi.mapping.theme.MapTheme;

/**
 * Application-supplied cylindrical equal-area projection.
 *
 * <p>
 * This projection is the Lambert cylindrical equal-area projection, sometimes
 * historically associated with Archimedes' theorem that the area of a zone on
 * a sphere equals the area of the corresponding zone on a circumscribed
 * cylinder.
 * </p>
 *
 * <p>
 * The projection equations are:
 * </p>
 *
 * <pre>
 * x = wrap(lambda - lambda0)
 * y = sin(phi)
 * </pre>
 *
 * <p>
 * where {@code lambda} is longitude, {@code phi} is latitude, and
 * {@code lambda0} is the central longitude. The inverse equations are:
 * </p>
 *
 * <pre>
 * lambda = wrap(x + lambda0)
 * phi    = asin(y)
 * </pre>
 *
 * <p>
 * Coordinates follow the MDI map convention: geographic points use
 * {@code x = longitude}, {@code y = latitude}, in radians.
 * </p>
 */
public class ArchimedesLambertCylindricalProjection implements IMapProjection {

    /** Projection x minimum. */
    private static final double XMIN = -Math.PI;

    /** Projection x maximum. */
    private static final double XMAX = Math.PI;

    /** Projection y minimum. */
    private static final double YMIN = -1.0;

    /** Projection y maximum. */
    private static final double YMAX = 1.0;

    /** Central longitude in radians. */
    private double centralLongitude;

    /** Active map theme. */
    private MapTheme theme;

    /**
     * Creates the projection with central longitude zero.
     *
     * @param theme the map theme
     */
    public ArchimedesLambertCylindricalProjection(MapTheme theme) {
        this(0.0, theme);
    }

    /**
     * Creates the projection.
     *
     * @param centralLongitude central longitude in radians
     * @param theme the map theme
     */
    public ArchimedesLambertCylindricalProjection(double centralLongitude, MapTheme theme) {
        this.centralLongitude = wrapLongitude(centralLongitude);
        setTheme(theme == null ? MapTheme.light() : theme);
    }

    @Override
    public void latLonToXY(Point2D.Double latLon, Point2D.Double xy) {
        double lon = latLon.x;
        double lat = clampLatitude(latLon.y);

        xy.x = wrapLongitude(lon - centralLongitude);
        xy.y = Math.sin(lat);
    }

    @Override
    public void latLonFromXY(Point2D.Double latLon, Point2D.Double xy) {
        if (!isPointOnMap(xy)) {
            latLon.x = Double.NaN;
            latLon.y = Double.NaN;
            return;
        }

        latLon.x = wrapLongitude(xy.x + centralLongitude);
        latLon.y = Math.asin(clamp(xy.y, -1.0, 1.0));
    }

    @Override
    public boolean isPointVisible(Point2D.Double latLon) {
        return latLon != null
                && Double.isFinite(latLon.x)
                && Double.isFinite(latLon.y)
                && latLon.y >= -Math.PI / 2.0
                && latLon.y <= Math.PI / 2.0;
    }

    @Override
    public boolean isPointOnMap(Point2D.Double xy) {
        return xy != null
                && Double.isFinite(xy.x)
                && Double.isFinite(xy.y)
                && xy.x >= XMIN
                && xy.x <= XMAX
                && xy.y >= YMIN
                && xy.y <= YMAX;
    }

    @Override
    public void drawMapOutline(Graphics2D g2, IContainer container) {
        Shape outline = createClipShape(container);
        if (outline == null) {
            return;
        }

        java.awt.Color oldColor = g2.getColor();
        g2.setColor(theme.getOutlineColor());
        g2.draw(outline);
        g2.setColor(oldColor);
    }

    @Override
    public void drawLatitudeLine(Graphics2D g2, IContainer container, double latitude) {
        double lat = clampLatitude(latitude);
		MapGraphics.drawHorizontalLatitudeLine(g2, container, this, lat,
				getCentralLongitude(), theme);
    }

    @Override
	public void drawLongitudeLine(Graphics2D g2, IContainer container, double longitude) {
		double lon = wrapLongitude(longitude);
		MapGraphics.drawVerticalLongitudeLine(g2, container, this, lon, theme);
	}

    @Override
    public Shape createClipShape(IContainer container) {
        Point p0 = new Point();
        Point p1 = new Point();
        Point p2 = new Point();
        Point p3 = new Point();

        container.worldToLocal(p0, new Point2D.Double(XMIN, YMIN));
        container.worldToLocal(p1, new Point2D.Double(XMAX, YMIN));
        container.worldToLocal(p2, new Point2D.Double(XMAX, YMAX));
        container.worldToLocal(p3, new Point2D.Double(XMIN, YMAX));

        Path2D.Double path = new Path2D.Double();
        path.moveTo(p0.x, p0.y);
        path.lineTo(p1.x, p1.y);
        path.lineTo(p2.x, p2.y);
        path.lineTo(p3.x, p3.y);
        path.closePath();

        return path;
    }

    /**
     * This is an application-supplied projection, not an MDI built-in enum
     * projection.
     *
     * @return {@code null}
     */
    @Override
    public EProjection getProjection() {
        return null;
    }

    @Override
    public String name() {
        return "Archimedes / Lambert Cylindrical";
    }

    @Override
    public Rectangle2D.Double getXYBounds() {
        return new Rectangle2D.Double(XMIN, YMIN, XMAX - XMIN, YMAX - YMIN);
    }

    @Override
    public MapTheme getTheme() {
        return theme;
    }

    @Override
    public void setTheme(MapTheme theme) {
        if (theme == null) {
            throw new IllegalArgumentException("theme must not be null.");
        }
        this.theme = theme;
    }

    @Override
    public boolean supportsRecenter() {
        return true;
    }

    @Override
    public boolean recenterOn(Point2D.Double latLon) {
        if (latLon == null || !Double.isFinite(latLon.x)) {
            return false;
        }

        centralLongitude = wrapLongitude(latLon.x);
        return true;
    }

    @Override
    public boolean crossesSeam(double lon1, double lon2) {
        double x1 = wrapLongitude(lon1 - centralLongitude);
        double x2 = wrapLongitude(lon2 - centralLongitude);
        return Math.abs(x2 - x1) > Math.PI;
    }

    /**
     * Gets the central longitude.
     *
     * @return the central longitude in radians
     */
    public double getCentralLongitude() {
        return centralLongitude;
    }

    /**
     * Sets the central longitude.
     *
     * @param centralLongitude the central longitude in radians
     */
    public void setCentralLongitude(double centralLongitude) {
        this.centralLongitude = wrapLongitude(centralLongitude);
    }

    /**
     * Clamps a latitude to the valid range.
     *
     * @param latitude latitude in radians
     * @return clamped latitude
     */
    private static double clampLatitude(double latitude) {
        return clamp(latitude, -Math.PI / 2.0, Math.PI / 2.0);
    }

    /**
     * Clamps a value.
     *
     * @param value input value
     * @param min minimum value
     * @param max maximum value
     * @return clamped value
     */
    private static double clamp(double value, double min, double max) {
        return Math.max(min, Math.min(max, value));
    }
    

}
