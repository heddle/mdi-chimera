package edu.cnu.mdi.chimera.draw;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Point2D;
import java.util.List;
import java.util.Objects;

import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.Crossing;
import edu.cnu.mdi.chimera.patch.BasePatch;
import edu.cnu.mdi.chimera.patch.PrePatch;
import edu.cnu.mdi.chimera.patch.ThetaPatch;
import edu.cnu.mdi.chimera.util.SphericalVector;
import edu.cnu.mdi.chimera.util.ThetaPhi;
import edu.cnu.mdi.graphics.style.LineStyle;
import edu.cnu.mdi.graphics.style.SymbolType;
import edu.cnu.mdi.mapping.MapView2D;
import edu.cnu.mdi.mapping.container.MapContainer;
import edu.cnu.mdi.mapping.graphics.MapGraphics;
import edu.cnu.mdi.ui.colors.X11Colors;

public class DrawPatch {
	
	// Constant for converting theta to latitude in radians (latitude = π/2 - theta)
	private static final double PIOVER2 = Math.PI / 2.0;
		
	/**
	 * Draws a patch on the map using the provided Graphics2D object and MapContainer.
	 *
	 * @param g2         The Graphics2D object to draw on.
	 * @param container  The MapContainer that provides the mapping context.
	 * @param patch      The BasePatch to be drawn.
	 * @param fillColor  The color to fill the patch with.
	 * @param lineColor  The color to draw the patch outline with.
	 * @param lineWidth  The width of the patch outline.
	 * @param lineStyle  The style of the patch outline (e.g., solid, dashed).
	 */
	public static void drawPatch(Graphics2D g2, MapContainer container, BasePatch patch, 
			Color fillColor, Color lineColor, float lineWidth, LineStyle lineStyle) {
		Objects.requireNonNull(g2, "Graphics2D object cannot be null");
		Objects.requireNonNull(container, "MapContainer cannot be null");
		Objects.requireNonNull(patch, "BasePatch cannot be null");
		
		List<ThetaPhi> vertices = patch.getSphericalVertices();
		int num = (vertices != null) ? vertices.size() : 0;
		if (num < 3) {
			// Not enough vertices to draw a patch
			return;
		}
		
		Point2D.Double[] latLonPoints = new Point2D.Double[num];
		for (int i = 0; i < num; i++) {
			ThetaPhi tp = vertices.get(i);
			double lat = PIOVER2 - tp.getTheta(); // Convert theta to latitude
			double lon = tp.getPhi(); // Convert to longitude
			
			if (lat < -Math.PI / 2.0 || lat > Math.PI / 2.0 || lon < -Math.PI || lon > Math.PI) {
				// Latitude out of bounds, skip this vertex
				System.err.println("Warning: Vertex " + i + " has out-of-bounds latitude or longitude. Skipping.");
				continue;
			}
			latLonPoints[i] = new Point2D.Double(lon, lat);
		}
		
		MapGraphics.drawMapPolygon(g2, container, latLonPoints, fillColor, lineColor,
				lineWidth, lineStyle);
	}
	
	/**
	 * Draws the theta crossings of a patch on the map using the provided Graphics2D object and MapContainer.
	 *
	 * @param g2        The Graphics2D object to draw on.
	 * @param container The MapContainer that provides the mapping context.
	 * @param patch     The PrePatch whose theta crossings are to be drawn.
	 */
	public static void drawThetaCrossings(Graphics2D g2, MapContainer container, PrePatch patch) {
		Objects.requireNonNull(g2, "Graphics2D object cannot be null");
		Objects.requireNonNull(container, "MapContainer cannot be null");
		Objects.requireNonNull(patch, "BasePatch cannot be null");
		
		MapView2D view = (MapView2D)(container.getView());
		
		List<Crossing> crossings = patch.getAllThetaCrossings();
		if (crossings == null || crossings.isEmpty()) {
			// No crossings to draw
			return;
		}
		
		for (Crossing crossing : crossings) {
			BaseCurve curve = crossing.curve();
			SphericalVector pos = curve.getSphericalVector(crossing.t());
			double lat = PIOVER2 - pos.theta; // Convert theta to latitude
			double lon = pos.phi; // Convert to longitude
			
			view.drawSymbol(g2, lat, lon, SymbolType.CIRCLE, 8, Color.black, Color.yellow);
		}
	}
	
	public static void drawPhiCrossings(Graphics2D g2, MapContainer container, ThetaPatch patch) {
		Objects.requireNonNull(g2, "Graphics2D object cannot be null");
		Objects.requireNonNull(container, "MapContainer cannot be null");
		Objects.requireNonNull(patch, "BasePatch cannot be null");
		
		MapView2D view = (MapView2D)(container.getView());
		
		List<Crossing> crossings = patch.getAllPhiCrossings();
		if (crossings == null || crossings.isEmpty()) {
			// No crossings to draw
			return;
		}
		
		for (Crossing crossing : crossings) {
			BaseCurve curve = crossing.curve();
			SphericalVector pos = curve.getSphericalVector(crossing.t());
			double lat = PIOVER2 - pos.theta; // Convert theta to latitude
			double lon = pos.phi; // Convert to longitude
			
			view.drawSymbol(g2, lat, lon, SymbolType.SQUARE, 8, Color.black, Color.cyan);
		}
	}

}
