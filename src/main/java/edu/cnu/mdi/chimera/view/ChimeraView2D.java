package edu.cnu.mdi.chimera.view;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmResult;
import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.cell.KissGeometry;
import edu.cnu.mdi.chimera.draw.DrawPatch;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.map.ArchimedesLambertCylindricalProjection;
import edu.cnu.mdi.chimera.map.ChimeraMapControlPanel;
import edu.cnu.mdi.chimera.mc.MonteCarloPoint;
import edu.cnu.mdi.chimera.model.ChimeraModel;
import edu.cnu.mdi.chimera.patch.BasePatch;
import edu.cnu.mdi.component.OptionPanel;
import edu.cnu.mdi.container.IContainer;
import edu.cnu.mdi.graphics.SymbolDraw;
import edu.cnu.mdi.graphics.style.LineStyle;
import edu.cnu.mdi.graphics.style.SymbolType;
import edu.cnu.mdi.graphics.toolbar.ToolBits;
import edu.cnu.mdi.log.Log;
import edu.cnu.mdi.mapping.MapView2D;
import edu.cnu.mdi.mapping.container.MapContainer;
import edu.cnu.mdi.mapping.projection.IMapProjection;
import edu.cnu.mdi.ui.colors.ColorMapSelectorPanel;

import edu.cnu.mdi.ui.colors.ScientificColorMap;
import edu.cnu.mdi.util.PropertyUtils;
import edu.cnu.mdi.util.UnicodeUtils;

@SuppressWarnings("serial")
public class ChimeraView2D extends MapView2D implements OptionPanel.OptionPanelListener,
ColorMapSelectorPanel.ColorMapChangeListener {
	// Constant for π/2, used in projection calculations.
	private final double PIOVER2 = Math.PI / 2.0;

	// Threshold for snapping points to the poles when checking for pole 
	// involvement in patch classification. This is a
	// visualization-only tolerance, not a geometric alteration of stored patch geometry.
	private static final double ANGULAR_SNAP_TOL = 1.0e-5;
	

	/** for theta index display in feedback */
	public static final String NTHETA = "n" + UnicodeUtils.SMALL_THETA;

	/** for phi index display in feedback */
	public static final String NPHI = "n" + UnicodeUtils.SMALL_PHI;

	// Reusable point for projection calculations to avoid unnecessary object
	// creation.
	private int[] indexArray = new int[5];


	/** Shared Chimera model. */
	private final ChimeraModel model;


	/** Active color map for Monte Carlo points. */
	private ScientificColorMap monteCarloColorMap = ScientificColorMap.VIRIDIS;

	/** Dot size for Monte Carlo points, in pixels. */
	private int monteCarloDotSize = 2;
		
	// the dsplay options panel
	private final ChimeraOptionPanel optionPanel;

	/**
	 * Constructor: Initializes the Chimera 2D view with the given model.
	 *
	 * @param model the shared Chimera model (must not be null)
	 * @throws NullPointerException if the model is null
	 */
	public ChimeraView2D(ChimeraModel model) {
		super(PropertyUtils.TITLE, "Chimera 2D", 
				PropertyUtils.FRACTION, 0.78, 
				PropertyUtils.ASPECT, 1.1,
				PropertyUtils.TOOLBARBITS, (ToolBits.MAPTOOLS | ToolBits.ZOOMTOOLS) & ~ToolBits.STATUS);
		
		Objects.requireNonNull(model, "model must not be null.");
		this.model = model;
		
		/*
		 * These calls are intentionally made after super(...) returns, avoiding the
		 * superclass-constructor override/lifecycle problem.
		 */
		setMapControlPanel(new ChimeraMapControlPanel(this));
		
		ChimeraSidePanel sidePanel = new ChimeraSidePanel(this, monteCarloColorMap, this);
		addCustomSidePanelComponent(sidePanel);
		optionPanel = sidePanel.optionPanel;
		
		
		setProjection(new ArchimedesLambertCylindricalProjection(getCurrentMapTheme()));
		
		model.addModelChangedListener(event -> {
			System.out.println("Model changed: " + event.getType());

			switch (event.getType()) {
			case GRID_SPEC_CHANGED, MODEL_RESET, DISPLAY_OPTIONS_CHANGED, SELECTION_CHANGED, INTERSECTING_CELLS_CHANGED,
					ALGORITHM_RESULT_CHANGED, ALGORITHM_OPTIONS_CHANGED, PREPATCHES_CHANGED, THETA_PATCHES_CHANGED,
					PATCHES_CHANGED, MONTE_CARLO_CHANGED, MONTE_CARLO_CLEARED ->
				gridChange();
			default -> {
				// No redraw needed.
			}
			}
		});
		pack();
	}
	
	// Redraw when the grid changes or options change
	private void gridChange() {
		Log.getInstance().config("Grid changed ");
		refresh();
	}
	
	/**
	 * Override to draw custom map rendering
	 */
	@Override
	protected void drawCustomMapContent(Graphics2D g, IContainer container) {

		drawMonteCarloPoints(g, container);

		if (optionPanel.showSphericalGrid()) {
			// Draw phi lines (longitudes)
			drawPhiLines(g, container);

			// Draw theta lines (latitudes, sort of)
			drawThetaLines(g, container);
		}
		
		// Draw Kiss markers if enabled
		if (optionPanel.showKissMarkers()) {
			drawKissMarkers(g, container);
		}
		
		// Drawprepatches if enabled
		if (optionPanel.showPrepatches()) {
			ChimeraAlgorithmResult result = model.getAlgorithmResult();
			if (result != null) {
				drawPatchList(g, container, result.getPrePatches(), null, Color.red);
			}
		}
	}
	
	// Draw a list of patches with specified fill and line colors. 
	// This is used for prepatches, theta patches, and phi patches.
	private void drawPatchList(Graphics2D g, IContainer container, List<? extends BasePatch> patches, 
			Color fillColor, Color lineColor) {
		if (patches == null) {
			return;
		}
		MapContainer mapContainer = (MapContainer) container;

		for (BasePatch patch : patches) {
			DrawPatch.drawPatch(g, mapContainer, patch, fillColor, lineColor, 1.5f, LineStyle.SOLID);
		}
		
	}
	
	// Draw markers for cells identified as "Kiss" cells by the algorithm, 
	// if the option is enabled.
	private void drawKissMarkers(Graphics2D g, IContainer container) {
		if (!optionPanel.showKissMarkers()) {
			return;
		}
		
		ChimeraAlgorithmResult result = model.getAlgorithmResult();
		if (result == null || result.getKissCells() == null) {
			return;
		}
		
		List<Cell> kissCells = result.getKissCells();
		for (Cell cell : kissCells) {
			KissGeometry kissGeometry = cell.getKissGeometry();
			double[] latLon = kissGeometry.getCenter();
			double lat = latLon[0];
			double lon = latLon[1];
			drawSymbol(g, lat, lon, SymbolType.CIRCLE, 8, Color.RED, new Color(255, 0, 0, 128));
		}

	}

	/**
	 * Draws Monte Carlo points using their stored GSM angular coordinates and
	 * stored color-map values.
	 *
	 * @param g         graphics context
	 * @param container map container
	 */
	private void drawMonteCarloPoints(Graphics2D g, IContainer container) {
		if (!optionPanel.showMonteCarlo() || model.getMonteCarloPointCount() == 0) {
			return;
		}

		IMapProjection projection = getProjection();

		Point2D.Double latLon = new Point2D.Double();
		Point2D.Double xy = new Point2D.Double();
		Point screen = new Point();

		int s = monteCarloDotSize;
		int half = s / 2;

		for (MonteCarloPoint point : model.getMonteCarloPoints()) {
			latLon.x = point.phi();
			latLon.y = Math.PI / 2.0 - point.theta();

			projection.latLonToXY(latLon, xy);

			if (!Double.isFinite(xy.x) || !Double.isFinite(xy.y) || !projection.isPointOnMap(xy)) {
				continue;
			}

			container.worldToLocal(screen, xy);

			Color color = monteCarloColorMap.colorAt(point.colorValue());
			g.setColor(color);
			g.fillRect(screen.x - half, screen.y - half, s, s);
		}
	}

	// method to draw phi lines (longitudes) on the map.
	private void drawPhiLines(Graphics2D g, IContainer container) {
		SphericalGrid grid = model.getGridSpec().getSphericalGrid();
		Grid1D phiGrid = grid.getPhiGrid();
		IMapProjection projection = getProjection();

		for (double phi : phiGrid.getPoints()) {
			// Convert phi to map coordinates and draw the line
			// This is a placeholder; actual implementation would depend on the projection
			// and map scale
			projection.drawLongitudeLine(g, container, phi);
		}

	}
	
	// method to draw theta lines (latitudes) on the map.
	private void drawThetaLines(Graphics2D g, IContainer container) {
		SphericalGrid grid = model.getGridSpec().getSphericalGrid();
		Grid1D thetaGrid = grid.getThetaGrid();
		IMapProjection projection = getProjection();

		for (double theta : thetaGrid.getPoints()) {
			double latitude = PIOVER2 - theta; // Convert theta to latitude

			// Snap near-equator and near-pole values. The historical test grids use
			// rounded values such as 1.5708 instead of Math.PI / 2, which can create
			// near-antipode artifacts in azimuthal projections.
			if (Math.abs(latitude) < ANGULAR_SNAP_TOL) {
				latitude = 0.0;
			} else if (Math.abs(latitude - PIOVER2) < ANGULAR_SNAP_TOL) {
				latitude = PIOVER2;
			} else if (Math.abs(latitude + PIOVER2) < ANGULAR_SNAP_TOL) {
				latitude = -PIOVER2;
			}

			projection.drawLatitudeLine(g, container, latitude);
		}
	}


	// Override to disable standard graticules (latitude/longitude lines) if
	// desired.
	@Override
	protected boolean useStandardGraticules() {
		return false; // Disable standard graticules to avoid cluttering the map
	}
	/**
	 * Override to set a custom side panel width suitable for controls.
	 *
	 * @return the side panel width in pixels
	 */
	@Override
	protected int getSidePanelWidth() {
		return 300;
	}
	
	@Override
	protected boolean includeShapeFileMenu() {
		return false;
	}



	//callback for display toggles
	@Override
	public void optionStateChanged(OptionPanel source, String label, boolean selected) {
		refresh();
	}


	@Override
	public void colorMapChanged(ColorMapSelectorPanel source, ScientificColorMap map) {
		monteCarloColorMap = map;
		refresh();
	}

	
	/**
	 * Gets feedback strings for the current mouse position, including GSM
	 * coordinates and grid indices.
	 *
	 * @param container       map container
	 * @param pp              mouse position in panel coordinates
	 * @param wp              mouse position in world coordinates (GSM lat/lon)
	 * @param feedbackStrings list to add feedback strings to
	 */
	@Override
	public void getFeedbackStrings(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {

		IMapProjection projection = getProjection();

		if (projection.isPointOnMap(wp)) {
			projection.latLonFromXY(latLon, wp);
			double gsmTheta = 90 - Math.toDegrees(latLon.y);
			double gsmPhi = Math.toDegrees(latLon.x);

			double r = model.getGridSpec().getSphericalGrid().getRadius();
			double sinTheta = Math.sin(Math.toRadians(gsmTheta));
			double x = r * sinTheta * Math.cos(Math.toRadians(gsmPhi));
			double y = r * sinTheta * Math.sin(Math.toRadians(gsmPhi));
			double z = r * Math.cos(Math.toRadians(gsmTheta));

			String polarStr = String.format("(r, %s, %s) = (%.2fRe, %.2f%s, %.2f%s)", UnicodeUtils.SMALL_THETA,
					UnicodeUtils.SMALL_PHI, r, gsmTheta, DEG, gsmPhi, DEG);

			String carStr = String.format("(x, y, z) = (%.2fRe, %.2fRe, %.2fRe)", x, y, z);

			feedbackStrings.add(polarStr);
			feedbackStrings.add(carStr);
			
			model.getGridSpec().getPatchIndices(Math.toRadians(gsmTheta), Math.toRadians(gsmPhi), r, indexArray);
			feedbackStrings.add(String.format("nx=%d, ny=%d" + ", nz=%d, %s=%d, %s=%d", indexArray[0], indexArray[1],
					indexArray[2], NTHETA, indexArray[3], NPHI, indexArray[4]));

			
			addGridFeedback(container, pp, wp, feedbackStrings);
			addMonteCarloFeedback(container, pp, wp, feedbackStrings);
			addPrepatchFeedback(container, pp, wp, feedbackStrings);
			addThetaPatchFeedback(container, pp, wp, feedbackStrings);
			addPhiPatchFeedback(container, pp, wp, feedbackStrings);
			algorithmFeedback(container, pp, wp, feedbackStrings);
		}
	}
	
	// Add feedback about the Cartesian and spherical grids
	private void addGridFeedback(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {
		String colorPrefix = "$coral$";
		CartesianGrid cgrid = model.getGridSpec().getCartesianGrid();
		SphericalGrid sgrid = model.getGridSpec().getSphericalGrid();
		
		ArrayList<String> cgridFeedback = cgrid.feedbackStrings();
		for (String s : cgridFeedback) {
			feedbackStrings.add(colorPrefix + s);
		}
		
		ArrayList<String> sgridFeedback = sgrid.feedbackStrings();
		for (String s : sgridFeedback) {
			feedbackStrings.add(colorPrefix + s);
		}
		
	}
	
	// Add feedback about MonteCarlo points
	private void addMonteCarloFeedback(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {
		String colorPrefix = "$wheat$";
		int mcCount = model.getMonteCarloPointCount();
		feedbackStrings.add(colorPrefix + "Monte Carlo points: " + mcCount);
	}

	private void addPrepatchFeedback(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {
		String colorPrefix = "$magenta$";
	}

	private void addThetaPatchFeedback(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {
		String colorPrefix = "$brown$";
	}
	
	private void addPhiPatchFeedback(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {
		String colorPrefix = "$purple$";
	}
	
	private void algorithmFeedback(IContainer container, Point pp, Point2D.Double wp, List<String> feedbackStrings) {
		ChimeraAlgorithmResult result = model.getAlgorithmResult();
		if (result != null) {
			String colorPrefix = "$orange$";
			result.feedbackSummary(colorPrefix, feedbackStrings);
		}
	}




}
