package edu.cnu.mdi.chimera.view;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmResult;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.model.ChimeraModel;
import edu.cnu.mdi.view.AbstractViewInfo;

/**
 * Provides metadata and help text for the {@link ChimeraView2D} "Info" dialog.
 *
 * <p>This object is constructed by {@link ChimeraView2D#getViewInfo()} and
 * passed to the framework's info-dialog machinery. It replaces the generic,
 * demo-oriented {@code MapViewInfo} (cities/countries on a sample map) that
 * {@code ChimeraView2D} would otherwise inherit from {@link
 * edu.cnu.mdi.mapping.MapView2D MapView2D}, and instead describes what
 * Chimera itself does: intersecting a Cartesian grid with an enclosed
 * spherical shell.</p>
 *
 * <h2>Live technical notes</h2>
 * <p>{@link #getTechnicalNotes()} reads the owning view's current {@link
 * ChimeraModel} so the reported grid size, coordinate system, and patch
 * counts always reflect what is actually loaded, rather than a hardcoded
 * description that could drift out of date.</p>
 */
public class ChimeraView2DInfo extends AbstractViewInfo {

    /** Reference to the owning view, used to obtain live model state. */
    private final ChimeraView2D chimeraView;

    /**
     * Creates a view-info object for the given Chimera 2D view.
     *
     * @param chimeraView the owning view; must not be {@code null}
     */
    public ChimeraView2DInfo(ChimeraView2D chimeraView) {
        this.chimeraView = chimeraView;
    }

    /**
     * {@inheritDoc}
     *
     * @return {@code "Chimera 2D — Grid/Sphere Intersection"}
     */
    @Override
    public String getTitle() {
        return "Chimera 2D — Grid/Sphere Intersection";
    }

    /**
     * {@inheritDoc}
     *
     * @return a short paragraph describing the purpose of this view
     */
    @Override
    public String getPurpose() {
        return "Chimera computes and visualizes the exact geometric intersection of a "
             + "rectangular Cartesian grid with an enclosed spherical shell, expressed in "
             + "GSM (Geocentric Solar Magnetospheric) coordinates and Earth radii. This "
             + "view renders the spherical grid's θ/φ lines together with the "
             + "algorithm's intersecting cells, prepatches, θ-patches, and final "
             + "patches on an Archimedes–Lambert equal-area cylindrical projection of "
             + "the sphere.";
    }

    /**
     * {@inheritDoc}
     *
     * <p>Ordered steps covering the normal workflow: pick or build a grid, run
     * the algorithm, optionally validate with Monte Carlo sampling, then read
     * or export the results.</p>
     */
    @Override
    public List<String> getUsageSteps() {
        return List.of(
                "Pick a preset from the Test Grid Gallery menu, or build a custom grid "
              + "with Grid → Grid Specifications…, which defines the Cartesian "
              + "grid, the spherical shell, and its GSM/Earth-radii interpretation.",

                "Run Algorithm → Run Algorithm to intersect the Cartesian grid with "
              + "the sphere. The algorithm proceeds in stages — intersecting cells, "
              + "then prepatches, then θ-patches, then final patches — and each "
              + "stage becomes its own visible layer.",

                "Use Monte Carlo → Generate Monte Carlo… to scatter up to 50 "
              + "million sample points on the spherical shell, colored by the selected "
              + "scientific color map, as an independent visual check on patch coverage "
              + "and classification.",

                "Toggle layers in the side panel's Display Options: Spherical Grid, "
              + "Prepatches, θ Patches, Patches, Kiss Markers, and Monte Carlo.",

                "Hover over the map for live feedback: GSM (θ, φ) in degrees, "
              + "Cartesian (x, y, z) in Earth radii, the Cartesian and spherical grid cell "
              + "indices under the cursor, and the nearest patch's estimated area.",

                "Export, validate, or import final-patch results as JSON from the "
              + "Algorithm menu, e.g. to compare runs or hand results off to other tools."
        );
    }

    /**
     * {@inheritDoc}
     *
     * <p>Supplementary notes about the algorithm's intermediate concepts and
     * how the view reacts to model changes.</p>
     */
    @Override
    public List<String> getUsageBullets() {
        return List.of(
                "\"Kiss\" cells are cells that intersect the sphere only tangentially "
              + "— a special case the algorithm tracks separately from ordinary "
              + "intersecting cells. Enable Kiss Markers to see them highlighted.",

                "Hovering over a rendered patch highlights it and overlays its θ- or "
              + "φ-crossing construction lines, which is useful when debugging the "
              + "patch-building algorithm.",

                "Loading a new grid (gallery preset or custom) clears any existing "
              + "algorithm result and Monte Carlo point cloud, since both are tied to the "
              + "previous spherical grid.",

                "The Log and JSON views are available but hidden by default — enable "
              + "them from the View menu when diagnosing algorithm or model behavior."
        );
    }

    /**
     * {@inheritDoc}
     *
     * <p>Returns the standard map-container navigation bindings shared by all
     * {@code MapView2D}-based views.</p>
     */
    @Override
    public Map<String, String> getKeyboardShortcuts() {
        Map<String, String> shortcuts = new LinkedHashMap<>();
        shortcuts.put("Scroll Wheel",          "Zoom in / out");
        shortcuts.put("Right Click + Drag",    "Pan the map");
        return shortcuts;
    }

    /**
     * {@inheritDoc}
     *
     * <p>Reports the currently loaded grid's name, coordinate system, size,
     * and (if computed) the current algorithm and Monte Carlo state, read
     * live from the view's {@link ChimeraModel} so the numbers never go
     * stale.</p>
     *
     * @return technical notes describing the live grid and result state
     */
    @Override
    public String getTechnicalNotes() {
        ChimeraModel model = chimeraView.getModel();
        ChimeraGridSpec spec = model.getGridSpec();
        CartesianGrid cartesian = spec.getCartesianGrid();
        SphericalGrid spherical = spec.getSphericalGrid();
        ChimeraAlgorithmResult result = model.getAlgorithmResult();

        String grid = String.format(
                "Grid: \"%s\" [%s, %s]. Cartesian: %d×%d×%d cells. "
              + "Spherical: %d×%d cells, radius %.4f.",
                spec.getName(), spec.getCoordinateSystem(), spec.getLengthUnit(),
                cartesian.getXGrid().numCells(), cartesian.getYGrid().numCells(),
                cartesian.getZGrid().numCells(),
                spherical.getNumThetaCells(), spherical.getNumPhiCells(),
                spherical.getRadius());

        String algorithm = (result == null || result.getIntersectingCells() == null)
                ? "Algorithm not yet run."
                : String.format(
                        "Intersecting cells: %d (kiss: %d). Prepatches: %d. "
                      + "θ-patches: %d. Final patches: %d.",
                        result.getIntersectingCellCount(), result.getKissCellCount(),
                        result.getPrePatchCount(), result.getThetaPatchCount(),
                        result.getPatchCount());

        String monteCarlo = String.format("Monte Carlo points: %,d.",
                model.getMonteCarloPointCount());

        return grid + " " + algorithm + " " + monteCarlo;
    }

    /**
     * {@inheritDoc}
     *
     * <p>Uses a deep space blue accent, distinct from the framework's default
     * orange and appropriate for a magnetospheric grid tool.</p>
     */
    @Override
    protected String getAccentColorHex() {
        return "#1b4f72";
    }

    /**
     * {@inheritDoc}
     *
     * @return {@code "Chimera — MDI Framework"}
     */
    @Override
    public String getFooter() {
        return "Chimera — MDI Framework";
    }
}
