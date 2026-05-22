package edu.cnu.mdi.chimera.dialog;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Window;
import java.text.ParseException;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.grid.ThetaSpacing;
import edu.cnu.mdi.chimera.model.ChimeraGridPresets;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.model.CoordinateSystem;
import edu.cnu.mdi.chimera.model.LengthUnit;
import edu.cnu.mdi.dialog.SimpleDialog;

/**
 * Dialog for creating or replacing the active grid specification.
 * <p>
 * The first version intentionally supports generated grids and a few presets
 * rather than direct editing of arbitrary coordinate arrays. The underlying
 * model still supports fully nonuniform grids.
 * </p>
 */
@SuppressWarnings("serial")
public class GridSetupDialog extends SimpleDialog {

	/** High-precision display format for floating-point grid parameters. */
	private static final String DOUBLE_SPINNER_FORMAT = "0.###############";

	/** Small step for generated-grid floating point controls. */
	private static final double GRID_DOUBLE_STEP = 1.0e-5;

	/** Radius spinner step. */
	private static final double RADIUS_STEP = 1.0e-5;

	/** Grid specification supplied when the dialog was opened. */
	private ChimeraGridSpec initialGridSpec;

	/** True while controls are being populated programmatically. */
	private boolean populatingControls;

    /** OK command. */
    private static final String OK = " OK ";

    /** Cancel command. */
    private static final String CANCEL = " Cancel ";

    /** Preset choices. */
    private enum Preset {
        CURRENT_GRID("Current grid"),
        GENERATED("Generated grid"),
        PAPER_TEST("Paper Test Grid"),
        SMALL_DEBUG("Small Debug Grid"),
        COARSE_GSM("Coarse GSM Grid");

        private final String label;

        Preset(String label) {
            this.label = label;
        }

        @Override
        public String toString() {
            return label;
        }
    }

    /** Preset selector. */
    private JComboBox<Preset> presetCombo;

    /** Theta spacing selector. */
    private JComboBox<ThetaSpacing> thetaSpacingCombo;

    /** X min spinner. */
    private JSpinner xminSpinner;

    /** X max spinner. */
    private JSpinner xmaxSpinner;

    /** X cell count spinner. */
    private JSpinner nxSpinner;

    /** Y min spinner. */
    private JSpinner yminSpinner;

    /** Y max spinner. */
    private JSpinner ymaxSpinner;

    /** Y cell count spinner. */
    private JSpinner nySpinner;

    /** Z min spinner. */
    private JSpinner zminSpinner;

    /** Z max spinner. */
    private JSpinner zmaxSpinner;

    /** Z cell count spinner. */
    private JSpinner nzSpinner;

    /** Radius spinner. */
    private JSpinner radiusSpinner;

    /** Theta cell count spinner. */
    private JSpinner nthetaSpinner;

    /** Phi cell count spinner. */
    private JSpinner nphiSpinner;

    /** Summary area. */
    private JTextArea summaryArea;

    /** The selected grid specification, if the dialog closed with OK. */
    private ChimeraGridSpec selectedGridSpec;

    /** Whether the dialog was cancelled. */
    private boolean cancelled = true;

    /**
     * Creates the grid setup dialog.
     *
     * @param initialSpec the initial grid specification
     */
    public GridSetupDialog(ChimeraGridSpec initialSpec) {
        super("Grid Setup", true, OK, CANCEL);

        initialGridSpec = (initialSpec != null)
                ? initialSpec
                : ChimeraGridPresets.paperTestGrid();

        selectedGridSpec = initialGridSpec;

        populateFromGridSpec(initialGridSpec);
        updateSummary();
    }

    /**
     * Shows a modal grid setup dialog.
     *
     * @param parent parent component
     * @param initialSpec initial grid specification
     * @return the selected grid specification, or {@code null} if cancelled
     */
    public static ChimeraGridSpec showDialog(Component parent, ChimeraGridSpec initialSpec) {
        GridSetupDialog dialog = new GridSetupDialog(initialSpec);

        if (parent != null) {
            Window window = SwingUtilities.getWindowAncestor(parent);
            if (window != null) {
                dialog.setLocationRelativeTo(window);
            }
        }

        dialog.setVisible(true);
        return dialog.isCancelled() ? null : dialog.getSelectedGridSpec();
    }

    /**
     * Populates the editable controls from a grid specification.
     * <p>
     * If the grid is nonuniform, the dialog shows the grid's min/max and cell
     * counts. Pressing OK while the "Current grid" preset is selected preserves the
     * exact nonuniform grid; pressing OK while "Generated grid" is selected creates
     * a new uniform grid from the displayed min/max/count values.
     * </p>
     *
     * @param spec the grid specification to display
     */
    /**
     * Populates the editable controls from a grid specification.
     * <p>
     * If the grid is nonuniform, the dialog shows the grid's min/max and cell
     * counts. Pressing OK while the "Current grid" preset is selected preserves the
     * exact nonuniform grid; pressing OK after editing a field switches the dialog
     * to "Generated grid" and creates a new generated grid from the displayed
     * values.
     * </p>
     *
     * @param spec the grid specification to display
     */
    private void populateFromGridSpec(ChimeraGridSpec spec) {
        if (spec == null) {
            return;
        }

        populatingControls = true;

        try {
            CartesianGrid cart = spec.getCartesianGrid();
            SphericalGrid sphere = spec.getSphericalGrid();

            setSpinnerValue(xminSpinner, cart.getXMin());
            setSpinnerValue(xmaxSpinner, cart.getXMax());
            setSpinnerValue(nxSpinner, cart.getNumXCells());

            setSpinnerValue(yminSpinner, cart.getYMin());
            setSpinnerValue(ymaxSpinner, cart.getYMax());
            setSpinnerValue(nySpinner, cart.getNumYCells());

            setSpinnerValue(zminSpinner, cart.getZMin());
            setSpinnerValue(zmaxSpinner, cart.getZMax());
            setSpinnerValue(nzSpinner, cart.getNumZCells());

            setSpinnerValue(radiusSpinner, sphere.getRadius());
            setSpinnerValue(nthetaSpinner, sphere.getNumThetaCells());
            setSpinnerValue(nphiSpinner, sphere.getNumPhiCells());

            ThetaSpacing inferred = inferThetaSpacing(sphere);
            thetaSpacingCombo.setSelectedItem(inferred);
        } finally {
            populatingControls = false;
        }
    }

    /**
     * Handles a user edit to one of the generated-grid controls.
     * <p>
     * Any manual edit means the dialog is no longer returning the exact current
     * grid or one of the named presets. It is now building a generated grid from
     * the visible control values.
     * </p>
     */
    private void controlChanged() {
        if (populatingControls) {
            return;
        }

        if (presetCombo != null && presetCombo.getSelectedItem() != Preset.GENERATED) {
            presetCombo.setSelectedItem(Preset.GENERATED);
        } else {
            updateSummary();
        }
    }


    /**
     * Infers the closest theta-spacing mode for display.
     *
     * @param sphere the spherical grid
     * @return the inferred theta spacing
     */
    private static ThetaSpacing inferThetaSpacing(SphericalGrid sphere) {
        double[] theta = sphere.getThetaGrid().getPoints();

        if (theta.length < 3) {
            return ThetaSpacing.UNIFORM_THETA;
        }

        double thetaErr = maxSpacingError(theta);

        double[] mu = new double[theta.length];
        for (int i = 0; i < theta.length; i++) {
            mu[i] = Math.cos(theta[i]);
        }

        double muErr = maxSpacingError(mu);

        return (muErr < thetaErr) ? ThetaSpacing.UNIFORM_COS_THETA : ThetaSpacing.UNIFORM_THETA;
    }

    /**
     * Measures the maximum deviation from uniform spacing.
     *
     * @param values values to inspect
     * @return maximum spacing deviation
     */
    private static double maxSpacingError(double[] values) {
        if (values.length < 3) {
            return 0.0;
        }

        double expected = (values[values.length - 1] - values[0]) / (values.length - 1);
        double maxErr = 0.0;

        for (int i = 1; i < values.length; i++) {
            double actual = values[i] - values[i - 1];
            maxErr = Math.max(maxErr, Math.abs(actual - expected));
        }

        return maxErr;
    }

    /**
     * Sets a spinner value.
     *
     * @param spinner the spinner
     * @param value the new value
     */
    private static void setSpinnerValue(JSpinner spinner, double value) {
        spinner.setValue(value);
    }

    /**
     * Sets a spinner value.
     *
     * @param spinner the spinner
     * @param value the new value
     */
    private static void setSpinnerValue(JSpinner spinner, int value) {
        spinner.setValue(value);
    }

    @Override
    protected Component createCenterComponent() {
        JPanel main = new JPanel(new BorderLayout(6, 6));
        main.setBorder(BorderFactory.createEmptyBorder(8, 10, 4, 10));

        main.add(createCoordinatePanel(), BorderLayout.NORTH);
        main.add(createGridPanel(), BorderLayout.CENTER);
        main.add(createSummaryPanel(), BorderLayout.SOUTH);

        return main;
    }

    /**
     * Creates the coordinate-system panel.
     *
     * @return the coordinate panel
     */
    private Component createCoordinatePanel() {
        JTextArea area = new JTextArea();
        area.setEditable(false);
        area.setOpaque(false);
        area.setText(
                "Coordinate system: " + CoordinateSystem.GSM + "\n"
                        + "Origin: Earth center\n"
                        + "Length unit: " + LengthUnit.EARTH_RADII + "\n"
                        + "Axes: +X toward Sun, +Z tilted toward north magnetic pole, +Y completes GSM system");

        JPanel panel = new JPanel(new BorderLayout());
        panel.setBorder(BorderFactory.createTitledBorder("Physical Coordinates"));
        panel.add(area, BorderLayout.CENTER);
        return panel;
    }

    /**
     * Creates the editable grid panel.
     *
     * @return the grid panel
     */
    private Component createGridPanel() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));

        JPanel presetPanel = new JPanel(new GridLayout(1, 2, 8, 8));
        presetPanel.setBorder(BorderFactory.createTitledBorder("Preset"));

        presetCombo = new JComboBox<>(Preset.values());
        presetCombo.setSelectedItem(Preset.CURRENT_GRID);

        presetCombo.addActionListener(e -> {
            if (populatingControls) {
                return;
            }

            Preset preset = (Preset) presetCombo.getSelectedItem();

            if (preset == Preset.PAPER_TEST) {
                populateFromGridSpec(ChimeraGridPresets.paperTestGrid());
            } else if (preset == Preset.SMALL_DEBUG) {
                populateFromGridSpec(ChimeraGridPresets.smallDebugGrid());
            } else if (preset == Preset.COARSE_GSM) {
                populateFromGridSpec(ChimeraGridPresets.coarseGsmGrid());
            } else if (preset == Preset.CURRENT_GRID) {
                populateFromGridSpec(initialGridSpec);
            } else if (preset == Preset.GENERATED) {
                // Leave the current visible control values alone.
                // They now define the generated grid.
            }

            updateSummary();
        });

        presetPanel.add(new JLabel("Grid source"));
        presetPanel.add(presetCombo);

        JPanel cartPanel = createCartesianPanel();
        JPanel spherePanel = createSphericalPanel();

        JPanel editPanel = new JPanel(new GridLayout(1, 2, 8, 8));
        editPanel.add(cartPanel);
        editPanel.add(spherePanel);

        panel.add(presetPanel, BorderLayout.NORTH);
        panel.add(editPanel, BorderLayout.CENTER);

        return panel;
    }

    /**
     * Creates the Cartesian grid editing panel.
     *
     * @return the Cartesian grid panel
     */
    private JPanel createCartesianPanel() {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(BorderFactory.createTitledBorder("Cartesian GSM Grid"));

        xminSpinner = doubleSpinner(-6.0, -1.0e6, 1.0e6, GRID_DOUBLE_STEP);
        xmaxSpinner = doubleSpinner(6.0, -1.0e6, 1.0e6, GRID_DOUBLE_STEP);
        nxSpinner = intSpinner(36, 1, 10000);

        yminSpinner = doubleSpinner(-6.0, -1.0e6, 1.0e6, GRID_DOUBLE_STEP);
        ymaxSpinner = doubleSpinner(6.0, -1.0e6, 1.0e6, GRID_DOUBLE_STEP);
        nySpinner = intSpinner(36, 1, 10000);

        zminSpinner = doubleSpinner(-6.0, -1.0e6, 1.0e6, GRID_DOUBLE_STEP);
        zmaxSpinner = doubleSpinner(6.0, -1.0e6, 1.0e6, GRID_DOUBLE_STEP);
        nzSpinner = intSpinner(36, 1, 10000);

        addAxisRows(panel, "X", xminSpinner, xmaxSpinner, nxSpinner, 0);
        addAxisRows(panel, "Y", yminSpinner, ymaxSpinner, nySpinner, 3);
        addAxisRows(panel, "Z", zminSpinner, zmaxSpinner, nzSpinner, 6);

        return panel;
    }

    /**
     * Adds the three label/spinner rows for one Cartesian axis using GridBagLayout.
     *
     * @param panel  destination panel
     * @param axis   axis letter ("X", "Y", or "Z")
     * @param min    min spinner
     * @param max    max spinner
     * @param cells  cell-count spinner
     * @param startRow first GridBag row to use
     */
    private static void addAxisRows(JPanel panel, String axis,
                                    JSpinner min, JSpinner max, JSpinner cells,
                                    int startRow) {
        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(2, 4, 2, 6);
        lc.fill = GridBagConstraints.NONE;
        lc.weightx = 0;

        GridBagConstraints sc = new GridBagConstraints();
        sc.fill = GridBagConstraints.HORIZONTAL;
        sc.insets = new Insets(2, 0, 2, 4);
        sc.weightx = 1;

        String[][] rows = {
            { axis + " min",   null },
            { axis + " max",   null },
            { axis + " cells", null }
        };
        JSpinner[] spinners = { min, max, cells };

        for (int i = 0; i < 3; i++) {
            lc.gridx = 0; lc.gridy = startRow + i;
            sc.gridx = 1; sc.gridy = startRow + i;

            panel.add(new JLabel(rows[i][0]), lc);
            panel.add(spinners[i], sc);
        }
    }

    /**
     * Creates the spherical grid editing panel.
     *
     * @return the spherical grid panel
     */
    private JPanel createSphericalPanel() {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(BorderFactory.createTitledBorder("Spherical Grid"));

        radiusSpinner = doubleSpinner(4.60993, 1.0e-12, 1.0e6, RADIUS_STEP);
        nthetaSpinner = intSpinner(48, 1, 10000);
        nphiSpinner = intSpinner(32, 1, 10000);
        thetaSpacingCombo = new JComboBox<>(ThetaSpacing.values());

        thetaSpacingCombo.addActionListener(e -> {
            Component c = thetaSpacingCombo.getTopLevelAncestor();
            if (c instanceof GridSetupDialog dialog) {
                dialog.controlChanged();
            }
        });

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(2, 4, 2, 6);
        lc.fill = GridBagConstraints.NONE;
        lc.weightx = 0;

        GridBagConstraints fc = new GridBagConstraints();
        fc.fill = GridBagConstraints.HORIZONTAL;
        fc.insets = new Insets(2, 0, 2, 4);
        fc.weightx = 1;

        String[] labels = { "Radius, R\u2091", "Theta cells", "Phi cells", "Theta spacing" };
        Component[] fields = { radiusSpinner, nthetaSpinner, nphiSpinner, thetaSpacingCombo };

        for (int i = 0; i < labels.length; i++) {
            lc.gridx = 0; lc.gridy = i;
            fc.gridx = 1; fc.gridy = i;
            panel.add(new JLabel(labels[i]), lc);
            panel.add(fields[i], fc);
        }

        // Push rows to the top so they don't spread when the panel is taller than Cartesian
        GridBagConstraints filler = new GridBagConstraints();
        filler.gridx = 0; filler.gridy = labels.length;
        filler.weighty = 1;
        panel.add(Box.createGlue(), filler);

        return panel;
    }

    /**
     * Creates the summary panel.
     *
     * @return the summary panel
     */
    private Component createSummaryPanel() {
        summaryArea = new JTextArea(4, 70);
        summaryArea.setEditable(false);
        summaryArea.setLineWrap(true);
        summaryArea.setWrapStyleWord(false);
        summaryArea.setBorder(BorderFactory.createEmptyBorder(3, 4, 3, 4));
        // Slightly smaller monospaced font keeps long coordinate lines readable
        Font base = summaryArea.getFont();
        summaryArea.setFont(base.deriveFont(base.getSize2D() - 1f));

        JScrollPane scroll = new JScrollPane(summaryArea,
                ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scroll.setBorder(BorderFactory.createEmptyBorder());

        JPanel panel = new JPanel(new BorderLayout());
        panel.setBorder(BorderFactory.createTitledBorder("Grid Summary"));
        panel.add(scroll, BorderLayout.CENTER);
        return panel;
    }

    /**
     * Gets the selected grid specification.
     *
     * @return the selected grid specification
     */
    public ChimeraGridSpec getSelectedGridSpec() {
        return selectedGridSpec;
    }

    /**
     * Checks whether the dialog was cancelled.
     *
     * @return {@code true} if cancelled
     */
    public boolean isCancelled() {
        return cancelled;
    }

    @Override
    protected void handleCommand(String command) {
        if (command != null && command.trim().equals("OK")) {
            try {
            	commitSpinnerEdits();
            	selectedGridSpec = buildGridSpec();
            	cancelled = false;
                setVisible(false);
            } catch (RuntimeException e) {
                summaryArea.setText("Grid setup error:\n" + e.getMessage());
            }
        } else {
            cancelled = true;
            setVisible(false);
        }
    }

    /**
     * Builds the grid specification from the current dialog state.
     * <p>
     * Important: the Paper Test Grid is not exactly equivalent to a generated
     * uniform Cartesian grid with the same min/max/count values. Therefore, when
     * the Cartesian controls still match the initial grid's displayed Cartesian
     * extent and cell counts, this method preserves the exact initial Cartesian
     * grid arrays. This allows the user to change only the spherical angular grid
     * without silently changing the Cartesian cell geometry.
     * </p>
     *
     * @return the selected grid specification
     */
    private ChimeraGridSpec buildGridSpec() {
        Preset preset = (Preset) presetCombo.getSelectedItem();

        if (preset == Preset.CURRENT_GRID) {
            return initialGridSpec;
        }

        if (preset == Preset.PAPER_TEST) {
            return ChimeraGridPresets.paperTestGrid();
        }

        if (preset == Preset.SMALL_DEBUG) {
            return ChimeraGridPresets.smallDebugGrid();
        }

        if (preset == Preset.COARSE_GSM) {
            return ChimeraGridPresets.coarseGsmGrid();
        }

        double xmin = doubleValue(xminSpinner);
        double xmax = doubleValue(xmaxSpinner);
        int nx = intValue(nxSpinner);

        double ymin = doubleValue(yminSpinner);
        double ymax = doubleValue(ymaxSpinner);
        int ny = intValue(nySpinner);

        double zmin = doubleValue(zminSpinner);
        double zmax = doubleValue(zmaxSpinner);
        int nz = intValue(nzSpinner);

        double radius = doubleValue(radiusSpinner);
        int ntheta = intValue(nthetaSpinner);
        int nphi = intValue(nphiSpinner);
        ThetaSpacing spacing = (ThetaSpacing) thetaSpacingCombo.getSelectedItem();

        CartesianGrid cartesianGrid;

        if (cartesianControlsMatchInitialGrid()) {
            /*
             * Preserve the exact existing Cartesian grid. This is essential for the
             * Paper Test Grid, whose coordinate arrays are not identical to a newly
             * generated uniform min/max/count grid.
             */
            cartesianGrid = initialGridSpec.getCartesianGrid();
        } else {
            cartesianGrid = CartesianGrid.uniform(
                    xmin, xmax, nx,
                    ymin, ymax, ny,
                    zmin, zmax, nz);
        }


        SphericalGrid sphericalGrid = SphericalGrid.generated(radius, ntheta, nphi, spacing);

        String name;

        if (cartesianControlsMatchInitialGrid()) {
            name = initialGridSpec.getName() + " with generated spherical grid";
        } else {
            name = "Generated GSM Grid";
        }

        return new ChimeraGridSpec(name, CoordinateSystem.GSM,
                LengthUnit.EARTH_RADII, cartesianGrid, sphericalGrid);
    }

    /**
     * Checks whether the Cartesian controls still describe the initial grid's
     * displayed Cartesian extent and cell counts.
     * <p>
     * Matching these values does not prove that a generated uniform grid would be
     * identical to the initial grid. In fact, for the Paper Test Grid it is not.
     * Therefore this method is used to decide when to preserve the exact initial
     * Cartesian grid arrays instead of rebuilding them from min/max/count.
     * </p>
     *
     * @return {@code true} if the Cartesian controls match the initial grid
     */
    private boolean cartesianControlsMatchInitialGrid() {
        if (initialGridSpec == null) {
            return false;
        }

        CartesianGrid initialCartesian = initialGridSpec.getCartesianGrid();

        return nearlyEqual(doubleValue(xminSpinner), initialCartesian.getXMin())
                && nearlyEqual(doubleValue(xmaxSpinner), initialCartesian.getXMax())
                && intValue(nxSpinner) == initialCartesian.getNumXCells()

                && nearlyEqual(doubleValue(yminSpinner), initialCartesian.getYMin())
                && nearlyEqual(doubleValue(ymaxSpinner), initialCartesian.getYMax())
                && intValue(nySpinner) == initialCartesian.getNumYCells()

                && nearlyEqual(doubleValue(zminSpinner), initialCartesian.getZMin())
                && nearlyEqual(doubleValue(zmaxSpinner), initialCartesian.getZMax())
                && intValue(nzSpinner) == initialCartesian.getNumZCells();
    }

    /**
     * Floating-point comparison for dialog round trips.
     *
     * @param a first value
     * @param b second value
     * @return {@code true} if the values are close enough to be considered equal
     */
    private static boolean nearlyEqual(double a, double b) {
        double scale = Math.max(1.0, Math.max(Math.abs(a), Math.abs(b)));
        return Math.abs(a - b) <= 1.0e-12 * scale;
    }

    /**
     * Updates the summary display.
     */
    private void updateSummary() {
        if (summaryArea == null) {
            return;
        }

        try {
            summaryArea.setText(buildGridSpec().summary());
        } catch (RuntimeException e) {
            summaryArea.setText("Grid setup error:\n" + e.getMessage());
        }
    }

    /**
     * Creates a high-precision floating-point spinner.
     * <p>
     * The default Swing number editor may display values such as {@code 4.60993}
     * as {@code 4.61}. If the dialog then rebuilds a generated grid from the
     * visible controls, that display rounding can silently change the geometry.
     * This method installs a high-precision editor so that grid coordinates and
     * the sphere radius survive a round trip through the dialog.
     * </p>
     *
     * @param value initial value
     * @param min minimum value
     * @param max maximum value
     * @param step step size
     * @return the spinner
     */
    private static JSpinner doubleSpinner(double value, double min, double max, double step) {
        JSpinner spinner = new JSpinner(new SpinnerNumberModel(value, min, max, step));

        JSpinner.NumberEditor editor =
                new JSpinner.NumberEditor(spinner, DOUBLE_SPINNER_FORMAT);

        editor.getTextField().setColumns(12);
        spinner.setEditor(editor);

        // Pin height so the spinner arrow buttons never grow oversized when the
        // parent panel stretches vertically.
        Dimension pref = spinner.getPreferredSize();
        spinner.setMaximumSize(new Dimension(Integer.MAX_VALUE, pref.height));

        spinner.addChangeListener(e -> {
            Component c = spinner.getTopLevelAncestor();
            if (c instanceof GridSetupDialog dialog) {
                dialog.controlChanged();
            }
        });

        return spinner;
    }
    /**
     * Creates an integer spinner.
     *
     * @param value initial value
     * @param min minimum value
     * @param max maximum value
     * @return the spinner
     */
    private static JSpinner intSpinner(int value, int min, int max) {
        JSpinner spinner = new JSpinner(new SpinnerNumberModel(value, min, max, 1));

        Dimension pref = spinner.getPreferredSize();
        spinner.setMaximumSize(new Dimension(Integer.MAX_VALUE, pref.height));

        spinner.addChangeListener(e -> {
            Component c = spinner.getTopLevelAncestor();
            if (c instanceof GridSetupDialog dialog) {
                dialog.controlChanged();
            }
        });
        return spinner;
    }

    /**
     * Adds a label/component pair to a grid panel.
     *
     * @param panel destination panel
     * @param label label text
     * @param component editor component
     */
    private static void addLabeled(JPanel panel, String label, Component component) {
        panel.add(new JLabel(label));
        panel.add(component);
    }

    /**
     * Gets a double value from a spinner.
     *
     * @param spinner the spinner
     * @return the double value
     */
    private static double doubleValue(JSpinner spinner) {
        return ((Number) spinner.getValue()).doubleValue();
    }

    /**
     * Gets an integer value from a spinner.
     *
     * @param spinner the spinner
     * @return the integer value
     */
    private static int intValue(JSpinner spinner) {
        return ((Number) spinner.getValue()).intValue();
    }

    /**
     * Commits any text currently being edited in the spinner text fields.
     * <p>
     * Without this, a user can type a new value and press OK while the spinner
     * model still contains the previous value.
     * </p>
     */
    private void commitSpinnerEdits() {
        commitSpinnerEdit(xminSpinner);
        commitSpinnerEdit(xmaxSpinner);
        commitSpinnerEdit(nxSpinner);

        commitSpinnerEdit(yminSpinner);
        commitSpinnerEdit(ymaxSpinner);
        commitSpinnerEdit(nySpinner);

        commitSpinnerEdit(zminSpinner);
        commitSpinnerEdit(zmaxSpinner);
        commitSpinnerEdit(nzSpinner);

        commitSpinnerEdit(radiusSpinner);
        commitSpinnerEdit(nthetaSpinner);
        commitSpinnerEdit(nphiSpinner);
    }

    /**
     * Commits one spinner edit.
     *
     * @param spinner spinner
     */
    private static void commitSpinnerEdit(JSpinner spinner) {
        if (spinner == null) {
            return;
        }

        try {
            spinner.commitEdit();
        } catch (ParseException e) {
            throw new IllegalArgumentException("Invalid numeric value: " + spinner.getValue(), e);
        }
    }


}