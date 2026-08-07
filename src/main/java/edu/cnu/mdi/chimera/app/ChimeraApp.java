package edu.cnu.mdi.chimera.app;

import java.io.IOException;
import java.nio.file.Path;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.SwingWorker;

import edu.cnu.mdi.app.BaseMDIApplication;
import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmController;
import edu.cnu.mdi.chimera.dialog.GridSetupDialog;
import edu.cnu.mdi.chimera.dialog.MonteCarloDialog;
import edu.cnu.mdi.chimera.export.ChimeraPatchJsonExporter;
import edu.cnu.mdi.chimera.export.ChimeraPatchJsonReader;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Document;
import edu.cnu.mdi.chimera.export.ChimeraPatchJsonReader.ImportedRun;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.model.ChimeraGridPresets;
import edu.cnu.mdi.chimera.model.ChimeraGridPresets.GalleryPreset;
import edu.cnu.mdi.chimera.model.ChimeraModel;
import edu.cnu.mdi.chimera.model.ModelChangedEvent;
import edu.cnu.mdi.chimera.view.ChimeraView2D;
import edu.cnu.mdi.dialog.FileDialogs;
import edu.cnu.mdi.dialog.FileType;
import edu.cnu.mdi.log.Log;
import edu.cnu.mdi.ui.colors.X11Colors;
import edu.cnu.mdi.ui.menu.MenuManager;
import edu.cnu.mdi.util.PropertyUtils;
import edu.cnu.mdi.view.JsonView;
import edu.cnu.mdi.view.LogView;
import edu.cnu.mdi.view.VirtualView;

/**
 * Main application class for Chimera.
 * <p>
 * Chimera computes and visualizes the exact intersection of a rectangular
 * Cartesian grid with an enclosed spherical grid. The application uses MDI for
 * the desktop/view framework.
 * </p>
 */
@SuppressWarnings("serial")
public class ChimeraApp extends BaseMDIApplication {

	/** Singleton instance of the application. */
	private static ChimeraApp INSTANCE;

	/** Main 2D map view. */
	private ChimeraView2D chimeraView2D;

	/** the log view */
	private LogView logView;

	/** JSON view for debugging. Not always visible. */
	private JsonView jsonView;

	private ChimeraAlgorithmController algorithmController;

	/** File type used by the final-patch interchange exporter. */
	private static final FileType PATCH_JSON = FileType.of(
			"Chimera final patches (*.json)", "json");

	/** the main data model for Chimera with the grid information */
	private final static ChimeraModel chimeraModel;
	static {
		chimeraModel = new ChimeraModel();
	}

	/**
	 * Creates the application.
	 */
	private ChimeraApp(Object... keyVals) {
		super(keyVals);
		algorithmController = new ChimeraAlgorithmController(chimeraModel);
		modifyMenus();
	}

	/**
	 * Returns the singleton instance of the application.
	 *
	 * @return the application instance
	 */
	public static ChimeraApp getInstance() {
		if (INSTANCE == null) {
			INSTANCE = new ChimeraApp(PropertyUtils.TITLE, "Chimera", PropertyUtils.CONSOLELOG, true,
					PropertyUtils.BACKGROUND, X11Colors.getX11Color("dark orange"), PropertyUtils.FRACTION, 0.8);
		}
		return INSTANCE;
	}

	@Override
	protected int getVirtualDesktopColumns() {
		return 3;
	}

	@Override
	protected void addInitialViews() {
		chimeraView2D = new ChimeraView2D(chimeraModel);

		// Log view is useful but not always visible.
		logView = new LogView();
		logView.setVisible(false);

		// JSON view is useful for debugging but not always visible.
		jsonView = new JsonView();
		jsonView.setVisible(false);

		chimeraModel.addModelChangedListener(event -> {
			String message = event.getMessage();
			if (message != null && !message.isBlank()) {
				System.out.println("[Chimera] " + message);
			}
			
			switch (event.getType()) {
			case ALGORITHM_STARTED:
				Log.getInstance().info("Algorithm started.");
				break;
			case ALGORITHM_COMPLETED:
				Log.getInstance().info("Algorithm completed.");
				chimeraView2D.refresh();
				break;
			case ALGORITHM_FAILED:
				Log.getInstance().error("Algorithm failed. ");
				break;
			case GRID_SPEC_CHANGED:
				Log.getInstance().config("Grid specification changed. ");
				Log.getInstance().config(chimeraModel.getGridSpec().summary());
				break;

			default:
				// For other event types, no specific action is needed here.
				break;
			}

		});

	}

	// Add menu items for app-specific commands, such as grid setup. This is called
	// from the constructor after the base
	private void modifyMenus() {
		addGridMenu();
		addTestGridGalleryMenu();
		addMonteCarloMenu();
		addAlgorithmMenu();
	}

	/**
	 * Adds a gallery of progressively larger deterministic test grids. Selecting
	 * an entry loads the grid but does not automatically run the algorithm.
	 */
	private void addTestGridGalleryMenu() {
		JMenu galleryMenu = new JMenu("Test Grid Gallery");
		getJMenuBar().add(galleryMenu);

		for (GalleryPreset preset : GalleryPreset.values()) {
			JMenuItem item = new JMenuItem(preset.toString());
			item.addActionListener(event -> selectGalleryGrid(preset));
			galleryMenu.add(item);
		}
	}

	/** Loads one gallery grid after protecting any current computed result. */
	private void selectGalleryGrid(GalleryPreset preset) {
		var currentResult = chimeraModel.getAlgorithmResult();
		if (currentResult != null && currentResult.getPatchCount() > 0) {
			int answer = JOptionPane.showConfirmDialog(this,
					"Loading another grid will clear the current algorithm result.\nContinue?",
					"Load Test Grid", JOptionPane.OK_CANCEL_OPTION,
					JOptionPane.QUESTION_MESSAGE);
			if (answer != JOptionPane.OK_OPTION) {
				return;
			}
		}

		ChimeraGridSpec spec = preset.createGridSpec();
		chimeraModel.setGridSpec(spec);
		Log.getInstance().config("Loaded test-grid gallery entry: " + spec.getName());
	}

	/// Adds the "Grid" menu with options to set up the grid.
	private void addGridMenu() {
		JMenuItem gridItem = new JMenuItem("Grid Specifications...");
		gridItem.addActionListener(e -> gridDialog());

		JMenu fileMenu = MenuManager.getInstance().getFileMenu();
		int quitIndex = fileMenu.getItemCount() - 1;

		// Insert the new item just above the separator before Quit.
		fileMenu.insert(gridItem, quitIndex);

	}

	// Adds the "Monte Carlo" menu with options to generate and clear Monte Carlo
	// points.
	private void addMonteCarloMenu() {
		JMenu mcMenu = new JMenu("Monte Carlo");
		getJMenuBar().add(mcMenu);

		JMenuItem generateMonteCarloItem = new JMenuItem("Generate Monte Carlo...");
		generateMonteCarloItem.addActionListener(e -> MonteCarloDialog.showDialog(this, chimeraModel));
		mcMenu.add(generateMonteCarloItem);

		JMenuItem clearMonteCarloItem = new JMenuItem("Clear Monte Carlo");
		clearMonteCarloItem.addActionListener(e -> chimeraModel.clearMonteCarloPoints());
		mcMenu.add(clearMonteCarloItem);
	}

	/**
	 * Adds the algorithm menu.
	 */
	private void addAlgorithmMenu() {
		JMenu algorithmMenu = new JMenu("Algorithm");
		getJMenuBar().add(algorithmMenu);

		JMenuItem runAlgorithmItem = new JMenuItem("Run Algorithm");
		runAlgorithmItem.addActionListener(e -> algorithmController.runAlgorithm());
		algorithmMenu.add(runAlgorithmItem);

		JMenuItem exportJsonItem = new JMenuItem("Export Final Patches JSON...");
		exportJsonItem.addActionListener(e -> exportFinalPatchesJson());
		algorithmMenu.add(exportJsonItem);

		JMenuItem validateJsonItem = new JMenuItem("Validate Final Patches JSON...");
		validateJsonItem.addActionListener(e -> validateFinalPatchesJson());
		algorithmMenu.add(validateJsonItem);

		JMenuItem importJsonItem = new JMenuItem("Import Final Patches JSON...");
		importJsonItem.addActionListener(e -> importFinalPatchesJson());
		algorithmMenu.add(importJsonItem);

		algorithmMenu.addSeparator();

		JMenuItem clearAlgorithmItem = new JMenuItem("Clear Algorithm Result");
		clearAlgorithmItem.addActionListener(e -> algorithmController.clearAlgorithmResult());
		algorithmMenu.add(clearAlgorithmItem);
	}

	/** Selects and validates an interchange file without changing the model. */
	private void validateFinalPatchesJson() {
		var input = FileDialogs.openFile(this, "chimera-final-patches",
				"Validate Final Patches", PATCH_JSON);
		if (input.isEmpty()) {
			return;
		}
		Path path = input.get();
		chimeraModel.setStatusMessage("Validating final patches in " + path + "...");
		new SwingWorker<Document, Void>() {
			@Override
			protected Document doInBackground() throws IOException {
				return ChimeraPatchJsonReader.validate(path);
			}

			@Override
			protected void done() {
				try {
					Document document = get();
					String message = "Valid Chimera patch file: "
							+ document.patches().size() + " final patches.";
					chimeraModel.setStatusMessage(message);
					JOptionPane.showMessageDialog(ChimeraApp.this, message,
							"Validation Successful", JOptionPane.INFORMATION_MESSAGE);
				} catch (Exception exception) {
					showJsonFailure("Validation Failed", exception);
				}
			}
		}.execute();
	}

	/** Selects, validates, reconstructs, and displays an exported patch run. */
	private void importFinalPatchesJson() {
		var input = FileDialogs.openFile(this, "chimera-final-patches",
				"Import Final Patches", PATCH_JSON);
		if (input.isEmpty()) {
			return;
		}
		Path path = input.get();
		chimeraModel.setStatusMessage("Importing final patches from " + path + "...");
		new SwingWorker<ImportedRun, Void>() {
			@Override
			protected ImportedRun doInBackground() throws IOException {
				return ChimeraPatchJsonReader.read(path);
			}

			@Override
			protected void done() {
				try {
					ImportedRun imported = get();
					chimeraModel.setGridSpec(imported.gridSpec());
					chimeraModel.setAlgorithmResult(imported.algorithmResult());
					chimeraModel.algorithmCompleted("Imported "
							+ imported.algorithmResult().getPatchCount()
							+ " final patches from " + path + ".");
				} catch (Exception exception) {
					showJsonFailure("Import Failed", exception);
				}
			}
		}.execute();
	}

	/** Reports the root cause from a background JSON task. */
	private void showJsonFailure(String title, Exception exception) {
		Throwable cause = exception.getCause() == null ? exception : exception.getCause();
		String message = cause.getMessage() == null ? cause.toString() : cause.getMessage();
		chimeraModel.setStatusMessage(message);
		Log.getInstance().error(message);
		JOptionPane.showMessageDialog(this, message, title, JOptionPane.ERROR_MESSAGE);
	}

	/**
	 * Prompts for a destination and exports the current final patches without
	 * blocking the Swing event-dispatch thread during JSON generation.
	 */
	private void exportFinalPatchesJson() {
		var result = chimeraModel.getAlgorithmResult();
		if (result == null || result.getPatches() == null || result.getPatches().isEmpty()) {
			JOptionPane.showMessageDialog(this,
					"Run the algorithm before exporting final patches.",
					"No Final Patches", JOptionPane.INFORMATION_MESSAGE);
			return;
		}

		String baseName = chimeraModel.getGridSpec().getName()
				.replaceAll("[^A-Za-z0-9._-]+", "-")
				.replaceAll("^-+|-+$", "");
		if (baseName.isBlank()) {
			baseName = "chimera";
		}

		var destination = FileDialogs.saveFile(this, "chimera-final-patches",
				"Export Final Patches", baseName + "-patches.json", PATCH_JSON);
		if (destination.isEmpty()) {
			return;
		}

		Path output = destination.get();
		ChimeraGridSpec spec = chimeraModel.getGridSpec();
		chimeraModel.setStatusMessage("Exporting final patches to " + output + "...");

		new SwingWorker<Void, Void>() {
			@Override
			protected Void doInBackground() throws IOException {
				ChimeraPatchJsonExporter.write(output, spec, result);
				return null;
			}

			@Override
			protected void done() {
				try {
					get();
					String message = "Exported " + result.getPatchCount()
							+ " final patches to " + output + ".";
					chimeraModel.setStatusMessage(message);
					Log.getInstance().info(message);
				} catch (Exception exception) {
					Throwable cause = exception.getCause() == null
							? exception : exception.getCause();
					String message = "Could not export final patches: " + cause.getMessage();
					chimeraModel.setStatusMessage(message);
					Log.getInstance().error(message);
					JOptionPane.showMessageDialog(ChimeraApp.this, message,
							"Export Failed", JOptionPane.ERROR_MESSAGE);
				}
			}
		}.execute();
	}

	// Show the grid setup dialog and update the model if the user accepts a new
	private void gridDialog() {
		ChimeraGridSpec spec = GridSetupDialog.showDialog(null, chimeraModel.getGridSpec());
		if (spec != null) {
			chimeraModel.setGridSpec(spec);
		}
	}

	/**
	 * Place the views in the virtual desktop in a reasonable default layout.
	 *
	 * <p>
	 * Note: this placement will be ignored if the user has a persisted
	 * layout/config.
	 * </p>
	 */
	@Override
	protected void defaultViewLayout() {
		virtualViewMove(chimeraView2D, 0, VirtualView.CENTER);
		virtualViewMove(logView, 2, VirtualView.UPPERLEFT);
		virtualViewMove(jsonView, 2, VirtualView.BOTTOMRIGHT);
	}

	// -------- convenience methods to access views and model --------

	/**
	 * Get the model
	 * 
	 * @return the ChimeraModel instance
	 */
	public ChimeraModel getChimeraModel() {
		return chimeraModel;
	}

	/**
	 * Get the GridSpec from the model
	 * 
	 * @return the GridSpec from the model
	 */
	public ChimeraGridSpec getGridSpec() {
		return chimeraModel.getGridSpec();
	}

	/**
	 * Get the CartesianGrid from the model's GridSpec
	 * 
	 * @return the CartesianGrid from the model's GridSpec
	 */
	public CartesianGrid getCartesianGrid() {
		return getGridSpec().getCartesianGrid();
	}

	/**
	 * Get the SphericalGrid from the model's GridSpec
	 * 
	 * @return the SphericalGrid from the model's GridSpec
	 */
	public SphericalGrid getSphericalGrid() {
		return getGridSpec().getSphericalGrid();
	}
	
	/**
	 * Get the radius of the spherical grid from the model's GridSpec
	 * 
	 * @return the radius of the spherical grid from the model's GridSpec
	 */
	public double getRadius() {
		return getSphericalGrid().getRadius();
	}

	/**
	 * Application entry point.
	 *
	 * @param args command-line arguments
	 */
	public static void main(String[] args) {
		BaseMDIApplication.launch(ChimeraApp::getInstance);
	}

}
