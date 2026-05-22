package edu.cnu.mdi.chimera.app;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

import edu.cnu.mdi.app.BaseMDIApplication;
import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmController;
import edu.cnu.mdi.chimera.dialog.GridSetupDialog;
import edu.cnu.mdi.chimera.dialog.MonteCarloDialog;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.model.ChimeraModel;
import edu.cnu.mdi.chimera.model.ModelChangedEvent;
import edu.cnu.mdi.chimera.view.ChimeraView2D;
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
			INSTANCE = new ChimeraApp(PropertyUtils.TITLE, "Chimera",
					PropertyUtils.CONSOLELOG, true,
					PropertyUtils.BACKGROUND, X11Colors.getX11Color("dark orange"),
					PropertyUtils.FRACTION, 0.8);
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

	        if (event.getType() == ModelChangedEvent.Type.GRID_SPEC_CHANGED) {
	            Log.getInstance().config(chimeraModel.getGridSpec().summary());
	        }
	    });

	}

	// Add menu items for app-specific commands, such as grid setup. This is called from the constructor after the base
	private void modifyMenus() {
		addGridMenu();
		addMonteCarloMenu();
		addAlgorithmMenu();
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

	// Adds the "Monte Carlo" menu with options to generate and clear Monte Carlo points.
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
//
//	    JMenuItem exportJsonItem = new JMenuItem("Export Final Patches JSON...");
//	    exportJsonItem.addActionListener(e -> exportFinalPatchesJson());
//	    algorithmMenu.add(exportJsonItem);
//
//	    algorithmMenu.addSeparator();
//
	    JMenuItem clearAlgorithmItem = new JMenuItem("Clear Algorithm Result");
	    clearAlgorithmItem.addActionListener(e -> algorithmController.clearAlgorithmResult());
	    algorithmMenu.add(clearAlgorithmItem);
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


	/**
	 * Application entry point.
	 *
	 * @param args command-line arguments
	 */
	public static void main(String[] args) {
		BaseMDIApplication.launch(ChimeraApp::getInstance);
	}


}