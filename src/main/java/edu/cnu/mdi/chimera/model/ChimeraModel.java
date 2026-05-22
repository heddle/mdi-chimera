package edu.cnu.mdi.chimera.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

import javax.swing.SwingUtilities;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmResult;
import edu.cnu.mdi.chimera.mc.MonteCarloPoint;

/**
 * Shared application model.
 * <p>
 * The model stores the active grid specification and broadcasts model
 * change events to any interested views or controllers. The model should be the
 * shared center of the application: dialogs update it, algorithm runners
 * update it, and views listen to it.
 * </p>
 */
public final class ChimeraModel {

	/** Active grid specification. */
	private ChimeraGridSpec gridSpec;

	/** Current selected object, such as a cell, prepatch, patch, or patch id. */
	private Object selectedObject;

	/** Current status message. */
	private String statusMessage = "";

	/** Monte Carlo sample points on the current spherical grid. */
	private final List<MonteCarloPoint> monteCarloPoints = new ArrayList<>();

	/** Registered model listeners. */
	private final List<ModelChangedListener> listeners = new ArrayList<>();

	// Current result of the algorithm,
    private ChimeraAlgorithmResult algorithmResult = ChimeraAlgorithmResult.empty();

	 /**
	 * Gets the current algorithm result, or {@code null} if no result is available.
	 *
	 * @return the current algorithm result
	 */
	/**
	 * Creates a model using the default paper test grid.
	 */
	public ChimeraModel() {
		gridSpec = ChimeraGridPresets.paperTestGrid();
		statusMessage = "Loaded default grid: " + gridSpec.getName();
	}

	/**
	 * Gets the active grid specification.
	 *
	 * @return the active grid specification
	 */
	public ChimeraGridSpec getGridSpec() {
		return gridSpec;
	}

	/**
	 * Sets the active grid specification.
	 * <p>
	 * This should be called by the grid setup dialog after the user accepts a new
	 * grid. Changing the grid specification invalidates any computed
	 * algorithm results, so future result fields should be cleared here.
	 * </p>
	 *
	 * @param newGridSpec the new grid specification
	 */
	public void setGridSpec(ChimeraGridSpec newGridSpec) {
		Objects.requireNonNull(newGridSpec, "newGridSpec must not be null");

		ChimeraGridSpec oldSpec = gridSpec;
		if (oldSpec == newGridSpec) {
			return;
		}

		// Clear any computed results that are now invalid. This includes the
		// algorithm result and any Monte Carlo points, which are associated with
		// the spherical grid.

		int oldMonteCarloCount = monteCarloPoints.size();
		monteCarloPoints.clear();

		ChimeraAlgorithmResult oldAlgorithmResult = algorithmResult;
		algorithmResult = ChimeraAlgorithmResult.empty();

		
		gridSpec = newGridSpec;
		selectedObject = null;
		statusMessage = "Grid changed: " + newGridSpec.getName();

		fireModelChanged(ModelChangedEvent.Type.GRID_SPEC_CHANGED, oldSpec, newGridSpec, statusMessage);

		fireModelChanged(ModelChangedEvent.Type.SELECTION_CHANGED, null, null,
				"Selection cleared because the grid changed.");

		if (oldMonteCarloCount > 0) {
			fireModelChanged(ModelChangedEvent.Type.MONTE_CARLO_CLEARED, oldMonteCarloCount, 0,
					"Monte Carlo points cleared because the grid changed.");
		}
		
//		if (oldAlgorithmResult != null && oldAlgorithmResult.getIntersectingCellCount() > 0) {
//			fireModelChanged(ModelChangedEvent.Type.ALGORITHM_RESULT_CHANGED, oldAlgorithmResult, algorithmResult,
//					"Algorithm result cleared because the grid changed.");
//
//			fireModelChanged(ModelChangedEvent.Type.INTERSECTING_CELLS_CHANGED,
//					oldAlgorithmResult.getIntersectingCells(), algorithmResult.getIntersectingCells(),
//					"Intersecting cells cleared because the grid changed.");
//		}


	}

	/**
	 * Resets the model to the default paper test grid.
	 */
	public void resetToDefaultGrid() {
		ChimeraGridSpec oldSpec = gridSpec;
		Object oldSelection = selectedObject;

		gridSpec = ChimeraGridPresets.paperTestGrid();
		selectedObject = null;
		statusMessage = "Model reset to default paper test grid.";

		fireModelChanged(ModelChangedEvent.Type.MODEL_RESET, oldSpec, gridSpec, statusMessage);

		if (oldSelection != null) {
			fireModelChanged(ModelChangedEvent.Type.SELECTION_CHANGED, oldSelection, null, "Selection cleared.");
		}
	}

	/**
	 * Gets the currently selected object.
	 *
	 * @return the selected object, or {@code null}
	 */
	public Object getSelectedObject() {
		return selectedObject;
	}

	/**
	 * Sets the currently selected object.
	 * <p>
	 * The selected object can initially be informal. Later this can become a
	 * stronger type, such as a selected cell, selected prepatch, selected patch, or
	 * selected patch id.
	 * </p>
	 *
	 * @param newSelection the new selected object, or {@code null}
	 */
	public void setSelectedObject(Object newSelection) {
		Object oldSelection = selectedObject;

		if (oldSelection == newSelection) {
			return;
		}

		selectedObject = newSelection;

		fireModelChanged(ModelChangedEvent.Type.SELECTION_CHANGED, oldSelection, newSelection, "Selection changed.");
	}

	/**
	 * Clears the current selection.
	 */
	public void clearSelection() {
		setSelectedObject(null);
	}

	/**
	 * Gets the current status message.
	 *
	 * @return the status message
	 */
	public String getStatusMessage() {
		return statusMessage;
	}

	/**
	 * Sets the current status message.
	 *
	 * @param message the new status message
	 */
	public void setStatusMessage(String message) {
		String oldMessage = statusMessage;
		statusMessage = (message == null) ? "" : message;

		if (!Objects.equals(oldMessage, statusMessage)) {
			fireModelChanged(ModelChangedEvent.Type.STATUS_CHANGED, oldMessage, statusMessage, statusMessage);
		}
	}

	/**
	 * Signals that an algorithm run has started.
	 */
	public void algorithmStarted() {
		setStatusMessage("Algorithm started.");
		fireModelChanged(ModelChangedEvent.Type.ALGORITHM_STARTED, null, null, statusMessage);
	}

	/**
	 * Signals that an algorithm run completed.
	 *
	 * @param message completion message
	 */
	public void algorithmCompleted(String message) {
		setStatusMessage(message == null ? "Algorithm completed." : message);
		fireModelChanged(ModelChangedEvent.Type.ALGORITHM_COMPLETED, null, null, statusMessage);
	}
	
	/**
	 * Gets the most recent algorithm result.
	 *
	 * @return the algorithm result
	 */
	public ChimeraAlgorithmResult getAlgorithmResult() {
		return algorithmResult;
	}
	
	/**
	 * Sets the current algorithm result.
	 * <p>
	 * This is the normal endpoint of an algorithm run. The model owns the result
	 * and notifies views that exact algorithm data has changed.
	 * </p>
	 *
	 * @param result the new algorithm result
	 */
	public void setAlgorithmResult(ChimeraAlgorithmResult result) {
		if (result == null) {
			result = ChimeraAlgorithmResult.empty();
		}

		ChimeraAlgorithmResult oldResult = algorithmResult;
		algorithmResult = result;

//		String message = String.format("Algorithm result updated: %,d intersecting cells.",
//				algorithmResult.getIntersectingCellCount());
//
//		setStatusMessage(message);
//
//		fireModelChanged(ModelChangedEvent.Type.ALGORITHM_RESULT_CHANGED, oldResult, algorithmResult, message);
//
//		fireModelChanged(ModelChangedEvent.Type.INTERSECTING_CELLS_CHANGED,
//				oldResult == null ? null : oldResult.getIntersectingCells(), algorithmResult.getIntersectingCells(),
//				message);
//		
//		fireModelChanged(ModelChangedEvent.Type.PREPATCHES_CHANGED,
//		        oldResult == null ? null : oldResult.getPrepatches(),
//		        algorithmResult.getPrepatches(),
//		        String.format("Built %,d ordinary prepatches.",
//		                algorithmResult.getPrepatchCount()));
//		
//		fireModelChanged(ModelChangedEvent.Type.THETA_PATCHES_CHANGED,
//		        oldResult == null ? null : oldResult.getThetaPatches(),
//		        algorithmResult.getThetaPatches(),
//		        String.format("Built %,d theta patches.",
//		                algorithmResult.getThetaPatchCount()));
//		
//		fireModelChanged(ModelChangedEvent.Type.PATCHES_CHANGED,
//		        oldResult == null ? null : oldResult.getPhiPatches(),
//		        algorithmResult.getPhiPatches(),
//		        String.format("Built %,d final phi patches.",
//		                algorithmResult.getPhiPatchCount()));
	}
	
	/**
	 * Clears the current algorithm result.
	 */
	public void clearAlgorithmResult() {
		ChimeraAlgorithmResult oldResult = algorithmResult;
		algorithmResult = ChimeraAlgorithmResult.empty();
	}


	/**
	 * Signals that an algorithm run failed.
	 *
	 * @param throwable the exception or error that caused the failure
	 */
	public void algorithmFailed(Throwable throwable) {
		String message = (throwable == null) ? "Algorithm failed."
				: "Algorithm failed: " + throwable.getMessage();

		setStatusMessage(message);
		fireModelChanged(ModelChangedEvent.Type.ALGORITHM_FAILED, null, throwable, message);
	}

	/**
	 * Signals that intersecting cells changed.
	 * <p>
	 * This is a placeholder hook for the first algorithm step. 
	 * </p>
	 *
	 * @param newValue the new intersecting-cell data, or {@code null}
	 */
	public void intersectingCellsChanged(Object newValue) {
		fireModelChanged(ModelChangedEvent.Type.INTERSECTING_CELLS_CHANGED, null, newValue,
				"Intersecting cells changed.");
	}

	/**
	 * Signals that prepatches changed.
	 *
	 * @param newValue the new prepatch data, or {@code null}
	 */
	public void prepatchesChanged(Object newValue) {
		fireModelChanged(ModelChangedEvent.Type.PREPATCHES_CHANGED, null, newValue, "Prepatches changed.");
	}

	/**
	 * Signals that theta patches changed.
	 *
	 * @param newValue the new theta-patch data, or {@code null}
	 */
	public void thetaPatchesChanged(Object newValue) {
		fireModelChanged(ModelChangedEvent.Type.THETA_PATCHES_CHANGED, null, newValue, "Theta patches changed.");
	}

	/**
	 * Signals that final patches changed.
	 *
	 * @param newValue the new patch data, or {@code null}
	 */
	public void patchesChanged(Object newValue) {
		fireModelChanged(ModelChangedEvent.Type.PATCHES_CHANGED, null, newValue, "Final patches changed.");
	}

	/**
	 * Signals that display options changed.
	 * <p>
	 * Use this when a view option changes but the grid and computed geometry do
	 * not.
	 * </p>
	 *
	 * @param message optional message
	 */
	public void displayOptionsChanged(String message) {
		fireModelChanged(ModelChangedEvent.Type.DISPLAY_OPTIONS_CHANGED, null, null, message);
	}

	/**
	 * Adds a model changed listener.
	 *
	 * @param listener the listener to add
	 */
	public void addModelChangedListener(ModelChangedListener listener) {
		if (listener == null) {
			return;
		}

		synchronized (listeners) {
			if (!listeners.contains(listener)) {
				listeners.add(listener);
			}
		}
	}

	/**
	 * Removes a model changed listener.
	 *
	 * @param listener the listener to remove
	 */
	public void removeModelChangedListener(ModelChangedListener listener) {
		if (listener == null) {
			return;
		}

		synchronized (listeners) {
			listeners.remove(listener);
		}
	}

	/**
	 * Gets an immutable view of the current Monte Carlo point list.
	 *
	 * @return immutable view of Monte Carlo points
	 */
	public List<MonteCarloPoint> getMonteCarloPoints() {
		return Collections.unmodifiableList(monteCarloPoints);
	}

	/**
	 * Gets the current number of Monte Carlo points.
	 *
	 * @return the number of Monte Carlo points
	 */
	public int getMonteCarloPointCount() {
		return monteCarloPoints.size();
	}

	/**
	 * Adds Monte Carlo points to the current point cloud.
	 *
	 * @param points points to add
	 */
	public void addMonteCarloPoints(List<MonteCarloPoint> points) {
		if (points == null || points.isEmpty()) {
			return;
		}

		monteCarloPoints.addAll(points);

		String message = String.format("Added %,d Monte Carlo points. Total: %,d.", points.size(),
				monteCarloPoints.size());

		setStatusMessage(message);
		fireModelChanged(ModelChangedEvent.Type.MONTE_CARLO_CHANGED, null, points, message);
	}

	/**
	 * Replaces the current Monte Carlo point cloud.
	 *
	 * @param points new points
	 */
	public void setMonteCarloPoints(List<MonteCarloPoint> points) {
		int oldCount = monteCarloPoints.size();

		monteCarloPoints.clear();
		if (points != null) {
			monteCarloPoints.addAll(points);
		}

		String message = String.format("Set Monte Carlo point cloud. Old: %,d. New: %,d.", oldCount,
				monteCarloPoints.size());

		setStatusMessage(message);
		fireModelChanged(ModelChangedEvent.Type.MONTE_CARLO_CHANGED, oldCount, monteCarloPoints.size(), message);
	}

	/**
	 * Clears all Monte Carlo points.
	 */
	public void clearMonteCarloPoints() {
		int oldCount = monteCarloPoints.size();
		if (oldCount == 0) {
			return;
		}

		monteCarloPoints.clear();

		String message = String.format("Cleared %,d Monte Carlo points.", oldCount);
		setStatusMessage(message);
		fireModelChanged(ModelChangedEvent.Type.MONTE_CARLO_CLEARED, oldCount, 0, message);
	}

	/**
	 * Fires a model changed event.
	 * <p>
	 * Events are delivered on the Swing event dispatch thread. If this method is
	 * called from a worker thread, delivery is posted back to the EDT. This keeps
	 * MDI/Swing views safe when the algorithm eventually runs in the
	 * background.
	 * </p>
	 *
	 * @param type     the event type
	 * @param oldValue optional old value
	 * @param newValue optional new value
	 * @param message  optional message
	 */
	public void fireModelChanged(ModelChangedEvent.Type type, Object oldValue, Object newValue, String message) {

		ModelChangedEvent event = new ModelChangedEvent(this, type, oldValue, newValue, message);

		if (SwingUtilities.isEventDispatchThread()) {
			notifyListeners(event);
		} else {
			SwingUtilities.invokeLater(() -> notifyListeners(event));
		}
	}

	/**
	 * Notifies a stable snapshot of the listener list.
	 *
	 * @param event the event to deliver
	 */
	private void notifyListeners(ModelChangedEvent event) {
		List<ModelChangedListener> snapshot;

		synchronized (listeners) {
			snapshot = new ArrayList<>(listeners);
		}

		for (ModelChangedListener listener : snapshot) {
			try {
				listener.modelChanged(event);
			} catch (RuntimeException e) {
				e.printStackTrace();
			}
		}
	}
}