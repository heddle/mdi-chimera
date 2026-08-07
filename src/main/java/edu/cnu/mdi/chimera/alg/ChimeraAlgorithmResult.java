package edu.cnu.mdi.chimera.alg;

import java.util.Collections;
import java.util.List;

import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.cell.IntersectionType;
import edu.cnu.mdi.chimera.patch.Patch;
import edu.cnu.mdi.chimera.patch.PrePatch;
import edu.cnu.mdi.chimera.patch.ThetaPatch;

public class ChimeraAlgorithmResult {

    // The list of ALL cells INCLUDING KISS that intersect the sphere.
    private List<Cell> intersectingCells;

    // The list of KISS cells that intersect the sphere.
    private List<Cell> kissCells;

    // The list of prepatches built from non-Kiss intersecting cells.
    private List<PrePatch> prePatches;

    // The list of theta-patches produced by the theta splice.
    private List<ThetaPatch> thetaPatches;
    
    // The final, complete patches after theta and phi slicing
    private List<Patch> patches;
    
    // The total normalized area of the prepatches, used for feedback and debugging.
    private double prePatchArea = 0.0;
    
    // The total normalized area of the theta patches, used for feedback and debugging.
    private double thetaPatchArea = 0.0;
    
    // The total normalized area of the final patches, used for feedback and debugging.
    private double patchArea = 0.0;

    // Private constructor to enforce use of factory method for creating empty results.
    private ChimeraAlgorithmResult() {
        intersectingCells = null;
        prePatches        = null;
    }

    /**
     * Create an empty result.
     *
     * @return an empty result
     */
    public static ChimeraAlgorithmResult empty() {
        return new ChimeraAlgorithmResult();
    }

    // -----------------------------------------------------------------------
    // Intersecting cells
    // -----------------------------------------------------------------------

    /**
     * Get the list of all intersecting cells (including Kiss).
     *
     * @return the list of intersecting cells, or {@code null} if not yet set
     */
    public List<Cell> getIntersectingCells() {
        return intersectingCells;
    }

    /**
     * Get the list of Kiss cells.
     *
     * @return the list of Kiss cells, or {@code null} if not yet set
     */
    public List<Cell> getKissCells() {
        return kissCells;
    }

    /**
     * Get the count of intersecting cells (excluding Kiss).
     *
     * @return the count, or 0 if not yet set
     */
    public int getIntersectingCellCount() {
        return intersectingCells == null ? 0 : intersectingCells.size() - getKissCellCount();
    }

    /**
     * Get the count of Kiss cells.
     *
     * @return the count, or 0 if not yet set
     */
    public int getKissCellCount() {
        return kissCells == null ? 0 : kissCells.size();
    }

    /**
     * Set the list of intersecting cells. Automatically partitions Kiss cells
     * into a separate list.
     *
     * @param intersectingCells the list of intersecting cells to set
     */
    public void setIntersectingCells(List<Cell> intersectingCells) {
        this.intersectingCells = intersectingCells;
        this.kissCells = intersectingCells.stream()
                .filter(cell -> cell.getIntersectionType() == IntersectionType.KISS)
                .toList();
    }

    // -----------------------------------------------------------------------
    // PrePatches
    // -----------------------------------------------------------------------

    /**
     * Get the list of prepatches.
     *
     * @return the list of prepatches, or {@code null} if not yet built
     */
    public List<PrePatch> getPrePatches() {
        return prePatches;
    }

    /**
     * Get the count of prepatches.
     *
     * @return the count, or 0 if not yet built
     */
    public int getPrePatchCount() {
        return prePatches == null ? 0 : prePatches.size();
    }

    /**
     * Set the list of prepatches.
     *
     * @param prePatches the prepatch list to store
     */
    public void setPrePatches(List<PrePatch> prePatches) {
		this.prePatches = prePatches;
		Collections.sort(this.prePatches);

		// compute total area for feedback and debugging
		prePatchArea = 0;
		if (prePatches != null) {

			for (PrePatch prePatch : prePatches) {
				prePatchArea += prePatch.areaEstimate();
			}
		}
	}

	// -----------------------------------------------------------------------
	// ThetaPatches
    // -----------------------------------------------------------------------

    public List<ThetaPatch> getThetaPatches() { return thetaPatches; }

    public int getThetaPatchCount() { return thetaPatches == null ? 0 : thetaPatches.size(); }

    /**
	 * Set the list of theta patches.
	 *
	 * @param thetaPatches the theta patch list to store
	 */
	public void setThetaPatches(List<ThetaPatch> thetaPatches) {
		this.thetaPatches = thetaPatches;
        Collections.sort(this.thetaPatches);

		// compute total area for feedback and debugging
		thetaPatchArea = 0;
		if (thetaPatches != null) {
			for (ThetaPatch thetaPatch : thetaPatches) {
				thetaPatchArea += thetaPatch.areaEstimate();
			}
		}
	}
	
	// -----------------------------------------------------------------------
	// Final Patches
    // -----------------------------------------------------------------------

		/** @return sorted final patches, or {@code null} before phi splicing */
		public List<Patch> getPatches() { return patches; }

		/** @return number of final patches, or zero before phi splicing */
		public int getPatchCount() { return patches == null ? 0 : patches.size(); }
		
		/**
		 * Set the list of final patches.
		 * 
		 * @param patches the final patch list to store
		 */
		public void setPatches(List<Patch> patches) {
			this.patches = patches;
			Collections.sort(this.patches);
			
			// compute total area for feedback and debugging
			patchArea = 0;
			if (patches != null) {
				for (Patch patch : patches) {
					patchArea += patch.areaEstimate();
				}
			}
		}

    // -----------------------------------------------------------------------
    // Feedback
    // -----------------------------------------------------------------------

    /**
     * Adds a summary of the result to the provided feedback list.
     *
     * @param colorString a string representing the color to be used in the
     *                    feedback (pass {@code ""} for no color tag)
     * @param feedbackList the list to which the feedback summary will be added
     */
    public void feedbackSummary(String colorString, List<String> feedbackList) {
        int total    = intersectingCells == null ? 0 : intersectingCells.size();
        int kissCount = getKissCellCount();
        int nonKiss  = total - kissCount;

        /*
         * Imported interchange files intentionally contain final patches but not
         * the intermediate intersecting-cell scan. Each available algorithm stage
         * must therefore contribute feedback independently; absence of cells must
         * not suppress the imported final-patch count and area.
         */
        if (intersectingCells != null) {
            feedbackList.add(String.format(
                    "%sintersecting cells  non-kiss: %d  kiss: %d  total: %d",
                    colorString, nonKiss, kissCount, total));
        }

        if (prePatches != null) {
            feedbackList.add(String.format(
                    "%sprepatches: %d normalized area: %.12f", colorString,
                    prePatches.size(), prePatchArea));
        }
        if (thetaPatches != null) {
            feedbackList.add(String.format(
                    "%stheta patches: %d normalized area: %.12f", colorString,
                    thetaPatches.size(), thetaPatchArea));
        }
        if (patches != null) {
			feedbackList.add(String.format(
					"%sfinal patches: %d normalized area: %.12f", colorString, patches.size(), patchArea));
		}
    }
}
