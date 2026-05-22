package edu.cnu.mdi.chimera.view;

import java.awt.Color;

import edu.cnu.mdi.component.CommonBorder;
import edu.cnu.mdi.component.OptionPanel;
import edu.cnu.mdi.ui.fonts.Fonts;
import edu.cnu.mdi.util.UnicodeUtils;

/**
 * Option panel for the Chimera view. This is a simple wrapper around {@link OptionPanel}
 * that provides the default options for the Chimera view.
 */
@SuppressWarnings("serial")
public class ChimeraOptionPanel extends OptionPanel  {
    public ChimeraOptionPanel(OptionPanel.OptionPanelListener listener) {
		 super(listener, 3, Fonts.plainFontDelta(-2),
				 Color.black, null, AS_ORDERED, displayDefaults());
		 setBorder(new CommonBorder("Display Options"));
	}
        
    // default values for the options table
    private static String[][] displayDefaults(){
    	String[][] defStr = {
    			{"Spherical Grid", "true", "true"}, 
       			{"Prepatches", "true", "false"},
    			{"Monte Carlo", "true", "true"},
 			    {UnicodeUtils.SMALL_THETA + " Patches", "true", "false"},
    			{"Kiss Markers", "true", "false"},
    			{"Patches", "true", "true"}};
     	return defStr;
    }
    
    /**
     * Returns whether the "Spherical Grid" option is selected.
     * @return true if the "Spherical Grid" option is selected, false otherwise
     */
    public boolean showSphericalGrid() {
		return isOptionSelected("Spherical Grid");
	}
    
    /**
	 * Returns whether the "Prepatches" option is selected.
	 * @return true if the "Prepatches" option is selected, false otherwise
	 */
    public boolean showPrepatches() {
		return isOptionSelected("Prepatches");
	}
	
    /**
	 * Returns whether the "Monte Carlo" option is selected.
	 * @return true if the "Monte Carlo" option is selected, false otherwise
	 */
	public boolean showMonteCarlo() {
		return isOptionSelected("Monte Carlo");
	}
	
	/**
	 * Returns whether the "Theta Patches" option is selected.
	 * @return true if the "Theta Patches" option is selected, false otherwise
	 */
	public boolean showThetaPatches() {
		return isOptionSelected(UnicodeUtils.SMALL_THETA + " Patches");
	}
	
	/**
	 * Returns whether the "Kiss Markers" option is selected.
	 * @return true if the "Kiss Markers" option is selected, false otherwise
	 */
	public boolean showKissMarkers() {
		return isOptionSelected("Kiss Markers");
	}
	
	/**
	 * Returns whether the "Patches" option is selected.
	 * @return true if the "Patches" option is selected, false otherwise
	 */
	public boolean showPatches() {
		return isOptionSelected("Patches");
	}
    
}
