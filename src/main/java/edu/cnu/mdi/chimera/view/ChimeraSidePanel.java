package edu.cnu.mdi.chimera.view;

import javax.swing.JPanel;

import edu.cnu.mdi.component.CommonBorder;
import edu.cnu.mdi.component.OptionPanel;
import edu.cnu.mdi.ui.colors.ColorMapSelectorPanel;
import edu.cnu.mdi.ui.colors.ScientificColorMap;

@SuppressWarnings("serial")
public class ChimeraSidePanel extends JPanel {
	
	/** the option panel for the chimera algorithm */
	public final ChimeraOptionPanel optionPanel;

	public ChimeraSidePanel(OptionPanel.OptionPanelListener optionListener,
			ScientificColorMap defaultMap,
			ColorMapSelectorPanel.ColorMapChangeListener colorMapListener) {
		// give a vetical layout to stack the option and colormap panels
		setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.Y_AXIS));
		
		optionPanel = new ChimeraOptionPanel(optionListener);
		add(optionPanel);
		
		ColorMapSelectorPanel colorMapPanel = new ColorMapSelectorPanel(colorMapListener, defaultMap);
		colorMapPanel.setBorder(new CommonBorder("Monte Carlo Color Map"));
		add(colorMapPanel);
	}
	
	
}
