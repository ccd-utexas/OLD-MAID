# ProEMOnline
Software for online reduction and analysis of time series photometry obtained with the ProEM Camera (in SPE format with LightField)

### GUI Development

We are just beginning development of the software.  Initially, the "master" branch will contain a general layout of the GUI, the individual aspects of which will be developed in their own branches and incoporated once they reach basic functionality.  The GUI is powered by the DockArea system of PyQtGraph (http://www.pyqtgraph.org/).  The GUI will contain the indicated docks:

 - Observing Log: Here the user can fill out mandatory information about the observing run (name, filter used, etc.) in a way that will standardize the format for later automated reference.

 - Plots: To visulalize the results of the online photometry in a variety of useful ways.  These include the time series raw light curves of the target and selected comparison stars, the divided light curve of the target (optionally smoothed), the FT of the target (neglecting any bad points selected by the observer), and time series of sky brightness and seeing.  These are all highly pannable, scalable.
 
- Image Display: Show the images as they come in and are reduced, indicating the centroids used in the photometric measurements of selected stars. Locations of target/comparison stars are mouse-selected in this pane.

- Process Log: This dock shows a process log, listing all actions made by the software (data loading/saving, errors/warnings, etc.).

- Menu bar: full access to useful/advanced settings and functions.

### Behind the GUI

No online reduction would be complete without reduction.  We need software to do the following behind the scenes:

- Automatic dark subtraction and flat fielding of science frames.
- Robust, python-based aperture photometry
- Automatic saving and management of auxiliary files including records of clock adjustments and weather data.
- Robust timestamp validation.

