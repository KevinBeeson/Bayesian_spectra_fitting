
# Fitting Galah spectra using a bayesian procedure


The code in this repository is used to fit [Galah spectra](https://www.galah-survey.org/) using a bayesian fitting procedure with priors. This procedure is outlined in the paper [The GALAH survey: elemental abundances in open clusters using joint effective temperature and surface gravity photometric priors](https://academic.oup.com/mnras/article/529/3/2483/7607381?login=false). The main fitting program is contained in Payne_machine_solar_finder.py. The rest of the repositorycontains some diagnostic tools and some plotting prorgams used to make the plots in the paper.

## Capabilities and functions of Payne_machine_solar_finder.py


This program is designed to facilitate the analysis of stellar spectra from the GALAH DR4 survey. It provides a comprehensive set of features for working with observed and synthetic spectra, including data retrieval, normalization, and plotting. Below is a summary of the key features and functionalities available in this class.

### Data Retrieval and Storage
- **Local Storage of Spectra**: Copies fits files containing un-normalized spectra, uncertainties, and resolution profiles from the computing cluster, enabling easy and fast access.
- **Resolution Profile Handling**: If the requested resolution profile is unavailable or invalid, the class retrieves the profile of the closest observation in time on the same fiber. The original profile is replaced with this one if necessary.

### Spectral Data Management
- **Spectral Bands**: Users can specify which bands of the spectra to load from the fits file.
- **Initial Stellar Parameters**: The class obtains initial values for $T_{\rm{eff}}$, $\log(g)$, $v_{\rm{mic}}$, $v_{\rm{broad}}$, and 26 elemental abundances from the GALAH DR4 catalogue or reduction process. Alternatively, users can provide an `astropy` table with initial values. Default abundances are set to solar values in the absence of measurements.

### Synthetic Spectra Generation
- **Parameter Input**: Users can create synthetic spectra by inserting a Python dictionary with the desired values. Unspecified parameters default to stored class values, maintaining a clear and readable code structure.
- **Sample Rate Specification**: The sample rate for synthesized spectra can be specified in the PySME implementation.
- **Re-normalization**: The class can re-normalize spectra using generated synthetic spectra.
- **Wavelength Range Specification**: Users can modify a single line of code to specify wavelength ranges for synthesis.

### Plotting and Visualization
- **Spectra Plotting**: The class can plot both synthetic and observed spectra. Users can annotate the plot with theoretical positions of strong lines and indicate areas masked during fitting.
- **Wavelength Masking**: A wavelength mask can be created, defining areas excluded from fitting.

### Output Files
- **Radial Velocities**: A `.npy` file containing fitted radial velocities and their associated bands.
- **Fitted Parameters**: A `.npy` file with all fitted parameters (e.g., $T_{\rm{eff}}$, [Li/Fe], etc.) processed by `emcee`.
- **Sampler State**: The class can save the state of the `emcee` sampler to an HDF5 file, allowing sampling to resume from the last saved point. Two HDF5 files are produced: one from initial sampling to create the wavelength mask, and another from fitting spectra using the mask.
- **Fitted Values and Diagnostics**: A `.fits` file containing an `astropy` table with mean fitted values, their precision, and the auto-correlation time of each parameter, including whether a photometric prior was used.
- **Diagnostic Plot**: A plot comparing the synthesized spectra with mean parameters, observed spectra, residual differences, and important line locations.




