# metasurface-enabled-lensless-endoscope
MATLAB code and supporting data for metasurface-assisted two-photon lensless imaging through a multicore fiber.

# \# Metasurface-Enabled Lensless Endoscope

# 

# MATLAB code and supporting data for metasurface-assisted two-photon lensless imaging through a multicore fiber.

# 

# This repository contains MATLAB scripts for mode-overlap simulations, distal focusing analysis, and experimental core-intensity measurements. Each analysis can be run independently using the input files and helper functions provided in the repository.

# 

# \## Requirements

# 

# \- MATLAB with the plotting functions used by the scripts, including `tiledlayout`, `turbo`, `clim`, and `exportgraphics`.

# \- Image Processing Toolbox for core detection and z-scan analysis.

# \- Curve Fitting Toolbox for Gaussian fitting in the PSF analysis and focusing-efficiency simulations.

# 

# \## Getting Started

# 

# 1\. Clone or download the repository and select the repository folder as the MATLAB \*\*Current Folder\*\*.

# 2\. Keep the scripts, helper functions, and input files in their supplied locations.

# 3\. Open the desired script and click \*\*Run\*\*, or enter its name in the MATLAB Command Window.

# 

# For example:

# 

# ```matlab

# pillar\_lookup\_interpolation

# ```

# 

# This loads the supplied GaN lookup table and displays transmission and phase maps, together with cuts at a pillar height of 1.5 µm. Parameters for each analysis are defined near the beginning of its script.

# 

# \## Available Analyses

# 

# | Script                                    | Purpose                                                               | Required input                                    |

# | ----------------------------------------- | --------------------------------------------------------------------- | ------------------------------------------------- |

# | `pillar\_lookup\_interpolation.m`           | Interpolate and plot nanopillar transmission and phase                | `pillars\_GaNonSaph\_wv920\_prd400\_thk1500nm.mat`    |

# | `mode\_overlap\_vs\_axial\_distance.m`        | Calculate mode overlap versus propagation distance                    | `PhaseVector400um.mat`                            |

# | `mode\_overlap\_z0\_sweep.m`                 | Compare mode overlap for different metasurface design focal distances | `PhaseVector400um.mat`                            |

# | `mode\_overlap\_vs\_lateral\_offset.m`        | Calculate mode overlap versus lateral displacement                    | `PhaseVector400um.mat`                            |

# | `mode\_overlap\_vs\_rotation.m`              | Calculate mode overlap versus relative rotation                       | `PhaseVector400um.mat`                            |

# | `mode\_overlap\_vs\_sampling\_phase\_levels.m` | Compare spatial sampling and phase discretization                     | `PhaseVector400um.mat`                            |

# | `focusing\_efficiency\_vs\_phase\_levels.m`   | Simulate distal focusing with different numbers of phase levels       | Parameters in the script                          |

# | `focusing\_efficiency\_vs\_phase\_error.m`    | Simulate distal focusing with random inter-core phase errors          | Parameters in the script                          |

# | `psf\_fwhm\_comparison.m`                   | Compare MS and SLM focal images, profiles, and FWHM                   | MS/SLM focus images and corresponding dark frames |

# | `exp\_coupling\_eff\_slm\_ms.m`               | Calculate the distal-facet core intensity fraction                    | `DistalEnd\_MS.fig` and `DistalEnd\_SLM.fig`        |

# | `z\_scan\_analysis.m`                       | Detect lenslet foci and display experimental X-Z sections             | `data.hdf5`, to be supplied separately            |

# 

# The mode-overlap scripts use `phase\_mask\_mcf.m` and `local\_voronoi\_segmentation.m`.

# 

# The focusing simulations use `FiberClassGenerator.m`, `fermatspiral.m`, `CalculatePhaseCore.m`, and `GetStrehl\_useMax\_Generic.m`.

# 

# These helper functions are called by the main scripts and do not need to be run individually. The supplied simulations use the `Spiral` fiber configuration.

# 

# \## Experimental Focal Analysis

# 

# Run `psf\_fwhm\_comparison.m` to compare the measured MS and SLM focal images.

# 

# The matching dark frames are subtracted before analysis, and the script displays the focal images, intensity profiles, and fitted FWHM values for direct comparison.

# 

# \## Experimental Core-Intensity Fraction

# 

# Run `exp\_coupling\_eff\_slm\_ms.m`.

# 

# For the MS image, select three points on the fiber boundary when prompted. The SLM boundary is predefined as:

# 

# ```matlab

# \[590, 615, 430]

# ```

# 

# corresponding to:

# 

# ```text

# \[x centre, y centre, radius]

# ```

# 

# To repeat the same MS analysis, enter the printed circle coordinates in the first row of `fiberCircleByCase`.

# 

# The script subtracts the image-border offset, detects core regions, and integrates their intensity using circular apertures with a diameter equal to twice the equivalent FWHM diameter.

# 

# The core fraction is calculated as the integrated intensity within the detected core regions divided by the total intensity inside the fiber boundary, including the cladding.

# 

# The results table reports the detected core count and intensity fractions. The script also writes a diagnostic PNG, a MAT results file, and detected-centre CSV files to the `output/` folder. Existing files in this folder are example results from the supplied repository and are overwritten when the script is rerun.

# 

# \## Input Data

# 

# \### Phase Vector

# 

# `PhaseVector400um.mat` contains `phase\_vec`, the calibrated phase values for the 120 fiber cores in radians. These phase pistons were experimentally measured using the SLM with the co-propagating reference method.

# 

# \### Nanopillar Lookup Table

# 

# The supplied MAT file contains:

# 

# \- `R`: pillar radii

# \- `z\_max\_p`: pillar heights

# \- `Tr`: transmission

# \- `DeltaX`: phase

# 

# Lengths are expressed in metres and phase in radians. The wavelength is 920 nm and the metasurface lattice period is 400 nm.

# 

# \### Focus Images

# 

# The `Focus\_\*.mat` files contain the variable `data\_img`.

# 

# Their corresponding dark-frame files contain the variable `DarkFrame`.

# 

# \### Distal-Facet Images

# 

# `DistalEnd\_MS.fig` and `DistalEnd\_SLM.fig` contain the experimentally measured MS and SLM distal-facet camera-intensity images.

# 

# \### Axial Stack

# 

# The experimental z-scan dataset `data.hdf5` is not included in this GitHub repository because of its large file size (\~3.5 GB).

# 

# To run `z\_scan\_analysis.m`, place `data.hdf5` in the same folder as the script. The dataset is currently available from the authors upon request.

# 

# The expected dataset is:

# 

# ```text

# /frames```

# 

