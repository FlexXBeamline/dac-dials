# `dac-dials`

[![DOI](https://zenodo.org/badge/722199920.svg)](https://doi.org/10.5281/zenodo.17058003)

A [DIALS](https://dials.github.io/) extension for processing high pressure macromolecular crystallography data collected using a diamond anvil cell (DAC). 

In this method, crystals are placed in a metal gasket sandwiched between two diamonds, and the (hydrostatic) pressure is controlled by the force pushing the diamonds together. The X-ray beam passes through the diamonds and gasket. Because the gasket is relatively thick, it can block high-angle diffraction, casting a "shadow" on the detector. The shadowed region changes during data collection as the DAC is rotated, and it depends on the three-dimensional location of the crystal within the gasket. To properly account for this, we refine a geometric model of the DAC that best accounts for the Bragg peaks we observe.

DIALS already includes a program called `dials.anvil_correction` that corrects for the diamond absorption, described here: <https://dials.github.io/documentation/programs/dials_anvil_correction.html>. 

`Dac-dials` adds two additional command-line programs:

- `dac.fit_gasket_mask` - Refine the gasket geometry to account for occlusion of diffracted rays during rotation.
- `dac.trim_image_range` - Detect and remove bad frame ranges from the beginning and the end of data collection.

With these programs, it is possible to process DAC datasets in a mostly automated fashion. For a complete data processing walk-through, see [docs/user_guide.md](docs/user_guide.md).

## Version history

0.1.0 -- first versioned release

- pip-installable into DIALS conda environment
- BSD-2 license
- add `dac.trim_image_range` utility
- more robust and performant gasket-fitting algorithm `dac.fit_gasket_mask` (old version still available as `dac.fit_gasket_mask_v0`).

0.0.1 -- initial script-based version

- download and run stand-alone script using `dials.python fit_gasket_mask.py`

## Installation

We recommend using creating a stand-alone conda environment with DIALS and dac-dials using micromamba (or equivalent conda-like package manager)

```bash
micromamba create -n dac-dials -c conda-forge
micromamba activate dac-dials
micromamba install python=3.10 dials xia2 pip
pip install git+https://github.com/FlexXBeamline/dac-dials.git
```

If processing data collected using the Eiger2 X 16M at CHESS beamline 7b2, you should also install our detector format plugin:

```bash
pip install git+https://github.com/FlexXBeamline/dials-extensions.git
```

### Usage

For a typical workflow, you would first process as normal through integration:

```bash
dials.import <path-to-master-h5-file> [options]
dials.find_spots imported.expt
dials.index imported.expt strong.refl [options]
dials.refine indexed.expt indexed.refl
dials.integrate refined.expt refined.refl
```

Then, the DAC-specific programs are run to mask bad peaks and correct diamond absorption 

```bash
dac.fit_gasket_mask integrated.expt integrated.refl
dac.trim_image_range integrated.expt masked.refl
dials.anvil_correction integrated.expt trimmed.refl [options]
```

Finally, the data are scaled and merged using `dials.scale` or `xia2.multiplex` (for multi-crystal workflows).

### Documentation

See [docs/user_guide.md](docs/user_guide.md).

### Authors/Credits

The software was written by Steve Meisburger at the Cornell High Energy Synchrotron Source (CHESS). The work was funded by the National Institutes of Health (NIH) and the National Science Foundation (NSF). If you have questions, or find a bug, please get in touch!