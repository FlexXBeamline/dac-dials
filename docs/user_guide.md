# `dac-dials`

How to process diamond anvil cell (DAC) data collected at CHESS using DIALS.

## Introduction

The DAC experiment is different from "normal" crystallography in several ways. First, the crystal sits inside a metal gasket, and the gasket will absorb diffracted X-rays casting a "shadow" on the detector. The region of the detector that is shadowed changes as the goniometer rotates, and it also depends on the position of the crystal relative to the center of the gasket. Second, the diamond absorbs diffracted X-rays in a manner depending on the path length through the diamond. Finally, the usable rotation range of the DAC is limited by the gasket diameter and the exit window, and consequently obtaining a complete dataset can be challenging for low-symmetry space groups. It is especially difficult to achieve complete data at high resolution due to the angle-dependent shadowing effect mentioned above. Typically, datasets are assembled from multiple crystals loaded in the DAC at once.

We solve the gasket "shadow" problem by fitting a geometric model to the distribution of indexed strong reflections on each frame. The script for doing this is included in this repository.

DIALS includes a program called `dials.anvil_correction` that corrects integrated intensities (and their standard errors) for the diamond absorption. You can read more about it here: <https://dials.github.io/documentation/programs/dials_anvil_correction.html>. 

Merging multi-crystal data is done using the DIALS program `xia2.multiplex`. For a tutorial, see: <https://dials.github.io/documentation/tutorials/br_lyso_multi.html>

Please note that DAC data processing is an active area of development at CHESS, and we are still refining best practices. This repository will be updated frequently. If you find any issues with our code, please get in touch!

## Getting started

First, you'll need DIALS. 

If you are processing data on a MacCHESS computer or CLASSE compute farm node, you can activate dials as follows:
```
source /nfs/chess/sw/macchess/dials/dials-current/dials_env.sh 
```

If you are processing at home, installation instructions are here: <https://dials.github.io/installation.html>

## Organizing directories

First, create a directory for your project. If processing locally on a CHESS computer, this should be in the "aux" folder for your beamtime:

```
/nfs/chess/aux/cycles/<run-cycle>/id7b2/<investigator>/<date>/
```

where the names in brackets should be replaced with correct values, i.e. `<run-cycle>` would be `2023-3` if data were collected in the Fall run of 2023. To keep the output files of DIALS straight, I recommend a directory structure like this:

```
.
├── dials_processing_notes.md
├── fit_gasket_mask.py
├── sample1
│   ├── ambient_crystal1
│   ├── ambient_crystal2
│   ├── ambient_crystal3
│   ├── 2kbar_crystal1
│   ├── 2kbar_crystal2
└── sample2
    ├── ambient_crystal1
    ├── ambient_crystal2
    ├── 2kbar_crystal1
    ├── 2kbar_crystal2
    ├── 3kbar_crystal1
    └── 3kbar_crystal2
```

The DIALS extension `fit_gasket_mask.py` should be downloaded from this repository and placed in the top level directory (or git clone the repository and use that as your base directory data for processing).

I recommend keeping a log of the commands you run in a text file, as well as your notes about decisions made, so that you can easily reproduce processing later: in the example above, it's called `dials_processing_notes.md`.

Here, a "sample" refers to a particular kind of crystal (for instance, a space group or particular ligand soak), and the subfolders refer to single datasets collected on that crystal, possibly under different conditions such as pressure. For the purposes of this tutorial, I'll assume the above directory structure.

## Processing each dataset

A typical DIALS workflow looks like this:

```
cd <dataset-directory>
dials.import <path-to-master-h5-file> [options]
dials.find_spots imported.expt
dials.index imported.expt strong.refl [options]
dials.refine indexed.expt indexed.refl
dials.integrate refined.expt refined.refl
dials.python ../../fit_gasket_mask.py integrated.expt integrated.refl
dials.anvil_correction integrated.expt masked.refl [options]
```

Each step and various options are described below.

### 1. Importing

The command `dials.import` will extract geometry information stored in the `_master.h5` image and create a text file, `imported.expt`. For the DAC experiments at ID7B2 performed prior to November 2023, it's necessary to override the detector distance, rotation angle, and spindle axis information. 

```
cd sample1/ambient_crystal1
dials.import <path-to-master-h5-file> distance=<distance> axis=0,1,0 geometry.scan.oscillation=<start-angle>,<angle-increment> 
```

On the CHESS network, the path will be something like `/nfs/chess/raw/<run-cycle>/id7b2/<investigator>/<date>/<...>/<prefix>_master.h5`.

The scan axis should always be specified as `0,1,0` (the DAC goniometer axis is vertical). 

The detector distance and oscillation ranges are stored in the image header, but may not be correct depending on when the data were collected. The distance is set in mm, for instance `distance=370`. The oscillation range is given by two numbers, the first is the angle at the start of the scan, and the second is the angle increment per image. For instance, `geometry.scan.oscillation=-25,0.5` means a scan starting at -25 degrees with 0.5 degrees per frame.

If import is successful, it will generate a file `imported.expt`. You can view information about the experiment like this:

```
dials.show imported.expt
```

This is helpful for making sure the geometry is correct. You can also inspect the images using:

```
dials.image_viewer imported.expt
```

### 2. Spotfinding, indexing, and integration

First, find all of the intense spots

```
dials.find_spots imported.expt
```

Then, index the reflections (find a lattice consistent with the data) as follows:

```
dials.index imported.expt strong.refl
```

I find that for multi-crystal data collected using the DAC, it's important to specify the space group here so that the cell angle restraints can be used in geometry refinement:

```
dials.index imported.expt strong.refl space_group=<sg>
```

For instance, `space_group=P212121`. If indexing fails, this is probably a bad dataset and not worth processing further. Or the initial geometry might be incorrect (see "Importing" section, above).

Finally, refine the geometry using the indexed data:

```
dials.refine indexed.expt indexed.refl
```

### 3. Integration

The default parameters for integration usually work just fine.

```
dials.integrate refined.expt refined.refl
```

If you know the resolution cutoff, you can specify it here to avoid integrating non-existent reflections:

```
dials.integrate refined.expt refined.refl d_min=<resolution>
```

For instance, to cut off at 2.1 Å, use `d_min=2.1`. This is optional (a resolution cutoff will be imposed later at the scaling step).

### 4. Masking the gasket shadow

Next, we will fit a model for the gasket aperture. The included script is run using `dials.python`.

```
dials.python ../../fit_gasket_mask.py integrated.expt integrated.refl
```

The `../../fit_gasket_mask.py` assumes that the script is two directory levels above the current directory, which should always be the case if using the directory structure above.

The script saves a file `masked.refl` where the integrated reflections that are occluded have been removed.

There are additional advanced options that can be specified if needed. To see the list of options, type:

```
dials.python ../fit_gasket_mask.py -c -a1
```

The script assumes that the diamond windows are normal to the beam when the rotation angle is zero. If not, the correct angle should be set using the option `phi0=<angle-in-degrees>`.

### 5. Correcting diamond absorption

To run the anvil correction script included with DIALS, you have to specify the diamond thickness and the vector normal to the plane of the diamond window when the rotation angle is zero. 

```
dials.anvil_correction integrated.expt masked.refl thickness=0.8 normal=0,0,1
```

The script saves a file `corrected.refl`. If the anvil correction is working, the R-factors should improve. I recommend processing `corrected.refl` and `masked.refl` separately to check whether the anvil correction was helpful or not.

## Merging multiple datasets

After you have performed the above steps for several isomorphous crystals, you can scale and merge them collectively. In this example, let's assume we're going to merge the sample1 data collected at 2kbar. The procedure is to go to the `sample1` directory, and then make a directory for processing and `cd` into it, and run `xia2.multiplex`. This program automatically merges multi-crystal data and outputs a lot of statistics and reports, which you should scrutinize. The merged and unmerged intensities will be output as an mtz file for further processing.

First, try merging without the anvil correction:

```
mkdir 2kbar_multiplex
cd 2kbar_multiplex
xia2.multiplex ../2kbar_*/masked.refl ../2kbar_*/integrated.expt
```

Then, try merging with the anvil correction:

```
mkdir 2kbar_multiplex_corrected
cd 2kbar_multiplex_corrected
xia2.multiplex ../2kbar_*/corrected.refl ../2kbar_*/integrated.expt
```

Compare the processing statistics for both the corrected and uncorrected data. A good place to start is the table called "Summary of merging statistics" in xia2.multiplex.log. You can also find more information in the various html reports.